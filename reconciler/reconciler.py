#!/Applications/Anaconda3/anaconda/bin/python3
# Author: Manzini Stefano; stefano.manzini@gmail.com
# coding=utf-8

__version__ = "0.4.7.160518"

"""This program:

- is an engine that reads and parses through two dataframes (of defined
    structure), and pulls values from specific locations, defined as diverse
    experimental conditions, iteratively.
    
- performs calculations on pulled values
    
- logs results into a dataframe, with relevant informations about the analyses
    being performed and all the values related to those analyses.
    
- ADVANCED! Reconciler Advanced uses permutations instead of averages to run
    Pearson's correlation.

- STRICT! Of the full mir x lip dataset, this generates LOADS of results
    even with low pvalues.
    Strict introduces additional checks other than Pearson's test pvalue
    to accept or reject correlations.

Changelog
=========
- 16/05/2018 Compiled some regexes for fragments of performance gains

- 23/11/2017 Added @functools.rlu_cache and cleaned up functions

- 25/11/2016 Reconciler drops NaN-containing rows after producing the
    dataframe (those can arise from very rare cases when miRNAs all have
    the same expression levels in all genotypes, causing some tests to
    pass but without any biological significance).

- 22/11/2016 now added the option of running statfriendly without
    having to merge the diets

- 01/04/2016 now handled through argparse module.

- 12/04/2016 now fully customizable through command line

- 27/05/2016 implementing new runmodes to comply with statisticians
    suggestions. Now "statfriendly" does a lot of stuff to generate
    better datapoints.


dependencies
============

> This program requires numpy, scipy and related scientific packages.
All of it is packaged into Anaconda - https://www.continuum.io/downloads

Also, it's better/required that you additionally install:

> Seaborn - https://stanford.edu/~mwaskom/software/seaborn/

pip install seaborn

> Plotly - https://plot.ly/python/getting-started/

pip install plotly

> Alternate regex module - https://pypi.python.org/pypi/regex

Download it, then cd into the dir and:
<your python> setup.py install

"""

# ==================================================
# ============     arguments parsing    ============
# ==================================================

import argparse
import sys
parser = argparse.ArgumentParser(
    description = (
        "Reconciler compares miRNomics and lipidomics datasets of "
        "Cariplo 2.1 Project."
    )
)

parser.add_argument(
    "mirna_dataframe",
    help = "Filename of miRNA DataFrame."
)

parser.add_argument(
    "lipid_dataframe",
    help = "Filename of lipid DataFrame."
)
        
parser.add_argument(
    "-p",
    "--pval",
    help = ("Sets the threshold for the p value. Higher values "
             " are rejected. Defaults to 0.01."),
    dest = "pear_pval",
    metavar = "",
)

parser.add_argument(
    "-c", 
    "--corr",
    help = """Sets the threshold for the correlation: 0 = no correlation, 
            1 = very strong correlation. 
            Lower values are rejected. Defaults to 0.7.""",
    metavar = "",
    dest = "pear_corr_threshold"
)

parser.add_argument(
    "-s", 
    "--sep",
    help = ("Sets the preferred csv separator. Defaults to <tab>."),
    metavar = "",
    dest = "outfile_separator"
)

parser.add_argument(
    "-n",
    "--nan",
    help = """Changes unavailable values in dataframes to specified
            numeric value. Defaults to 0.""",
    metavar = "",
    dest = "nan_value"
)

parser.add_argument(
    "-t",
    "--stdev-threshold",
    help = """Evaluates the difference of each biological replicate from 
            the standard deviation of all biological replicates. 
            Percentage deviations higher than specified will be rejected. 
            Defauts to 150 percent.""",
    metavar = "",
    dest = "stdev_test_threshold"
)

parser.add_argument(
    "-a",
    "--allowed-zero-lipid-averages",
    help = """Sets the number of averages in the lipid dataframe 
            that are allowed to be 0. Defaults to 1 (if more than one 
            average within the three genotypes is 0, tests are not 
            performed).""",
    metavar = "",
    dest = "allowed_zero_lipid_averages"
)

parser.add_argument(
    "-T",
    "--statistics-test",
    help = """Sets the correlation to use. Two statistics tests are 
            selectable: \"pearson\" and \"spearman\", running Pearson's 
            or Spearman's correlations respectively. Defaults to 
            spearman.""",
    metavar = "",
    dest = "statistics_test"
)

parser.add_argument(
    "-r",
    "--runmode",
    help = """Sets the way tests are performed. Three runmodes are 
            selectable: \"strict\", \"average\" or \"statfriendly\". 
            When running in 
            strict, Reconciler runs the tests on all possible 
            permutations of mir-to-lip values (genotype-wise); whereas 
            when running average Reconciler just runs the tests on 
            plain averages. In this case, sd percentage (-t argument) 
            is still evaluated - before averaging bio reps. 
            Last, statfriendly tries every 
            concievable optimization to produce unbiased mir-to-lip 
            values pairings. Defaults to statfriendly.""",
    metavar = "",
    dest = "runmode"
)

parser.add_argument(
    "-d",
    "--merge-diets",
    help = ("Only effective in statfriendly runmode. If selected, diets"
            " are merged into one comparison."
            ),
    action = "store_true"
)

parser.add_argument(
    "--leave-untidy",
    help = """By default, Reconciler tidies the produced table (removes
           spurious tests that come from very rare flat miRNA 
           quantifications), but the table is being loaded whole into
           memory. This option turns this off.""",
    metavar = "",
    dest = "leave_untidy"
)

parser.add_argument(
    "-v",
    "--version",
    action = "version",
    version = __version__
)

argv = parser.parse_args()

# =========================================================
# =========///     end of arguments parsing    ///=========
# =========================================================

# ==============================
# =====     imports     ========
# ==============================
print("Importing modules..")

from numpy import (
    mean,
    std,
    average
)
                  
from scipy.stats import (
    pearsonr, # (<pearson’s coefficient>, <2-tail-p-value)
    spearmanr
)

import regex as re
import pandas as pd
import time
from random import randint
import functools

from manzutils import (
    get_current_time,
    flush_to_file,
    stdev_test,
    gen_cond,
    mirna_ids,
    lipi_ids,
    strip_mirname,
    list_permute_all
)

# =====================================
# ==///     end of imports     ///=====
# =====================================

# ==============================
# ====     FUNCTIONS     =======
# ==============================

# Reconciler StatsFriendly Mergediets Functions
# =============================================

#=========================
def tidy_dataframe(filename,
                   sep = "\t",
                   verbose = True):
# ========================
    
    if verbose == True:
        print("Tidying {} up..".format(filename))
        print("(If this hangs, next time run with --leave-untidy, "
              "and Ctrl+C now.)\r", end = "")
    dataframe = pd.read_csv(filename, sep = sep)
    dataframe = dataframe.dropna()
    dataframe.to_csv(filename, sep = sep, index = False)
    print(" " * 78 + "\r", end = "")



# =========================
@functools.lru_cache(maxsize=None)
def inspect_lipi_id(string):
# =========================

    """It reads and parses the headers of adjusted lipidomics table.
    These look like: AO05_LDLR_1f_chow
    It returns a tuple containing translated and ordered elements, likewise:
    (sample, genotype, diet, animal, sex)
       0         1       2      3     4
    """

    betternames = {
        "LI" : "liv",
        "AO" : "aor",
        "PL" : "pla",
        "PCSK9" : "p",
        "BL6" : "b",
        "LDLR" : "l",
        "WD" : "wd",
        "chow" : "cd"
    }

    try:
        #160518 match = re.search(r"(..)_(.*)_(\d*)(.)_([^_]*)", string)
        match = match_lipi_cond.search(string)
        sample = betternames[match.group(1)] # returned 0
        genotype = betternames[match.group(2)]   # returned 1
        number = str(match.group(3)) # returned 3
        sex = match.group(4)
        diet = betternames[match.group(5)]  # returned 2
        return (sample, genotype, diet, number, sex)
    except:
        raise TypeError("Something like AO_LDLR_1f_chow expected, instead:",
            string)
        
        
# ===========================
@functools.lru_cache(maxsize=None)
def inspect_mirna_id(string):
# ===========================

    """It reads and parses the headers of adjusted mirnomics table.
    These look like: l6_ile_wd_desV2, or like b31f_aor_chow_desV2
    It returns a tuple containing translated and ordered elements, likewise:
    (organ, genotype, diet, number, sex, algorithm)
    
    If there is no sex in the name, male is assumed (or at least this is
    what happens in the Cariplo 2.1 mirna dataframe).
    """
    
    betternames = {
        "chow" : "cd",
        "wd" : "wd"
    }
    
    try:
        #160518 match = re.search(r"(.)(\d+)(.*)_(.[^_]*)_(.[^_]*)_(.[^_]*)", string)
        match = match_mirna_cond.search(string)
        
        genotype = match.group(1)   # returned 1
        number = str(match.group(2)) # returned 4
        organ = match.group(4)  #returned 0 | liv aor bra wat jej duo ile
        diet = betternames[match.group(5)]  # returned 2
        sex = match.group(3) if match.group(3) else "m"   # all mirna animals are male!
        algorithm = match.group(6)  #returned 6 | either UQ or desV2
        return (organ, genotype, diet, number, sex, algorithm)
    except:
        raise TypeError("Something like b3m_aor_chow_desV2 expected.")



# =========================
def wisely_attribute_x_to_y(
    mhs, mvs, lhs, lvs, debug = False
    ):
# =========================
    
    """each is a <list>.
    mhs = a <list> containing mirna headers, paired with
    mvs = their <values> (mirna values)
    lhs = a <list> containing lipid headers, paired with
    lvs = their <values> (lipid values)
    
    returns two <lists> of values: x and y
    """
    # We now make a two-pass screening of the lists.
    # First, we look if there is a match of the same, exact animal(s).
    # We then map the remaining ys to the remaining xs choosing randomly
    # from the y pool a sex-matched y.
    
    # these will hold values as soon as we progress in the wise choices
    mhs_pass1 = []
    mvs_pass1 = []
    lhs_pass1 = []
    lvs_pass1 = []
    
    # these lists will be returned
    x_vals = []
    y_vals = []
    
    x_heads = []
    y_heads = []
    # first pass
    
    # these will contain the indexes to be processed
    m_index = []
    l_index = []
    
#     print("\n**debug** wisely_attribute_x_to_y() has been called with values:")
#     print("mhs =", mhs)
#     print("mvs =", mvs)
#     print("lhs =", lhs)
#     print("lvs =", lvs)
    
    for mh in enumerate(mhs):
        mirna_details = inspect_mirna_id(mh[1])
        # mh[0] index, mh[1] mirna header
        
        for lh in enumerate(lhs):
            lipid_details = inspect_lipi_id(lh[1])
            # lh[0] index, lh[1] lipid header
            
            if (
                mirna_details[3] == lipid_details[3] and
                mirna_details[4] == lipid_details[4]
                ):  # if sex and number match == it's the same animal
                
                # taking notes of indices
                m_index.append(mh[0])
                l_index.append(lh[0])
                
            else:
                continue # not the same animal, see if next lipid header is
        
    # /// end of first pass.
    # at this point, if the same animal is found in mirna and lipid,
    # we have indexes. We start to fill the lists to be returned with
    # those values, and deleting values-headers pairs from input lists
    
    if debug:
        print("**debug** wisely_attribute_x_to_y() has finished pass 1")
        print("m_index = ", m_index)
        print("l_index = ", l_index)
    
    # gathering headers and values, if any
    if len(m_index) > 0:
    
        for index in m_index:            
            x_vals.append(mvs[index])
            x_heads.append(mhs[index])
    
        for index in l_index:
            y_vals.append(lvs[index])
            y_heads.append(lhs[index])
    
    # producing pass1 lists: values are kept if their headers have not
    # already been exploited
    for m_header in enumerate(mhs):
        if m_header[0] not in m_index:
            mhs_pass1.append(m_header[1])
    
    for m_value in enumerate(mvs):
        if m_value[0] not in m_index:
            mvs_pass1.append(m_value[1])
    
    for l_header in enumerate(lhs):
        if l_header[0] not in l_index:
            lhs_pass1.append(l_header[1])
    
    for l_value in enumerate(lvs):
        if l_value[0] not in l_index:
            lvs_pass1.append(l_value[1])
            
    # second pass
    
    # now we try and match (any) leftover x with a randomly chosen
    # leftover y. In the worst case, we still have 3 mirna values
    # (all males) to pair with three randomly chosen male values.
    # there should be n = 3 males and n = 3 females, but in case
    # there's just two I'll have the program pull a female values
    # in case
    
    # clearing indexes
    m_index = []
    l_index = []        
    
    # this stores the already taken indexes
    already_taken = []
    if len(mhs_pass1) > 0:
        
        already_taken = [] # stores the already probed indexes
        
        for mh in enumerate(mhs_pass1):
            
            mirna_details = inspect_mirna_id(mh[1])
            # mh[0] index, mh[1] mirna header
            
            # we DON'T cycle now through lipid headers, we pick
            # then random!
            success = False
            count = 0
            already_chosen = []
            while (count < len(lvs_pass1) and success == False):
                
                pick = randint(0, (len(lvs_pass1) - 1))
                # pick is the y index now
                # the x index is mh[0]
                
                if pick in already_chosen:
                    continue
                else:
                    already_chosen.append(pick)
                
                    lipid_details = inspect_lipi_id(lhs_pass1[pick])
                    # lhs_pass1[pick]: lipid header, picked index
                    
                    if mirna_details[4] == lipid_details[4]: # sex match
                    
                        if pick not in already_taken:
                        
                            # taking notes of indices
                            m_index.append(mh[0])
                            l_index.append(pick)
                            already_taken.append(pick)
                            success = True
                    
                    count += 1
                    
                    # count is updated even if the sex is not matched.
                    # this way, after having tried all possible tries
                    # if there is no sex match, we can trigger the else
                    # statement where we will give the remaining x(es)
                    # remaining lipid values at random
            else:
                while success == False:
                    
                    # any y is going to work now, brute force pulling
                    # values out of y list
                    
                    pick = randint(0, (len(lvs) - 1))
                    if pick not in already_taken:
                        
                        # taking notes of indices
                        m_index.append(mh[0])
                        l_index.append(pick)
                        already_taken.append(pick)
                        success = True
                    
                continue # proceed to next miRNA
                        
    else:   # all the mirna-animals were already matched (UNLIKELY!)
            # at pass-1
        pass
    
    # /// end of the second pass
    
    if debug:
        print("**debug** wisely_attribute_x_to_y() has finished pass 2")
        print("m_index = ", m_index)
        print("l_index = ", l_index)
        print("mhs_pass1 = ", mhs_pass1)
        print("mvs_pass1 = ", mvs_pass1)
        print("lhs_pass1 = ", lhs_pass1)
        print("lvs_pass1 = ", lvs_pass1)
            
    if len(m_index) > 0:
    
        for index in m_index:            
            x_vals.append(mvs_pass1[index])
            x_heads.append(mhs_pass1[index])
    
        for index in l_index:
            y_vals.append(lvs_pass1[index])
            y_heads.append(lhs_pass1[index])
    
    if debug:
        print("**debug** wisely_attribute_x_to_y() about to return:")
        print("x_vals = ", x_vals)
        print("y_vals = ", y_vals)
        print("x_heads = ", x_heads)
        print("y_heads = ", y_heads)
    
    return x_vals, y_vals, x_heads, y_heads


# ========================================================================
def do_magic(e,         # from extract_values()
             algo,      # from the calling code ("desV2", "UQ")
             sample):   # from the calling code ("pla", "aor", "liv")
# ========================================================================
    
    """This function deals with Reconciler --runmode statfriendly
    modes (-d and without).
    
    it does the magic we agreed to with Mika and the
    statisticians to return two lists of the x and y values that then
    are subjected to the statistical testing.
    
    This function is almost a wrapper for the underlying
    wisely_attribute_x_to_y() function.
    The <e> dict passed to the function relates to just one diet,
    so the diet is not to be specified here.
    
    it wants: <dict>(s) generated by extract_values()
    it returns: ([x1, ... xn], [y1, ... yn], [x headers..], [y headers..])
    
    The returned values are pulled from all genotypes for that particular
    mirna, algorithm, lipid, sample.
    
    example:
    x, y, x_headers, y_headers = do_magic(e, algo = "UQ", sample = "aor")
    
    """
    
    x = []  # x (mirna) values
    y = []  # y (lipid) values
    xh = [] # x (mirna) headers
    yh = [] # y (lipid) headers
    for genotype in ("b", "p", "l"):
        
        # "algo" and "sample" are mandatorily given to do_magic()
        # so we know who they are
        
        mir_value_handle = ("mirna_expression_levels_" +
                            genotype + "_" +
                            algo)
        
        mir_header_handle = ("mirna_headers_" +
                             genotype + "_" +
                             algo)
        
        lip_value_handle = ("lipid_amount_" +
                            genotype + "_" +
                            sample)
        
        lip_header_handle = ("lipid_headers_" +
                             genotype + "_" +
                             sample)
                             
        # HERE, we have four lists containing, for:
        # one genotype, sample and algorithm
        # - mirna and lipid values 
        # - mirna and lipid headers
        
#             print("\n**debug** do_magic() about to process:")
#             print("mhs = ", e[mir_header_handle], sep ="")
#             print("mvs = ", e[mir_value_handle], sep ="")
#             print("lhs = ", e[lip_header_handle], sep ="")
#             print("lvs = ", e[lip_value_handle], sep ="")
        
        x_curr, y_curr, xh_curr, yh_curr = wisely_attribute_x_to_y(
            mhs = e[mir_header_handle],
            mvs = e[mir_value_handle],
            lhs = e[lip_header_handle],
            lvs = e[lip_value_handle]
            )
        
#             print("**debug** do_magic() genotype", genotype)
#             print("**debug** do_magic() x_curr", x_curr)
#             print("**debug** do_magic() y_curr", y_curr)
        
        x = x + x_curr
        xh = xh + xh_curr
        
        y = y + y_curr
        yh = yh + yh_curr
        
#             print("**debug** do_magic() len(y)", len(y))
#             print("**debug** do_magic() len(x)", len(x))

    
#         print("\n**debug**#################################")
#         print("**debug** do_magic() final len(x)", len(x))
#         print("**debug** do_magic() final len(y)", len(y))
#         print("**debug**#################################\n")
    
    return x, y, xh, yh  # when all genotypes have been processed

# ======///                                                        ///   ======
# ======///   END OF Reconciler StatsFriendly Mergediets Functions ///   ======

# ==================================================
def extract_values(df_mirna,
                   df_lipi,
                   diet,
                   organ,
                   get_headers = False,
                   do_strip_mirname = True  # different name than function -.-'
                   ):
# ==================================================
    
    """ It digs into mirnomics and lipidomics dataframes, and yields, at
    every call, a dictionary that contains the values you can spot in the
    code below.
    Of note, this search is restricted to a particular organ (mirnomics) and
    dietary condition, so additional code in the outer scope is required to
    do the job for all combinations.
    This allows for some flexibility in the overall analysis.
    
    This is miRNA-dataframe-centric. That is, it considers as "conditions"
    those found in the miRNOmic dataset (organs, etc) while pulling all info
    for the lipid being considered in the lipidomic dataset.
    """

    if get_headers:    
        extracted_values = {
            "mirna_expression_levels_b_desV2" : [],
            "mirna_expression_levels_p_desV2" : [],
            "mirna_expression_levels_l_desV2" : [],
            "mirna_expression_levels_b_UQ" : [],
            "mirna_expression_levels_p_UQ" : [],
            "mirna_expression_levels_l_UQ" : [],
            "lipid_amount_b_pla" : [],
            "lipid_amount_p_pla" : [],
            "lipid_amount_l_pla" : [],
            "lipid_amount_b_aor" : [],
            "lipid_amount_p_aor" : [],
            "lipid_amount_l_aor" : [],
            "lipid_amount_b_liv" : [],
            "lipid_amount_p_liv" : [],
            "lipid_amount_l_liv" : [],
            "curr_mirna_ID" : "",
            "curr_lipid_ID" : "",
            "curr_diet" : "",
            "curr_organ" : "",
            "mirna_headers_b_desV2" : [],
            "mirna_headers_p_desV2" : [],
            "mirna_headers_l_desV2" : [],
            "mirna_headers_b_UQ" : [],
            "mirna_headers_p_UQ" : [],
            "mirna_headers_l_UQ" : [],
            "lipid_headers_b_pla" : [],
            "lipid_headers_p_pla" : [],
            "lipid_headers_l_pla" : [],
            "lipid_headers_b_aor" : [],
            "lipid_headers_p_aor" : [],
            "lipid_headers_l_aor" : [],
            "lipid_headers_b_liv" : [],
            "lipid_headers_p_liv" : [],
            "lipid_headers_l_liv" : []           
        }
    else:
        extracted_values = {
            "mirna_expression_levels_b_desV2" : [],
            "mirna_expression_levels_p_desV2" : [],
            "mirna_expression_levels_l_desV2" : [],
            "mirna_expression_levels_b_UQ" : [],
            "mirna_expression_levels_p_UQ" : [],
            "mirna_expression_levels_l_UQ" : [],
            "lipid_amount_b_pla" : [],
            "lipid_amount_p_pla" : [],
            "lipid_amount_l_pla" : [],
            "lipid_amount_b_aor" : [],
            "lipid_amount_p_aor" : [],
            "lipid_amount_l_aor" : [],
            "lipid_amount_b_liv" : [],
            "lipid_amount_p_liv" : [],
            "lipid_amount_l_liv" : [],
            "curr_mirna_ID" : "",
            "curr_lipid_ID" : "",
            "curr_diet" : "",
            "curr_organ" : ""
        }
      
    if diet == None:
        raise TypeError("No diet selected. Select either cd or wd.")
    elif organ == None:
        raise TypeError("No organ selected. "
                "Select either aor, liv, bra, wat, duo, jej, ile.")
    
    
    c = (organ, None, diet)    # current c (condition) parameters // not genotype - yet!
    # written like this to match indexes from mirna_ids() used later on
    
    extracted_values["curr_diet"] = diet
    extracted_values["curr_organ"] = organ
    
    mh_hits = []    # container for mirna headers
    for mh in df_mirna.columns[1:]:
        mh_ids = mirna_ids(mh)
        
        # we now have:
        # c = (o, None, d)
        # mh_ids = (organ, genotype, diet, animal, algorithm)
        
        # We're going to pick only those that match the current condition c:
        
        if c[0] == mh_ids[0] and c[2] == mh_ids[2]: # and c[1] == mh_ids[1]
            mh_hits.append(mh)
        
        # we now have something like:
        # mh_hits = ['b6_ile_wd_desV2', 'b7_ile_wd_desV2', ... ]
    
    #print("**debug** List of header hits:", mh_hits)
    # We now have all mirna headers that match current condition
    # (that is, organ and diet).
    # Now we need to sort the corresponding values into the dictionary

    for i in range(df_mirna.shape[0]):  # doing row by row (miRNA by miRNA)
        #print("**debug** : Current mirna slice coordinates:", i, i + 1)
        
        # clearing up previous entries
        extracted_values["mirna_expression_levels_b_desV2"] = []
        extracted_values["mirna_expression_levels_p_desV2"] = []
        extracted_values["mirna_expression_levels_l_desV2"] = []
        extracted_values["mirna_expression_levels_b_UQ"] = []
        extracted_values["mirna_expression_levels_p_UQ"] = []
        extracted_values["mirna_expression_levels_l_UQ"] = []

        if get_headers:
            extracted_values["mirna_headers_b_desV2"] = []
            extracted_values["mirna_headers_p_desV2"] = []
            extracted_values[ "mirna_headers_l_desV2"] = []
            extracted_values["mirna_headers_b_UQ"] = []
            extracted_values["mirna_headers_p_UQ"] = []
            extracted_values["mirna_headers_l_UQ"] = []            
        
        df_mirna_slice = df_mirna[i : i + 1]  # still has headers and all!
        
        # =============================================================
        # appending current miRNA ID name
        if do_strip_mirname:
            extracted_values["curr_mirna_ID"] = strip_mirname(df_mirna_slice["mirna"][i])
        else:
            extracted_values["curr_mirna_ID"] = df_mirna_slice["mirna"][i]
        # fetching expression levels and putting them into extracted_values dict
        
        
        if get_headers:
            for mh_hit in mh_hits:
                #160518 if re.match(r".*_desV2$", mh_hit):
                if match_desv2.match(mh_hit):
                    if mirna_ids(mh_hit)[1] == "b":
                        extracted_values["mirna_expression_levels_b_desV2"].append(df_mirna_slice[mh_hit][i])
                        extracted_values["mirna_headers_b_desV2"].append(mh_hit)
                    elif mirna_ids(mh_hit)[1] == "p":
                        extracted_values["mirna_expression_levels_p_desV2"].append(df_mirna_slice[mh_hit][i])
                        extracted_values["mirna_headers_p_desV2"].append(mh_hit)
                    elif mirna_ids(mh_hit)[1] == "l":
                        extracted_values["mirna_expression_levels_l_desV2"].append(df_mirna_slice[mh_hit][i])
                        extracted_values["mirna_headers_l_desV2"].append(mh_hit)
                #160518 elif re.match(r".*_UQ$", mh_hit):
                elif match_uq.match(mh_hit):
                    if mirna_ids(mh_hit)[1] == "b":
                        extracted_values["mirna_expression_levels_b_UQ"].append(df_mirna_slice[mh_hit][i])
                        extracted_values["mirna_headers_b_UQ"].append(mh_hit)
                    elif mirna_ids(mh_hit)[1] == "p":
                        extracted_values["mirna_expression_levels_p_UQ"].append(df_mirna_slice[mh_hit][i])
                        extracted_values["mirna_headers_p_UQ"].append(mh_hit)
                    elif mirna_ids(mh_hit)[1] == "l":
                        extracted_values["mirna_expression_levels_l_UQ"].append(df_mirna_slice[mh_hit][i])
                        extracted_values["mirna_headers_l_UQ"].append(mh_hit)
        else:
            for mh_hit in mh_hits:
                #160518 if re.match(r".*_desV2$", mh_hit):
                if match_desv2.match(mh_hit):
                    if mirna_ids(mh_hit)[1] == "b":
                        extracted_values["mirna_expression_levels_b_desV2"].append(df_mirna_slice[mh_hit][i])
                    elif mirna_ids(mh_hit)[1] == "p":
                        extracted_values["mirna_expression_levels_p_desV2"].append(df_mirna_slice[mh_hit][i])
                    elif mirna_ids(mh_hit)[1] == "l":
                        extracted_values["mirna_expression_levels_l_desV2"].append(df_mirna_slice[mh_hit][i])
                if match_uq.match(mh_hit):
                #160518 elif re.match(r".*_UQ$", mh_hit):
                    if mirna_ids(mh_hit)[1] == "b":
                        extracted_values["mirna_expression_levels_b_UQ"].append(df_mirna_slice[mh_hit][i])
                    elif mirna_ids(mh_hit)[1] == "p":
                        extracted_values["mirna_expression_levels_p_UQ"].append(df_mirna_slice[mh_hit][i])
                    elif mirna_ids(mh_hit)[1] == "l":
                        extracted_values["mirna_expression_levels_l_UQ"].append(df_mirna_slice[mh_hit][i])
        # ======== done fetching expression levels for current miRNA ==========
        
        lh_hits = []    # container for lipid headers
        for lh in df_lipi.columns[1:]:
            lh_ids = lipi_ids(lh) # (sample, genotype, diet, animal)
            
            # only diet matters now. We're going to pick aorta, liver and
            # plasma lipid levels regardless of the organ we're coming from
            # the miRNA table.
            # That would be 6(values per condition)*3(genotypes)*3(samples)=
            # = 54 total values (out of 108 total headers).
            
            if c[2] == lh_ids[2]:
                lh_hits.append(lh)
        
        for j in range(df_lipi.shape[0]):   # doing now LIPID by LIPID
            #print("**debug** : Current lipi slice coordinates:", j, j + 1)
            df_lipi_slice = df_lipi[j : j + 1]
            
            # =================================================================
            # appending current LIPID ID name
            extracted_values["curr_lipid_ID"] = df_lipi_slice["LIPID_NAME"][j]
            # fetching lipid amount quantification, putting them into dict
            
            # clearing up previous entries
            
            
            
            for lh_hit in lh_hits:
            
                if get_headers:
                    # need to double-check each item; better to regex it once then bool test it
                    lh_hit_ids = lipi_ids(lh_hit)   # (sample, genotype, diet, animal)
                    if lh_hit_ids[0] == "aor":
                        if lh_hit_ids[1] == "b":
                            extracted_values["lipid_amount_b_aor"].append(df_lipi_slice[lh_hit][j])
                            extracted_values["lipid_headers_b_aor"].append(lh_hit)
                        elif lh_hit_ids[1] == "p":
                            extracted_values["lipid_amount_p_aor"].append(df_lipi_slice[lh_hit][j])
                            extracted_values["lipid_headers_p_aor"].append(lh_hit)
                        elif lh_hit_ids[1] == "l":
                            extracted_values["lipid_amount_l_aor"].append(df_lipi_slice[lh_hit][j])
                            extracted_values["lipid_headers_l_aor"].append(lh_hit)
                    elif lh_hit_ids[0] == "liv":
                        if lh_hit_ids[1] == "b":
                            extracted_values["lipid_amount_b_liv"].append(df_lipi_slice[lh_hit][j])
                            extracted_values["lipid_headers_b_liv"].append(lh_hit)
                        elif lh_hit_ids[1] == "p":
                            extracted_values["lipid_amount_p_liv"].append(df_lipi_slice[lh_hit][j])
                            extracted_values["lipid_headers_p_liv"].append(lh_hit)
                        elif lh_hit_ids[1] == "l":
                            extracted_values["lipid_amount_l_liv"].append(df_lipi_slice[lh_hit][j])
                            extracted_values["lipid_headers_l_liv"].append(lh_hit)
                    elif lh_hit_ids[0] == "pla":
                        if lh_hit_ids[1] == "b":
                            extracted_values["lipid_amount_b_pla"].append(df_lipi_slice[lh_hit][j])
                            extracted_values["lipid_headers_b_pla"].append(lh_hit)
                        elif lh_hit_ids[1] == "p":
                            extracted_values["lipid_amount_p_pla"].append(df_lipi_slice[lh_hit][j])
                            extracted_values["lipid_headers_p_pla"].append(lh_hit)
                        elif lh_hit_ids[1] == "l":
                            extracted_values["lipid_amount_l_pla"].append(df_lipi_slice[lh_hit][j])
                            extracted_values["lipid_headers_l_pla"].append(lh_hit)
                else:
                    # need to double-check each item; better to regex it once then bool test it
                    lh_hit_ids = lipi_ids(lh_hit)   # (sample, genotype, diet, animal)
                    if lh_hit_ids[0] == "aor":
                        if lh_hit_ids[1] == "b":
                            extracted_values["lipid_amount_b_aor"].append(df_lipi_slice[lh_hit][j])
                        elif lh_hit_ids[1] == "p":
                            extracted_values["lipid_amount_p_aor"].append(df_lipi_slice[lh_hit][j])
                        elif lh_hit_ids[1] == "l":
                            extracted_values["lipid_amount_l_aor"].append(df_lipi_slice[lh_hit][j])
                    elif lh_hit_ids[0] == "liv":
                        if lh_hit_ids[1] == "b":
                            extracted_values["lipid_amount_b_liv"].append(df_lipi_slice[lh_hit][j])
                        elif lh_hit_ids[1] == "p":
                            extracted_values["lipid_amount_p_liv"].append(df_lipi_slice[lh_hit][j])
                        elif lh_hit_ids[1] == "l":
                            extracted_values["lipid_amount_l_liv"].append(df_lipi_slice[lh_hit][j])
                    elif lh_hit_ids[0] == "pla":
                        if lh_hit_ids[1] == "b":
                            extracted_values["lipid_amount_b_pla"].append(df_lipi_slice[lh_hit][j])
                        elif lh_hit_ids[1] == "p":
                            extracted_values["lipid_amount_p_pla"].append(df_lipi_slice[lh_hit][j])
                        elif lh_hit_ids[1] == "l":
                            extracted_values["lipid_amount_l_pla"].append(df_lipi_slice[lh_hit][j])
            # =================================================================    

            yield extracted_values
            # clearing up previous entries
            extracted_values["lipid_amount_b_aor"] = []
            extracted_values["lipid_amount_p_aor"] = []
            extracted_values["lipid_amount_l_aor"] = []
            extracted_values["lipid_amount_b_liv"] = []
            extracted_values["lipid_amount_p_liv"] = []
            extracted_values["lipid_amount_l_liv"] = []
            extracted_values["lipid_amount_b_pla"] = []
            extracted_values["lipid_amount_p_pla"] = []
            extracted_values["lipid_amount_l_pla"] = []
            
            if get_headers:
                extracted_values["lipid_headers_b_pla"] = []
                extracted_values["lipid_headers_p_pla"] = []
                extracted_values["lipid_headers_l_pla"] = []
                extracted_values["lipid_headers_b_aor"] = []
                extracted_values["lipid_headers_p_aor"] = []
                extracted_values["lipid_headers_l_aor"] = []
                extracted_values["lipid_headers_b_liv"] = []
                extracted_values["lipid_headers_p_liv"] = []
                extracted_values["lipid_headers_l_liv"] = []                
# /* END OF extract_values() */


# =============================================
def get_averages(dict,
                 zero_test_threshold = 0,
                 allowed_zero_lipid_averages = 1):
# =============================================
    
    """
    Differently from add_averages(), get_averages() is designed after
    get_permutations() to behave similarly. Instead of making permutations
    it returns averages (and the results of the "zeroesness" test).
    
    * doc for the original get_permutations(): *
    
    This function works on dicts generated by extract_values().
    This functions stores all possible permutations of genotype-specific
    miRNA - to - lipid expression levels.
    It returns, at each iteration, lists containing all permutations
    of the combinations of values found in two lists (this has been written for
    the combined analysis of mirna and lipid databases), that are the sum
    of individual permutations of b, p, and l genotypes (in that order).
    
    
    zero_test_threshold is also used to introduce the rule "do not allow more
    than one average in the lipid dataframe to be zero".
    If you want to SKIP this test, use None.
    If another value other than 0 was supplied when translating NaNs in the
    lipid dataframe, this can also be specified.
    
    allowed_zero_lipid_averages allows the specification of how many
    averages in the lipid DataFrame are considered acceptable. It defaults
    to 1, although 0 is way better in terms of Pearson's reliability.
    """
    
         
    for algoz in ("desV2", "UQ"):
        for samplz in ("pla", "aor", "liv"):
            x_masterlist = []   # will average values for b, p, l (mirna)
            y_masterlist = []   # will average values for b, p, l (lipid)
            lip_zeroes_test = 0
            test = ""
            for genoz in ("b", "p", "l"):
                mir_handle = "mirna_expression_levels_" + genoz + "_" + algoz
                lip_handle = "lipid_amount_" + genoz + "_" + samplz
                
                x, y = average(dict[mir_handle]), average(dict[lip_handle])

                #x_masterlist = x_masterlist + x
                x_masterlist.append(x)
                #y_masterlist = y_masterlist + y
                y_masterlist.append(y)
                
                if average(y) <= zero_test_threshold:
                    lip_zeroes_test += 1
            
            if lip_zeroes_test > allowed_zero_lipid_averages:
                test = "ZERO_TEST_FAILED"
            else:
                test = "ZERO_TEST_PASSED"
            
            yield (x_masterlist,    # [0]
                   y_masterlist,    # [1]
                   algoz,   # [2]
                   samplz,  # [3]
                   test # [4]
                   )
            
            # resetting values
            x_masterlist = []
            y_masterlist = []
            lip_zeroes_test = 0
            test = ""
# /* END OF get_averages() */



# =============================================
def get_permutations(dict,
                     zero_test_threshold = 0,
                     allowed_zero_lipid_averages = 1):
# =============================================
    
    """This function works on dicts generated by extract_values().
    This functions stores all possible permutations of genotype-specific
    miRNA - to - lipid expression levels.
    It returns, at each iteration, lists containing all permutations
    of the combinations of values found in two lists (this has been written for
    the combined analysis of mirna and lipid databases), that are the sum
    of individual permutations of b, p, and l genotypes (in that order).
    
    # Input modification ***    
    Please note that * THE INPUT IS RETURNED MODIFIED * and this is the
    INTENDED behavior.
    
    The input dict is added other key:value pairs and the existing ones
    are not modified (if you're not re-running this on the dict, which
    would be pointless).
    
    Please ALSO note that the above depends - or not - on the commenting
    of "CODE FOR STORING IN THE DICT" lines. So, check it.
    # ********************** 
    
    zero_test_threshold is also used to introduce the rule "do not allow more
    than one average in the lipid dataframe to be zero".
    If you want to SKIP this test, use None.
    If another value other than 0 was supplied when translating NaNs in the
    lipid dataframe, this can also be specified.
    
    allowed_zero_lipid_averages allows the specification of how many
    averages in the lipid DataFrame are considered acceptable. It defaults
    to 1, although 0 is way better in terms of Pearson's reliability.
    """
    
         
    for algoz in ("desV2", "UQ"):
        for samplz in ("pla", "aor", "liv"):
            x_masterlist = []   # will hold values for all genotypes
            y_masterlist = []   # will hold values for all genotypes
            lip_zeroes_test = 0
            test = ""
            for genoz in ("b", "p", "l"):
                mir_handle = "mirna_expression_levels_" + genoz + "_" + algoz
                lip_handle = "lipid_amount_" + genoz + "_" + samplz
                
                x, y = list_permute_all(dict[mir_handle], dict[lip_handle])

                # *** CODE FOR STORING IN THE DICT ***
                #x_key = genoz + "_" + algoz + "_" + samplz + "_x"
                # something like b_desV2_pla_x
                #y_key = genoz + "_" + algoz + "_" + samplz + "_y"
                # something like b_desV2_pla_y
                
                # storing in the dict
                #dict[x_key] = x
                #dict[y_key] = y
                # *** END OF CODE FOR STORING IN THE DICT ***
                
                x_masterlist = x_masterlist + x
                y_masterlist = y_masterlist + y
                
                if average(y) <= zero_test_threshold:
                    lip_zeroes_test += 1
                
            # *** CODE FOR STORING IN THE DICT ***
            # at this point we have cycled through the three genotypes
            # and stored x,y values for all three genotypes. We now have
            # six keys as these (corresponding to paired lists of values):
            
            # b_desV2_pla_x
            # b_desV2_pla_y
            # p_desV2_pla_x
            # p_desV2_pla_y
            # l_desV2_pla_x
            # l_desV2_pla_y
            
            # the above values are currently not meant to be included in the
            # final data table, but FIX!: the values they come from are, and
            # I'll write a function to return (properly formatted) all the
            # values to be readily fed to plotly.
            # *** END OF CODE FOR STORING IN THE DICT ***
            
            if lip_zeroes_test > allowed_zero_lipid_averages:
                test = "ZERO_TEST_FAILED"
            else:
                test = "ZERO_TEST_PASSED"
            
            yield (x_masterlist,
                   y_masterlist,
                   algoz,
                   samplz,
                   test
                   )
            
            # resetting values
            x_masterlist = []
            y_masterlist = []
            lip_zeroes_test = 0
            test = ""
# /* END OF get_permutations() */



# =====================================================================
def statsfriendly_singlediets_log(
    d, # yielded dict --> d <-- from extract_values()
    algorithm,  # string UQ, desV2
    sample, # pla, aor, liv
    pvalue, # must be given from correlation
    stats_corr, # must be given from correlation
    all_values,
    all_headers,
    sep = "\t"
    ):
# =====================================================================
    
    """Used to log a row of results found by statistics correlation,
    in a consistent way.
    """
    
    holder = (
        d["curr_mirna_ID"] + sep
        + d["curr_organ"] + sep
        + algorithm + sep
        + d["curr_lipid_ID"] + sep
        + sample + sep
        + pvalue + sep
        + stats_corr + sep
        )
    
    holder = holder + sep.join(all_headers) + sep
    
    for v in enumerate(all_values):
        if v[0] < (len(all_values) - 1):
            holder = holder + str(v[1]) + sep
        else:
            holder = holder + str(v[1]) + "\n"
    
    return holder




# =====================================================================
def statsfriendly_mergediets_log(
    d, # yielded dict --> d <-- from extract_values()
    algorithm,  # string UQ, desV2
    sample, # pla, aor, liv
    pvalue, # must be given from correlation
    stats_corr, # must be given from correlation
    all_values,
    all_headers,
    sep = "\t"
    ):
# =====================================================================
    
    """Used to log a row of results found by statistics correlation,
    in a consistent way.
    
    d: needed to get info about current mirna, lipid, etc. With
    --merge-diets there's actually two dicts, but they hold the same info.
    """
    
    holder = (
        d["curr_mirna_ID"] + sep
        + d["curr_organ"] + sep
        + algorithm + sep
        + d["curr_lipid_ID"] + sep
        + sample + sep
        + pvalue + sep
        + stats_corr + sep
        )
    
    holder = holder + sep.join(all_headers) + sep
    
    for v in enumerate(all_values):
        if v[0] < (len(all_values) - 1):
            holder = holder + str(v[1]) + sep
        else:
            holder = holder + str(v[1]) + "\n"
    
    return holder



# =====================================================================
def stat_log(d, # yielded dict --> d <-- from extract_values()
             algorithm,  # string UQ, desV2
             sample, # pla, aor, liv
             mir_expr_aves,
             mir_expr_stds,
             lip_expr_aves,
             lip_expr_stds,
             sep = "\t"
             ):
# =====================================================================
    
    """Used to log a row of results found by statistics correlation,
    in a consistent way.
    """
    
    holder = (
        d["curr_mirna_ID"] + sep
        + d["curr_organ"] + sep
        + d["curr_diet"] + sep + algorithm + sep
        + d["curr_lipid_ID"] + sep + sample + sep
        + str(d["pear_pval"]) + sep + str(d["pear_corr"]) + sep
        + str(mir_expr_aves[0]) + sep
        + str(mir_expr_aves[1]) + sep
        + str(mir_expr_aves[2]) + sep
        + str(lip_expr_aves[0]) + sep 
        + str(lip_expr_aves[1]) + sep
        + str(lip_expr_aves[2]) + sep
        + str(mir_expr_stds[0]) + sep
        + str(mir_expr_stds[1]) + sep
        + str(mir_expr_stds[2]) + sep
        + str(lip_expr_stds[0]) + sep 
        + str(lip_expr_stds[1]) + sep
        + str(lip_expr_stds[2]) + sep
        )
        
    # PLEASE NOTE! The following code EXPECTS to log three individual numbers
    # for the miRNAs and six individual numbers for the lipids, per genotype.
    # If there are MORE or LESS, this *will* result in a skewed dataframe, even
    # though the list joining, as it is written, will add (if there are less
    # values) the right number of separators nonetheless.
    
    m_keys = [("mirna_expression_levels_" + x + "_" + algorithm) for x in ("b", "p", "l")]
    l_keys = [("lipid_amount_" + x + "_" + sample) for x in ("b", "p", "l")]
    
    for m in m_keys:
        holder = holder + sep.join(map(str, d[m])) + sep
    
    for l in enumerate(l_keys):
        if l[0] < 2:
            holder = holder + sep.join(map(str, d[l[1]])) + sep
        else:
            holder = holder + sep.join(map(str, d[l[1]])) + "\n"
 
    return holder



# ==============================
# ==///  END OF FUNCTIONS  ///==
# ==============================

# ==============================
# =====     GLOBALS     ========
# ==============================

runmode = ""
if argv.runmode is None:
    runmode = "statfriendly"
elif argv.runmode == "strict":
    runmode = "strict"
elif argv.runmode == "average":
    runmode = "average"
elif argv.runmode == "statfriendly":
    runmode = "statfriendly"
else:
    print("Can't tell what \"", argv.runmode, "\" runmode is.", sep = "")
    sys.exit()

statistics_test = ""
if argv.statistics_test is None:
    statistics_test = "spearman"
elif argv.statistics_test == "pearson":
    statistics_test = "pearson"
elif argv.statistics_test == "spearman":
    statistics_test = "spearman"
else:
    print("Can't tell what \"",
          argv.statistics_test,
          "\" statistics test is.",
          sep = "")
    sys.exit()

mirna_dataframe = argv.mirna_dataframe
lipid_dataframe = argv.lipid_dataframe

if argv.outfile_separator is None:
    outfile_separator = "\t"
else:
    outfile_separator = argv.outfile_separator

# found stuff will be logged here in human readable form.. probably
logfile = ""    

# value to change NaN values in lipi dataframe to. Use a NUMBER! Not string.
if argv.nan_value is None:
    nan_value = 0
else:
    nan_value = float(argv.nan_value)

# allows to specify desired p-value threshold for correlation
if argv.pear_pval is None:
    pear_pval = 0.01
else:
    pear_pval = float(argv.pear_pval)

# and this is for drawing pictures (if it is being done)
pear_draw_pval = 0.0025

# allows to specify desired correlation threshold (0 - 1)
if argv.pear_corr_threshold is None:
    pear_corr_threshold = 0.7
else:
    pear_corr_threshold = float(argv.pear_corr_threshold)

# acceptable max % difference between all averages and their stdevs (like 70)
if argv.stdev_test_threshold is None:
    stdev_test_threshold = 150
else:
    stdev_test_threshold = float(argv.stdev_test_threshold)

# how many averages are allwowed to be 0 in the lipid dataframe?
if argv.allowed_zero_lipid_averages is None:
    allowed_zero_lipid_averages = 1
else:
    allowed_zero_lipid_averages = int(argv.allowed_zero_lipid_averages)
    
# ===============     defining output filename      ===============

if mirna_dataframe.endswith(".csv"):
    mirna_filename_prefix = mirna_dataframe.replace(".csv", "")
elif mirna_dataframe.endswith(".txt"):
    mirna_filename_prefix = mirna_dataframe.replace(".txt", "")
else:
    mirna_filename_prefix = mirna_dataframe

if lipid_dataframe.endswith(".csv"):
    lipi_filename_prefix = lipid_dataframe.replace(".csv", "")
elif lipid_dataframe.endswith(".txt"):
    lipi_filename_prefix = lipid_dataframe.replace(".txt", "")
else:
    lipi_filename_prefix = lipid_dataframe


outfile_pearson_filename = (
        runmode + "_")
        
if argv.merge_diets:
    outfile_pearson_filename = outfile_pearson_filename + "mergediets_"

outfile_pearson_filename = outfile_pearson_filename + (
        statistics_test + "_"
        + mirna_filename_prefix + "_"
        + lipi_filename_prefix + "_"
        + str(pear_pval) + "_"
        + str(pear_corr_threshold)
        + "_SD_" + str(stdev_test_threshold) + "p"
        + "_allowed0aves_" + str(allowed_zero_lipid_averages)
        + ".csv"
        )

# ============///     end of output filename      ///============

# =====================================
# ==///     END OF GLOBALS     ///=====
# =====================================


# =============================================================================
# =========            RECONCILER-SPECIFIC RUN MODES            ===============
# =============================================================================

# ====================================
def execute_statfriendly_singlediets(
    output_string,
    statistics_test = "pearson",
    debug = False,
    logfile = False
    ):
# ====================================
    
    if statistics_test == "pearson":
        print("\nEvaluating Pearson's correlations with these settings:")
    elif statistics_test == "spearman":
        print("\nEvaluating Spearman's correlations with these settings:")
    else:
        raise TypeError("Something went wrong with statistics tests checks.")

    print("Pval threshold: ", pear_pval,
          "\nPcorr threshold: ", pear_corr_threshold,
          "\nMax allowed standard deviation % from average: ",
          stdev_test_threshold, "%",
          "\nMax allowed zeroes in lipid DataFrame: ",
          allowed_zero_lipid_averages,
          sep = ""
          )
    
    max_hits = int(df_mirna.shape[0] * df_lipi.shape[0] * 2 * 7)
    # /2 because we're joining the diets
    counter = 0 # needed for the progress counter
    found = 0   # needed for the progress counter
    for organ in ("aor", "liv", "bra", "duo", "jej", "ile", "wat"):
        counter += 1    # needed for the progress counter
        # beginning to yield values cycling through miRNAs and lipids - 
        
        # FIXFIXFIX
        ### NOW!! here we need to differentiate cd and wd.
        # trying now to *just* write the stuff for cd,
        # then try also to pack it so that it repeats
        # for the other diet.
        # each test will remember how it's done,
        # so we *then* sort cd and wd.
        # REPEAT an OLD test to remember how 
        #the table headers look like !!!!!!!!!
        
        for diet in ("cd", "wd"):

            for e in extract_values(
                                    df_mirna,
                                    df_lipi,
                                    diet = diet,
                                    organ = organ,
                                    get_headers = True,
                                    do_strip_mirname = False
                                    ):
            
                counter += 1    # needed for the progress counter
                
                # e generator *already contain* lipid values for the three
                # lipid dataframe samples. So we now cycle through them (and
                # also check if they pass the allowed 0 averages tests).
                # next time, a class -.-'
                for lipi_sample in ("pla", "aor", "liv"):
                
                    failed_tests = 0
                
                    for lipi_genotype in ("b", "p", "l"):
                        l_handle = ("lipid_amount_" + 
                                  lipi_genotype +
                                  "_" +
                                  lipi_sample)
                              
                        #!!!
                                        
                        if average(e[l_handle]) == 0:
                            failed_tests += 1
                    
                    # cycled through all genotypes for this lipi sample.
                    if failed_tests > allowed_zero_lipid_averages:
                        continue    # to the next lipid sample
                    
                    # by the way, this next lipid sample is also yielded
                    # by the iterator by supplying another e (extracted values)
                    # dictionary
                
                    # here the current lipid sample has *passed* the tests
                    # for <allowed 0 averages>
                
                    # we now start cycling also @mirna_algo@ to perform tests
                    # <stdev_test_threshold>
                    # note: the test handled by stdev_test() further on checks
                    # *pairs* of average <--> SD, and tests if the *pair* passes
                    # the test. The pairs are given to the function as lists/
                    # tuples for convenience, but no relation is implied in the
                    # numbers. In this case, ALL pairs must pass the test, that
                    # is to say that the requirement to not exceed n times the SD
                    # must hold at both chow and Western diets
                    # 21/11/2016: but NOT different miRNA algorithms, of course.
                
                    for mirna_algo in ("desV2", "UQ"):
                    
                        mir_expr_aves = (
                            average(e[("mirna_expression_levels_b_" + mirna_algo)]),
                            average(e[("mirna_expression_levels_p_" + mirna_algo)]),
                            average(e[("mirna_expression_levels_l_" + mirna_algo)])
                            )
                    
                        lip_expr_aves = (
                            average(e[("lipid_amount_b_" + lipi_sample)]),
                            average(e[("lipid_amount_p_" + lipi_sample)]),
                            average(e[("lipid_amount_l_" + lipi_sample)])
                            )
                    
                        mir_expr_stds = (
                            std(e[("mirna_expression_levels_b_" + mirna_algo)]),
                            std(e[("mirna_expression_levels_p_" + mirna_algo)]),
                            std(e[("mirna_expression_levels_l_" + mirna_algo)])
                            )
                    
                        lip_expr_stds = (
                            std(e[("lipid_amount_b_" + lipi_sample)]),
                            std(e[("lipid_amount_p_" + lipi_sample)]),
                            std(e[("lipid_amount_l_" + lipi_sample)])
                            )
                    
                        if (stdev_test(
                            averages = (mir_expr_aves + lip_expr_aves),
                            standard_devs = (mir_expr_stds + lip_expr_stds),
                            threshold = stdev_test_threshold)
                            ) == False:
                        
                            continue
                    
                        # A little recap. At this point:
                        # - allowed 0 <lipi averages test> is *PASSED*
                        # - for the current miRNA, quantified with the current
                        #   algorithm, in the current organ,
                        #   associated to the current lipid, all biological
                        #   replicates *PASSED* the <standard deviation test>.
                    
                        # @lipi_sample@ ("pla", "aor", "liv")
                        # @mirna_algo@ ("desV2", "UQ")
                    
                        # I now must call (joined) the do_magic() generators,
                        # each with the different diets. Each will give *wise*
                        # x and y values from each diet, that I need to combine
                        # and, at last, perform the statistical testing.
                    
                        if debug:
                            print("**debug** before doing magic\n")
                            print("===============================================")
                            print(
                                "lengths test {}:".format(diet),
                                len(e["mirna_expression_levels_b_desV2"]),
                                len(e["mirna_expression_levels_p_desV2"]),
                                len(e["mirna_expression_levels_l_desV2"]),
                                len(e["mirna_expression_levels_b_UQ"]),
                                len(e["mirna_expression_levels_p_UQ"]),
                                len(e["mirna_expression_levels_l_UQ"]),
                                len(e["lipid_amount_b_aor"]),
                                len(e["lipid_amount_p_aor"]),
                                len(e["lipid_amount_l_aor"]),
                                len(e["lipid_amount_b_pla"]),
                                len(e["lipid_amount_p_pla"]),
                                len(e["lipid_amount_l_pla"]),
                                len(e["lipid_amount_b_liv"]),
                                len(e["lipid_amount_p_liv"]),
                                len(e["lipid_amount_l_liv"]),
                                len(e["mirna_headers_b_desV2"]),
                                len(e["mirna_headers_p_desV2"]),
                                len(e["mirna_headers_l_desV2"]),
                                len(e["mirna_headers_b_UQ"]),
                                len(e["mirna_headers_p_UQ"]),
                                len(e["mirna_headers_l_UQ"]),
                                len(e["lipid_headers_b_aor"]),
                                len(e["lipid_headers_p_aor"]),
                                len(e["lipid_headers_l_aor"]),
                                len(e["lipid_headers_b_pla"]),
                                len(e["lipid_headers_p_pla"]),
                                len(e["lipid_headers_l_pla"]),
                                len(e["lipid_headers_b_liv"]),
                                len(e["lipid_headers_p_liv"]),
                                len(e["lipid_headers_l_liv"])
                                )
                            print()
                            print("e =", e_cd, "\n")
                            print("algo = \"", mirna_algo, "\"", sep = "")
                            print("sample = \"", lipi_sample, "\"", sep = "")
                    
                        x, y, x_heads, y_heads = do_magic(
                                              e,
                                              algo = mirna_algo,
                                              sample = lipi_sample
                                              )
                    
                        if statistics_test == "pearson":
                            correlation = pearsonr(x, y)
                        elif statistics_test == "spearman":
                            correlation = spearmanr(x, y)
                        else:
                            raise TypeError("Something went wrong with "
                                            "statistics tests checks.")
                    
                        # returns (<stats test's value>, <p_value>)
                
                        if (correlation[1] >= pear_pval
                        or abs(correlation[0]) <= pear_corr_threshold):
                        
                            continue    # = statistics test failed

                        # test passed!! we're going to log it.
                        found += 1
                    
                        # adjusting some variables for logging
                    
                        all_values = x + y
                        all_headers = x_heads + y_heads
                    
                        output_string = (
                            output_string
                            + statsfriendly_singlediets_log(
                                d = e,
                                algorithm = mirna_algo,
                                sample = lipi_sample,
                                pvalue = str(correlation[1]),
                                stats_corr = str(correlation[0]),
                                all_values = all_values,
                                all_headers = all_headers,
                                sep = outfile_separator)
                            )
                    
                        if debug:
                            print("\n**debug**#######################################")
                            print("lipi_sample =", lipi_sample)
                            print("mirna_algo =", mirna_algo)
                            print("x_merged =", x_merged)
                            print("y_merged =", y_merged)
                            print("xh_merged =", xh_merged)
                            print("yh_merged =", yh_merged)
                            print("LEN x_merged ", len(x_merged))
                            print("LEN y_merged ", len(y_merged))
                            print("LEN all_values ", len(all_values))
                            print("LEN all_headers ", len(all_headers))
                            if correlation:
                                print("correlation = ", correlation)
                            print("**debug**#######################################\n")
                
                # done cycling in all lipid samples in e iterator
                # (each contains data from miRNA in one organ and
                # lipids from all three lipidomics)
                
                # faster progress display
                print("Processing ", counter, " of ", max_hits, 
                    "     found: ", found, "\r",
                    sep = "", end = "")
            
                #writes at each sample change in the iterator
                flush_to_file(output_string, outfile_pearson_filename)
                output_string = ""
    
    # iterators have now been extinguished
    print("Total", found, "comparisons passed the test.")
    
    global stop_seconds
    stop_seconds = time.time()
    
    if not argv.leave_untidy:
        tidy_dataframe(outfile_pearson_filename)        
    
    if logfile:
        logfile = logfile + "Total " + str(found) + "comparisons passed the test.\n"
        
    print("Analysis ended on", get_current_time().lower(), "         ")
    
    if logfile:
        logfile = logfile + "Analysis ended on: " + get_current_time().lower() + "\n"
# ////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////






# ====================================
def execute_statfriendly_mergediets(
    output_string,
    statistics_test = "pearson",
    debug = False,
    logfile = False
    ):
# ====================================
    
    if statistics_test == "pearson":
        print("\nEvaluating Pearson's correlations with these settings:")
    elif statistics_test == "spearman":
        print("\nEvaluating Spearman's correlations with these settings:")
    else:
        raise TypeError("Something went wrong with statistics tests checks.")

    print("Pval threshold: ", pear_pval,
          "\nPcorr threshold: ", pear_corr_threshold,
          "\nMax allowed standard deviation % from average: ",
          stdev_test_threshold, "%",
          "\nMax allowed zeroes in lipid DataFrame: ",
          allowed_zero_lipid_averages,
          sep = ""
          )
    
    max_hits = int((df_mirna.shape[0] * df_lipi.shape[0] * 2 * 7) / 2)
    # /2 because we're joining the diets
    counter = 0 # needed for the progress counter
    found = 0   # needed for the progress counter
    for organ in ("aor", "liv", "bra", "duo", "jej", "ile", "wat"):
        counter += 1    # needed for the progress counter
        # beginning to yield values cycling through miRNAs and lipids - 
        # at the same time for cd and wd !
        
        for e_cd, e_wd in zip(
            extract_values(
                df_mirna,
                df_lipi,
                diet = "cd",
                organ = organ,
                get_headers = True,
                do_strip_mirname = False
                ),
            extract_values(
                df_mirna,
                df_lipi,
                diet = "wd",
                organ = organ,
                get_headers = True,
                do_strip_mirname = False
                )
            ):
            
        
            counter += 1    # needed for the progress counter
            
            
            # e_cd and e_wd *already contain* lipid values for the three
            # lipid dataframe samples. So we now cycle through them (and
            # also check if they pass the allowed 0 averages tests).
            for lipi_sample in ("pla", "aor", "liv"):
                
                failed_tests = 0
                
                for lipi_genotype in ("b", "p", "l"):
                    l_handle = ("lipid_amount_" + 
                              lipi_genotype +
                              "_" +
                              lipi_sample)
                    
                    if average(e_cd[l_handle]) == 0:
                        failed_tests += 1
                        
                    if average(e_wd[l_handle]) == 0:
                        failed_tests += 1
                    
                # cycled through all genotypes for this lipi sample.
                if failed_tests > allowed_zero_lipid_averages:
                    continue    # to the next lipid sample
                
                # here the current lipid sample has *passed* the tests
                # for <allowed 0 averages>, and we're cycling through
                # note: test performed on *both* dicts (@cd, @wd)
                # @lipid_sample@
                
                # we now start cycling also @mirna_algo@ to perform tests
                # <stdev_test_threshold>
                # note: the test handled by stdev_test() further on checks
                # *pairs* of average <--> SD, and tests if the *pair* passes
                # the test. The pairs are given to the function as lists/
                # tuples for convenience, but no relation is implied in the
                # numbers. In this case, ALL pairs must pass the test, that
                # is to say that the requirement to not exceed n times the SD
                # must hold at both chow and Western diets.
                
                for mirna_algo in ("desV2", "UQ"):
                    
                    mir_expr_aves = (
                        average(e_cd[("mirna_expression_levels_b_" + mirna_algo)]),
                        average(e_cd[("mirna_expression_levels_p_" + mirna_algo)]),
                        average(e_cd[("mirna_expression_levels_l_" + mirna_algo)]),
                        average(e_wd[("mirna_expression_levels_b_" + mirna_algo)]),
                        average(e_wd[("mirna_expression_levels_p_" + mirna_algo)]),
                        average(e_wd[("mirna_expression_levels_l_" + mirna_algo)])
                        )
                    
                    lip_expr_aves = (
                        average(e_cd[("lipid_amount_b_" + lipi_sample)]),
                        average(e_cd[("lipid_amount_p_" + lipi_sample)]),
                        average(e_cd[("lipid_amount_l_" + lipi_sample)]),
                        average(e_wd[("lipid_amount_b_" + lipi_sample)]),
                        average(e_wd[("lipid_amount_p_" + lipi_sample)]),
                        average(e_wd[("lipid_amount_l_" + lipi_sample)])
                        )
                    
                    mir_expr_stds = (
                        std(e_cd[("mirna_expression_levels_b_" + mirna_algo)]),
                        std(e_cd[("mirna_expression_levels_p_" + mirna_algo)]),
                        std(e_cd[("mirna_expression_levels_l_" + mirna_algo)]),
                        std(e_wd[("mirna_expression_levels_b_" + mirna_algo)]),
                        std(e_wd[("mirna_expression_levels_p_" + mirna_algo)]),
                        std(e_wd[("mirna_expression_levels_l_" + mirna_algo)])
                        )
                    
                    lip_expr_stds = (
                        std(e_cd[("lipid_amount_b_" + lipi_sample)]),
                        std(e_cd[("lipid_amount_p_" + lipi_sample)]),
                        std(e_cd[("lipid_amount_l_" + lipi_sample)]),
                        std(e_wd[("lipid_amount_b_" + lipi_sample)]),
                        std(e_wd[("lipid_amount_p_" + lipi_sample)]),
                        std(e_wd[("lipid_amount_l_" + lipi_sample)])
                        )
                    
                    if (stdev_test(
                        averages = (mir_expr_aves + lip_expr_aves),
                        standard_devs = (mir_expr_stds + lip_expr_stds),
                        threshold = stdev_test_threshold)
                        ) == False:
                        
                        continue
                    
                    # A little recap. At this point:
                    # - allowed 0 <lipi averages test> is *PASSED*
                    # - for the current miRNA, quantified with the current
                    #   algorithm, in the current organ,
                    #   associated to the current lipid, all biological
                    #   replicates *PASSED* the <standard deviation test>.
                    #   ...for both diets, joined.
                    
                    # @lipi_sample@ ("pla", "aor", "liv")
                    # @mirna_algo@ ("desV2", "UQ")
                    
                    # I now must call (joined) the do_magic() generators,
                    # each with the different diets. Each will give *wise*
                    # x and y values from each diet, that I need to combine
                    # and, at last, perform the statistical testing.
                    
                    if debug:
                        print("**debug** before doing magic\n")
                        print("===============================================")
                        print(
                            "lengths test cd:",
                            len(e_cd["mirna_expression_levels_b_desV2"]),
                            len(e_cd["mirna_expression_levels_p_desV2"]),
                            len(e_cd["mirna_expression_levels_l_desV2"]),
                            len(e_cd["mirna_expression_levels_b_UQ"]),
                            len(e_cd["mirna_expression_levels_p_UQ"]),
                            len(e_cd["mirna_expression_levels_l_UQ"]),
                            len(e_cd["lipid_amount_b_aor"]),
                            len(e_cd["lipid_amount_p_aor"]),
                            len(e_cd["lipid_amount_l_aor"]),
                            len(e_cd["lipid_amount_b_pla"]),
                            len(e_cd["lipid_amount_p_pla"]),
                            len(e_cd["lipid_amount_l_pla"]),
                            len(e_cd["lipid_amount_b_liv"]),
                            len(e_cd["lipid_amount_p_liv"]),
                            len(e_cd["lipid_amount_l_liv"]),
                            len(e_cd["mirna_headers_b_desV2"]),
                            len(e_cd["mirna_headers_p_desV2"]),
                            len(e_cd["mirna_headers_l_desV2"]),
                            len(e_cd["mirna_headers_b_UQ"]),
                            len(e_cd["mirna_headers_p_UQ"]),
                            len(e_cd["mirna_headers_l_UQ"]),
                            len(e_cd["lipid_headers_b_aor"]),
                            len(e_cd["lipid_headers_p_aor"]),
                            len(e_cd["lipid_headers_l_aor"]),
                            len(e_cd["lipid_headers_b_pla"]),
                            len(e_cd["lipid_headers_p_pla"]),
                            len(e_cd["lipid_headers_l_pla"]),
                            len(e_cd["lipid_headers_b_liv"]),
                            len(e_cd["lipid_headers_p_liv"]),
                            len(e_cd["lipid_headers_l_liv"])
                            )
                        print(
                            "\nlengths test wd:",
                            len(e_wd["mirna_expression_levels_b_desV2"]),
                            len(e_wd["mirna_expression_levels_p_desV2"]),
                            len(e_wd["mirna_expression_levels_l_desV2"]),
                            len(e_wd["mirna_expression_levels_b_UQ"]),
                            len(e_wd["mirna_expression_levels_p_UQ"]),
                            len(e_wd["mirna_expression_levels_l_UQ"]),
                            len(e_wd["lipid_amount_b_aor"]),
                            len(e_wd["lipid_amount_p_aor"]),
                            len(e_wd["lipid_amount_l_aor"]),
                            len(e_wd["lipid_amount_b_pla"]),
                            len(e_wd["lipid_amount_p_pla"]),
                            len(e_wd["lipid_amount_l_pla"]),
                            len(e_wd["lipid_amount_b_liv"]),
                            len(e_wd["lipid_amount_p_liv"]),
                            len(e_wd["lipid_amount_l_liv"]),
                            len(e_wd["mirna_headers_b_desV2"]),
                            len(e_wd["mirna_headers_p_desV2"]),
                            len(e_wd["mirna_headers_l_desV2"]),
                            len(e_wd["mirna_headers_b_UQ"]),
                            len(e_wd["mirna_headers_p_UQ"]),
                            len(e_wd["mirna_headers_l_UQ"]),
                            len(e_wd["lipid_headers_b_aor"]),
                            len(e_wd["lipid_headers_p_aor"]),
                            len(e_wd["lipid_headers_l_aor"]),
                            len(e_wd["lipid_headers_b_pla"]),
                            len(e_wd["lipid_headers_p_pla"]),
                            len(e_wd["lipid_headers_l_pla"]),
                            len(e_wd["lipid_headers_b_liv"]),
                            len(e_wd["lipid_headers_p_liv"]),
                            len(e_wd["lipid_headers_l_liv"])
                            )
                        print()
                        print("e_cd =", e_cd, "\n")
                        print("e_wd =", e_wd, "\n")
                        print("algo = \"", mirna_algo, "\"", sep = "")
                        print("sample = \"", lipi_sample, "\"", sep = "")
                    
                    x_cd, y_cd, xh_cd, yh_cd = do_magic(e_cd,
                                                        algo = mirna_algo,
                                                        sample = lipi_sample
                                                        )
                    x_wd, y_wd, xh_wd, yh_wd = do_magic(e_wd,
                                                        algo = mirna_algo,
                                                        sample = lipi_sample
                                                        ) 
                    
                    x_merged = x_cd + x_wd
                    y_merged = y_cd + y_wd
                    xh_merged = xh_cd + xh_wd
                    yh_merged = yh_cd + yh_wd
                    
                    if statistics_test == "pearson":
                        correlation = pearsonr(x_merged, y_merged)
                    elif statistics_test == "spearman":
                        correlation = spearmanr(x_merged, y_merged)
                    else:
                        raise TypeError("Something went wrong with "
                                        "statistics tests checks.")
                    
                    # returns (<stats test's value>, <p_value>)
                
                    if (correlation[1] >= pear_pval
                    or abs(correlation[0]) <= pear_corr_threshold):
                        
                        continue    # = statistics test failed

                    # test passed!! we're going to log it.
                    found += 1
                    
                    # these will be used for logging
                    all_values = x_merged + y_merged
                    all_headers = xh_merged + yh_merged
                    
                    output_string = (
                        output_string
                        + statsfriendly_mergediets_log(
                            d = e_cd,
                            algorithm = mirna_algo,
                            sample = lipi_sample,
                            pvalue = str(correlation[1]),
                            stats_corr = str(correlation[0]),
                            all_values = all_values,
                            all_headers = all_headers,
                            sep = outfile_separator)
                        )
                    
                    if debug:
                        print("\n**debug**#######################################")
                        print("lipi_sample =", lipi_sample)
                        print("mirna_algo =", mirna_algo)
                        print("x_merged =", x_merged)
                        print("y_merged =", y_merged)
                        print("xh_merged =", xh_merged)
                        print("yh_merged =", yh_merged)
                        print("LEN x_merged ", len(x_merged))
                        print("LEN y_merged ", len(y_merged))
                        print("LEN all_values ", len(all_values))
                        print("LEN all_headers ", len(all_headers))
                        if correlation:
                            print("correlation = ", correlation)
                        print("**debug**#######################################\n")

            
            # faster progress display
            print("Processing ", counter, " of ", max_hits, 
                "     found: ", found, "\r",
                sep = "", end = "")
            
            #writes at each sample change in the iterator
            flush_to_file(output_string, outfile_pearson_filename)
            output_string = ""

    print("Total", found, "comparisons passed the test.")
    
    global stop_seconds
    stop_seconds = time.time()
    
    if not argv.leave_untidy:
        tidy_dataframe(outfile_pearson_filename)
    
    if logfile:
        logfile = logfile + "Total " + str(found) + "comparisons passed the test.\n"
        
    print("Analysis ended on", get_current_time().lower(), "         ")
    
    if logfile:
        logfile = logfile + "Analysis ended on: " + get_current_time().lower() + "\n"
# ////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////




# ====================================
def execute_pearsons_with_permutations(
    output_string,
    logfile = False
    ):
# ====================================
    
    print("\nEvaluating Pearson's correlations with these settings:\n"
            "Pval threshold: ", pear_pval,
            "\nPcorr threshold: ", pear_corr_threshold,
            "\nMax allowed standard deviation % from average: ",
            stdev_test_threshold, "%",
            "\nMax allowed zeroes in lipid DataFrame: ",
            allowed_zero_lipid_averages,
            sep = ""
            )
    
    max_hits = df_mirna.shape[0] * df_lipi.shape[0] * 2 * 7
    counter = 0 # needed for the progress counter
    found = 0   # needed for the progress counter
    for organ, diet in gen_cond():
        counter += 1    # needed for the progress counter
        # beginning to yield values cycling through miRNAs and lipids
        for e in extract_values(df_mirna,
                                df_lipi,
                                diet = diet,
                                organ = organ,
                                ):
            counter += 1    # needed for the progress counter
        
            # slow progress display
            #percent = round(((100 * counter)/max_hits), 2)
            #print("Processing.. ", counter, " of ", max_hits, "\t[", percent, "%]\t", "found: ", found, "\r", sep = "", end = "")
        
            # fast progress display
            #print("Processing ", counter, " of ", max_hits, "     found: ", found, "\r", sep = "", end = "")
        
        
            for perm in get_permutations(
                    e,
                    allowed_zero_lipid_averages = allowed_zero_lipid_averages
                    ):
            
                if perm[4] == "ZERO_TEST_FAILED":
                    continue
            
                pearsons_corr = pearsonr(perm[0], perm[1])
                # returns (<pearson's value>, <p_value>)
            
                if (pearsons_corr[1] >= pear_pval
                or abs(pearsons_corr[0]) <= pear_corr_threshold):
                    continue
        
                # preparing stuff for further testing:
            
                mir_expr_aves = (average(e[("mirna_expression_levels_b_" + perm[2])]),
                            average(e[("mirna_expression_levels_p_" + perm[2])]),
                            average(e[("mirna_expression_levels_l_" + perm[2])])
                            )
            
                lip_expr_aves = (average(e[("lipid_amount_b_" + perm[3])]),
                            average(e[("lipid_amount_p_" + perm[3])]),
                            average(e[("lipid_amount_l_" + perm[3])])
                            )
            
                mir_expr_stds = (std(e[("mirna_expression_levels_b_" + perm[2])]),
                            std(e[("mirna_expression_levels_p_" + perm[2])]),
                            std(e[("mirna_expression_levels_l_" + perm[2])])
                            )
            
                lip_expr_stds = (std(e[("lipid_amount_b_" + perm[3])]),
                            std(e[("lipid_amount_p_" + perm[3])]),
                            std(e[("lipid_amount_l_" + perm[3])])
                            )
            
                if (stdev_test(
                        averages = (mir_expr_aves + lip_expr_aves),
                        standard_devs = (mir_expr_stds + lip_expr_stds),
                        threshold = stdev_test_threshold)
                        ) == False:
                    continue
    
                found += 1
            
                # the following are needed in stat_log()
                e["pear_pval"] = pearsons_corr[1]
                e["pear_corr"] = pearsons_corr[0]
            
                # Current condition coordinates:
                # organ (mir)
                # diet (both)
                # perm[2] : either desV2 or UQ (algorithm, mir)
                # perm[3] : either pla, aor or liv (sample, lip)
            
                output_string = output_string + stat_log(d = e,
                    algorithm = perm[2],
                    sample = perm[3],
                    sep = outfile_separator,
                    mir_expr_aves = mir_expr_aves,
                    lip_expr_aves = lip_expr_aves,
                    mir_expr_stds = mir_expr_stds,
                    lip_expr_stds = lip_expr_stds,
                    )
            
            # faster progress display  
            print("Processing ", counter, " of ", max_hits, 
                    "     found: ", found, "\r",
                    sep = "", end = "")
        
            # writes at each diet or mirna organ change in the iterator
            flush_to_file(output_string, outfile_pearson_filename)
            output_string = ""

    print("Total", found, "comparisons passed the test.")
    
    global stop_seconds
    stop_seconds = time.time()

    if not argv.leave_untidy:
        tidy_dataframe(outfile_pearson_filename)
    
    if logfile:
        logfile = logfile + "Total " + str(found) + "comparisons passed the test.\n"
        
    print("Analysis ended on", get_current_time().lower(), "         ")
    
    if logfile:
        logfile = logfile + "Analysis ended on: " + get_current_time().lower() + "\n"
# ============================================================================
# ============================================================================


# ====================================
def execute_spearman_with_permutations(
    output_string,
    logfile = False
    ):
# ====================================
    print("\nEvaluating Spearman's correlations with these settings:\n"
            "Pval threshold: ", pear_pval,
            "\nPcorr threshold: ", pear_corr_threshold,
            "\nMax allowed standard deviation % from average: ",
            stdev_test_threshold, "%",
            "\nMax allowed zeroes in lipid DataFrame: ",
            allowed_zero_lipid_averages,
            sep = ""
            )
    
    max_hits = df_mirna.shape[0] * df_lipi.shape[0] * 2 * 7
    counter = 0 # needed for the progress counter
    found = 0   # needed for the progress counter
    for organ, diet in gen_cond():
        counter += 1    # needed for the progress counter
        # beginning to yield values cycling through miRNAs and lipids
        for e in extract_values(df_mirna,
                                df_lipi,
                                diet = diet,
                                organ = organ,
                                ):
            counter += 1    # needed for the progress counter
        
            # slow progress display
            #percent = round(((100 * counter)/max_hits), 2)
            #print("Processing.. ", counter, " of ", max_hits, "\t[", percent, "%]\t", "found: ", found, "\r", sep = "", end = "")
        
            # fast progress display
            #print("Processing ", counter, " of ", max_hits, "     found: ", found, "\r", sep = "", end = "")
        
        
            for perm in get_permutations(
                    e,
                    allowed_zero_lipid_averages = allowed_zero_lipid_averages
                    ):
            
                if perm[4] == "ZERO_TEST_FAILED":
                    continue
                
                #not anymore a pearsons, recycling old code only in this scope
                spearmans_corr = spearmanr(perm[0], perm[1])
                # returns (<pearson's value>, <p_value>)
            
                if (spearmans_corr[1] >= pear_pval
                or abs(spearmans_corr[0]) <= pear_corr_threshold):
                    continue
        
                # preparing stuff for further testing:
            
                mir_expr_aves = (average(e[("mirna_expression_levels_b_" + perm[2])]),
                            average(e[("mirna_expression_levels_p_" + perm[2])]),
                            average(e[("mirna_expression_levels_l_" + perm[2])])
                            )
            
                lip_expr_aves = (average(e[("lipid_amount_b_" + perm[3])]),
                            average(e[("lipid_amount_p_" + perm[3])]),
                            average(e[("lipid_amount_l_" + perm[3])])
                            )
            
                mir_expr_stds = (std(e[("mirna_expression_levels_b_" + perm[2])]),
                            std(e[("mirna_expression_levels_p_" + perm[2])]),
                            std(e[("mirna_expression_levels_l_" + perm[2])])
                            )
            
                lip_expr_stds = (std(e[("lipid_amount_b_" + perm[3])]),
                            std(e[("lipid_amount_p_" + perm[3])]),
                            std(e[("lipid_amount_l_" + perm[3])])
                            )
            
                if (stdev_test(
                        averages = (mir_expr_aves + lip_expr_aves),
                        standard_devs = (mir_expr_stds + lip_expr_stds),
                        threshold = stdev_test_threshold)
                        ) == False:
                    continue
    
                found += 1
            
                # the following are needed in stat_log()
                e["pear_pval"] = spearmans_corr[1]
                e["pear_corr"] = spearmans_corr[0]
            
                # Current condition coordinates:
                # organ (mir)
                # diet (both)
                # perm[2] : either desV2 or UQ (algorithm, mir)
                # perm[3] : either pla, aor or liv (sample, lip)
            
                output_string = output_string + stat_log(d = e,
                    algorithm = perm[2],
                    sample = perm[3],
                    sep = outfile_separator,
                    mir_expr_aves = mir_expr_aves,
                    lip_expr_aves = lip_expr_aves,
                    mir_expr_stds = mir_expr_stds,
                    lip_expr_stds = lip_expr_stds,
                    )
            
            # faster progress display  
            print("Processing ", counter, " of ", max_hits, 
                    "     found: ", found, "\r",
                    sep = "", end = "")
        
            # writes at each diet or mirna organ change in the iterator
            flush_to_file(output_string, outfile_pearson_filename)
            output_string = ""
    print("Total", found, "comparisons passed the test.")
    
    if logfile:
        logfile = logfile + "Total " + str(found) + "comparisons passed the test.\n"

    global stop_seconds
    stop_seconds = time.time()
    
    if not argv.leave_untidy:
        tidy_dataframe(outfile_pearson_filename)
    
    print("Analysis ended on", get_current_time().lower(), "         ")
    
    if logfile:
        logfile = logfile + "Analysis ended on: " + get_current_time().lower() + "\n"
# ============================================================================
# ============================================================================


# ====================================
def execute_pearson_with_averages(
    output_string,
    logfile = False
    ):
# ====================================
    print("\nEvaluating Pearson's correlations with these settings:\n"
            "Pval threshold: ", pear_pval,
            "\nPcorr threshold: ", pear_corr_threshold,
            "\nMax allowed standard deviation % from average: ",
            stdev_test_threshold, "%",
            "\nMax allowed zeroes in lipid DataFrame: ",
            allowed_zero_lipid_averages,
            sep = ""
            )
    
    max_hits = df_mirna.shape[0] * df_lipi.shape[0] * 2 * 7
    counter = 0 # needed for the progress counter
    found = 0   # needed for the progress counter
    for organ, diet in gen_cond():
        counter += 1    # needed for the progress counter
        # beginning to yield values cycling through miRNAs and lipids
        for e in extract_values(df_mirna,
                                df_lipi,
                                diet = diet,
                                organ = organ,
                                ):
            counter += 1    # needed for the progress counter
        
            # slow progress display
            #percent = round(((100 * counter)/max_hits), 2)
            #print("Processing.. ", counter, " of ", max_hits, "\t[", percent, "%]\t", "found: ", found, "\r", sep = "", end = "")
        
            # fast progress display
            #print("Processing ", counter, " of ", max_hits, "     found: ", found, "\r", sep = "", end = "")
        
        
            for avz in get_averages(
                    e,
                    allowed_zero_lipid_averages = allowed_zero_lipid_averages
                    ):
                
                if avz[4] == "ZERO_TEST_FAILED":
                    continue
                
                pearsons_corr = pearsonr(avz[0], avz[1])
                # returns (<spearman's value>, <p_value>)
                
                if (pearsons_corr[1] >= pear_pval
                or abs(pearsons_corr[0]) <= pear_corr_threshold):
                    continue
                
                # preparing stuff for further testing:
                
                # averages are already there
                mir_expr_aves = avz[0]
                lip_expr_aves = avz[1]
                
                # standard deviations must be calculated    
                mir_expr_stds = (std(e[("mirna_expression_levels_b_" + avz[2])]),
                            std(e[("mirna_expression_levels_p_" + avz[2])]),
                            std(e[("mirna_expression_levels_l_" + avz[2])])
                            )
            
                lip_expr_stds = (std(e[("lipid_amount_b_" + avz[3])]),
                            std(e[("lipid_amount_p_" + avz[3])]),
                            std(e[("lipid_amount_l_" + avz[3])])
                            )
            
                if (stdev_test(
                        averages = (mir_expr_aves + lip_expr_aves),
                        standard_devs = (mir_expr_stds + lip_expr_stds),
                        threshold = stdev_test_threshold)
                        ) == False:
                    continue
    
                found += 1
            
                # the following are needed in stat_log()
                e["pear_pval"] = pearsons_corr[1]
                e["pear_corr"] = pearsons_corr[0]
            
                # Current condition coordinates:
                # organ (mir)
                # diet (both)
                # avz[2] : either desV2 or UQ (algorithm, mir)
                # avz[3] : either pla, aor or liv (sample, lip)
            
                output_string = output_string + stat_log(d = e,
                    algorithm = avz[2],
                    sample = avz[3],
                    sep = outfile_separator,
                    mir_expr_aves = mir_expr_aves,
                    lip_expr_aves = lip_expr_aves,
                    mir_expr_stds = mir_expr_stds,
                    lip_expr_stds = lip_expr_stds,
                    )
            
            # faster progress display  
            print("Processing ", counter, " of ", max_hits, 
                    "     found: ", found, "\r",
                    sep = "", end = "")
        
            # writes at each diet or mirna organ change in the iterator
            flush_to_file(output_string, outfile_pearson_filename)
            output_string = ""
            
    print("Total", found, "comparisons passed the test.")
    
    global stop_seconds
    stop_seconds = time.time()
    
    if not argv.leave_untidy:
        tidy_dataframe(outfile_pearson_filename)
    
    if logfile:
        logfile = logfile + "Total " + str(found) + "comparisons passed the test.\n"
    
    print("Analysis ended on", get_current_time().lower(), "         ")
    
    if logfile:
        logfile = logfile + "Analysis ended on: " + get_current_time().lower() + "\n"
# ============================================================================
# ============================================================================


# ====================================
def execute_spearman_with_averages(
    output_string,
    logfile = False
    ):
# ====================================
    print("\nEvaluating Spearman's correlations with these settings:\n"
            "Pval threshold: ", pear_pval,
            "\nPcorr threshold: ", pear_corr_threshold,
            "\nMax allowed standard deviation % from average: ",
            stdev_test_threshold, "%",
            "\nMax allowed zeroes in lipid DataFrame: ",
            allowed_zero_lipid_averages,
            sep = ""
            )
    
    max_hits = df_mirna.shape[0] * df_lipi.shape[0] * 2 * 7
    counter = 0 # needed for the progress counter
    found = 0   # needed for the progress counter
    for organ, diet in gen_cond():
        counter += 1    # needed for the progress counter
        # beginning to yield values cycling through miRNAs and lipids
        for e in extract_values(df_mirna,
                                df_lipi,
                                diet = diet,
                                organ = organ,
                                ):
            counter += 1    # needed for the progress counter
        
            # slow progress display
            #percent = round(((100 * counter)/max_hits), 2)
            #print("Processing.. ", counter, " of ", max_hits, "\t[", percent, "%]\t", "found: ", found, "\r", sep = "", end = "")
        
            # fast progress display
            #print("Processing ", counter, " of ", max_hits, "     found: ", found, "\r", sep = "", end = "")
        
        
            for avz in get_averages(
                    e,
                    allowed_zero_lipid_averages = allowed_zero_lipid_averages
                    ):
                
                if avz[4] == "ZERO_TEST_FAILED":
                    continue
                
                spearmans_corr = spearmanr(avz[0], avz[1])
                # returns (<spearman's value>, <p_value>)
                
                if (spearmans_corr[1] >= pear_pval
                or abs(spearmans_corr[0]) <= pear_corr_threshold):
                    continue
                
                # preparing stuff for further testing:
                
                # averages are already there
                mir_expr_aves = avz[0]
                lip_expr_aves = avz[1]
                
                # standard deviations must be calculated    
                mir_expr_stds = (std(e[("mirna_expression_levels_b_" + avz[2])]),
                            std(e[("mirna_expression_levels_p_" + avz[2])]),
                            std(e[("mirna_expression_levels_l_" + avz[2])])
                            )
            
                lip_expr_stds = (std(e[("lipid_amount_b_" + avz[3])]),
                            std(e[("lipid_amount_p_" + avz[3])]),
                            std(e[("lipid_amount_l_" + avz[3])])
                            )
            
                if (stdev_test(
                        averages = (mir_expr_aves + lip_expr_aves),
                        standard_devs = (mir_expr_stds + lip_expr_stds),
                        threshold = stdev_test_threshold)
                        ) == False:
                    continue
    
                found += 1
            
                # the following are needed in stat_log()
                e["pear_pval"] = spearmans_corr[1]
                e["pear_corr"] = spearmans_corr[0]
            
                # Current condition coordinates:
                # organ (mir)
                # diet (both)
                # avz[2] : either desV2 or UQ (algorithm, mir)
                # avz[3] : either pla, aor or liv (sample, lip)
            
                output_string = output_string + stat_log(d = e,
                    algorithm = avz[2],
                    sample = avz[3],
                    sep = outfile_separator,
                    mir_expr_aves = mir_expr_aves,
                    lip_expr_aves = lip_expr_aves,
                    mir_expr_stds = mir_expr_stds,
                    lip_expr_stds = lip_expr_stds,
                    )
            
            # faster progress display  
            print("Processing ", counter, " of ", max_hits, 
                    "     found: ", found, "\r",
                    sep = "", end = "")
        
            # writes at each diet or mirna organ change in the iterator
            flush_to_file(output_string, outfile_pearson_filename)
            output_string = ""
            
    print("Total", found, "comparisons passed the test.")
    
    global stop_seconds
    stop_seconds = time.time()
    
    if not argv.leave_untidy:
        tidy_dataframe(outfile_pearson_filename)
    
    if logfile:
        logfile = logfile + "Total " + str(found) + "comparisons passed the test.\n"
    
    print("Analysis ended on", get_current_time().lower(), "         ")
    
    if logfile:
        logfile = logfile + "Analysis ended on: " + get_current_time().lower() + "\n"
# ============================================================================
# ============================================================================


# =============================================================================
# ======///          END OF RECONCILER-SPECIFIC RUN MODES         ///==========
# =============================================================================

# GLOBAL COMPILED REGEXES:
match_desv2 = re.compile(r".*_desV2$")
match_uq = re.compile(r".*_UQ$")
match_lipi_cond = re.compile(r"(..)_(.*)_(\d*)(.)_([^_]*)")
match_mirna_cond = re.compile(r"(.)(\d+)(.*)_(.[^_]*)_(.[^_]*)_(.[^_]*)")



# Initializing output string. These are the headers of the table.
if (argv.merge_diets and runmode == "statfriendly"):
    output_string = ("mirna" + outfile_separator +
                     "organ" + outfile_separator +
                     "algorithm" + outfile_separator +
                     "lipid" + outfile_separator +
                     "sample" + outfile_separator +
                     "pvalue" + outfile_separator +
                     "stats_corr" + outfile_separator +
                     "mir_b1_cd_ID" + outfile_separator +
                     "mir_b2_cd_ID" + outfile_separator +
                     "mir_b3_cd_ID" + outfile_separator +
                     "mir_p1_cd_ID" + outfile_separator +
                     "mir_p2_cd_ID" + outfile_separator +
                     "mir_p3_cd_ID" + outfile_separator +
                     "mir_l1_cd_ID" + outfile_separator +
                     "mir_l2_cd_ID" + outfile_separator +
                     "mir_l3_cd_ID" + outfile_separator +
                     "mir_b1_wd_ID" + outfile_separator +
                     "mir_b2_wd_ID" + outfile_separator +
                     "mir_b3_wd_ID" + outfile_separator +
                     "mir_p1_wd_ID" + outfile_separator +
                     "mir_p2_wd_ID" + outfile_separator +
                     "mir_p3_wd_ID" + outfile_separator +
                     "mir_l1_wd_ID" + outfile_separator +
                     "mir_l2_wd_ID" + outfile_separator +
                     "mir_l3_wd_ID" + outfile_separator +
                     "lip_b1_cd_ID" + outfile_separator +
                     "lip_b2_cd_ID" + outfile_separator +
                     "lip_b3_cd_ID" + outfile_separator +
                     "lip_p1_cd_ID" + outfile_separator +
                     "lip_p2_cd_ID" + outfile_separator +
                     "lip_p3_cd_ID" + outfile_separator +
                     "lip_l1_cd_ID" + outfile_separator +
                     "lip_l2_cd_ID" + outfile_separator +
                     "lip_l3_cd_ID" + outfile_separator +
                     "lip_b1_wd_ID" + outfile_separator +
                     "lip_b2_wd_ID" + outfile_separator +
                     "lip_b3_wd_ID" + outfile_separator +
                     "lip_p1_wd_ID" + outfile_separator +
                     "lip_p2_wd_ID" + outfile_separator +
                     "lip_p3_wd_ID" + outfile_separator +
                     "lip_l1_wd_ID" + outfile_separator +
                     "lip_l2_wd_ID" + outfile_separator +
                     "lip_l3_wd_ID" + outfile_separator +
                     "mir_b1_cd_val" + outfile_separator +
                     "mir_b2_cd_val" + outfile_separator +
                     "mir_b3_cd_val" + outfile_separator +
                     "mir_p1_cd_val" + outfile_separator +
                     "mir_p2_cd_val" + outfile_separator +
                     "mir_p3_cd_val" + outfile_separator +
                     "mir_l1_cd_val" + outfile_separator +
                     "mir_l2_cd_val" + outfile_separator +
                     "mir_l3_cd_val" + outfile_separator +
                     "mir_b1_wd_val" + outfile_separator +
                     "mir_b2_wd_val" + outfile_separator +
                     "mir_b3_wd_val" + outfile_separator +
                     "mir_p1_wd_val" + outfile_separator +
                     "mir_p2_wd_val" + outfile_separator +
                     "mir_p3_wd_val" + outfile_separator +
                     "mir_l1_wd_val" + outfile_separator +
                     "mir_l2_wd_val" + outfile_separator +
                     "mir_l3_wd_val" + outfile_separator +
                     "lip_b1_cd_val" + outfile_separator +
                     "lip_b2_cd_val" + outfile_separator +
                     "lip_b3_cd_val" + outfile_separator +
                     "lip_p1_cd_val" + outfile_separator +
                     "lip_p2_cd_val" + outfile_separator +
                     "lip_p3_cd_val" + outfile_separator +
                     "lip_l1_cd_val" + outfile_separator +
                     "lip_l2_cd_val" + outfile_separator +
                     "lip_l3_cd_val" + outfile_separator +
                     "lip_b1_wd_val" + outfile_separator +
                     "lip_b2_wd_val" + outfile_separator +
                     "lip_b3_wd_val" + outfile_separator +
                     "lip_p1_wd_val" + outfile_separator +
                     "lip_p2_wd_val" + outfile_separator +
                     "lip_p3_wd_val" + outfile_separator +
                     "lip_l1_wd_val" + outfile_separator +
                     "lip_l2_wd_val" + outfile_separator +
                     "lip_l3_wd_val" + "\n" 
                     )
elif (argv.merge_diets == False and runmode == "statfriendly"):
    output_string = ("mirna" + outfile_separator +
                     "organ" + outfile_separator +
                     "algorithm" + outfile_separator +
                     "lipid" + outfile_separator +
                     "sample" + outfile_separator +
                     "pvalue" + outfile_separator +
                     "stats_corr" + outfile_separator +
                     "mir_b1_ID" + outfile_separator +
                     "mir_b2_ID" + outfile_separator +
                     "mir_b3_ID" + outfile_separator +
                     "mir_p1_ID" + outfile_separator +
                     "mir_p2_ID" + outfile_separator +
                     "mir_p3_ID" + outfile_separator +
                     "mir_l1_ID" + outfile_separator +
                     "mir_l2_ID" + outfile_separator +
                     "mir_l3_ID" + outfile_separator +
                     "lip_b1_ID" + outfile_separator +
                     "lip_b2_ID" + outfile_separator +
                     "lip_b3_ID" + outfile_separator +
                     "lip_p1_ID" + outfile_separator +
                     "lip_p2_ID" + outfile_separator +
                     "lip_p3_ID" + outfile_separator +
                     "lip_l1_ID" + outfile_separator +
                     "lip_l2_ID" + outfile_separator +
                     "lip_l3_ID" + outfile_separator +
                     "mir_b1_val" + outfile_separator +
                     "mir_b2_val" + outfile_separator +
                     "mir_b3_val" + outfile_separator +
                     "mir_p1_val" + outfile_separator +
                     "mir_p2_val" + outfile_separator +
                     "mir_p3_val" + outfile_separator +
                     "mir_l1_val" + outfile_separator +
                     "mir_l2_val" + outfile_separator +
                     "mir_l3_val" + outfile_separator +
                     "lip_b1_val" + outfile_separator +
                     "lip_b2_val" + outfile_separator +
                     "lip_b3_val" + outfile_separator +
                     "lip_p1_val" + outfile_separator +
                     "lip_p2_val" + outfile_separator +
                     "lip_p3_val" + outfile_separator +
                     "lip_l1_val" + outfile_separator +
                     "lip_l2_val" + outfile_separator +
                     "lip_l3_val" + "\n" 
                     )
else:
    output_string = ("mirna" + outfile_separator +
                     "organ" + outfile_separator +
                     "diet" + outfile_separator +
                     "algorithm" + outfile_separator +
                     "lipid" + outfile_separator +
                     "sample" + outfile_separator +
                     "pvalue" + outfile_separator +
                     "stats_corr" + outfile_separator +
                     "mir_ave_b" + outfile_separator +
                     "mir_ave_p" + outfile_separator +
                     "mir_ave_l" + outfile_separator +
                     "lip_ave_b" + outfile_separator +
                     "lip_ave_p" + outfile_separator +
                     "lip_ave_l" + outfile_separator +
                     "mir_std_b" + outfile_separator +
                     "mir_std_p" + outfile_separator +
                     "mir_std_l" + outfile_separator +
                     "lip_std_b" + outfile_separator +
                     "lip_std_p" + outfile_separator +
                     "lip_std_l" + outfile_separator +
                     "mir_b_1" + outfile_separator +
                     "mir_b_2" + outfile_separator +
                     "mir_b_3" + outfile_separator +
                     "mir_p_1" + outfile_separator +
                     "mir_p_2" + outfile_separator +
                     "mir_p_3" + outfile_separator +
                     "mir_l_1" + outfile_separator +
                     "mir_l_2" + outfile_separator +
                     "mir_l_3" + outfile_separator +
                     "lip_b_1" + outfile_separator +
                     "lip_b_2" + outfile_separator +
                     "lip_b_3" + outfile_separator +
                     "lip_b_4" + outfile_separator +
                     "lip_b_5" + outfile_separator +
                     "lip_b_6" + outfile_separator +
                     "lip_p_1" + outfile_separator +
                     "lip_p_2" + outfile_separator +
                     "lip_p_3" + outfile_separator +
                     "lip_p_4" + outfile_separator +
                     "lip_p_5" + outfile_separator +
                     "lip_p_6" + outfile_separator +
                     "lip_l_1" + outfile_separator +
                     "lip_l_2" + outfile_separator +
                     "lip_l_3" + outfile_separator +
                     "lip_l_4" + outfile_separator +
                     "lip_l_5" + outfile_separator +
                     "lip_l_6" + "\n"
                     )

print("Hello, I'm Reconciler Advanced, version ",
         __version__, ".", sep = "")


# printing stuff, logging stuff
print("Starting analysis on", get_current_time().lower())
logfile = logfile + "Analysis started: " + get_current_time().lower() + "\n"
logfile = logfile + "Version: " + __version__ + "\n"
print("\nLoading lipi and mirna tables..")

df_mirna = pd.read_csv(
    mirna_dataframe, sep = "\t"
    ).replace(
        to_replace = "NaN", value = 0
        )

print("Mirna data table:", df_mirna.shape[0],
      "rows,", df_mirna.shape[1], "cols.")
print("Mirna data table filename:", mirna_dataframe)

logfile = (logfile + "Mirna data table: "
           + str(df_mirna.shape[0])
           + " rows, "
           + str(df_mirna.shape[1])
           + " cols.\n"
           )

df_lipi = pd.read_csv(
    lipid_dataframe, sep = "\t"
    ).replace(
        to_replace = "NaN", value = nan_value
        )
# nan_value was defined in globals

print("Lipi data table:", df_lipi.shape[0],
      "rows,", df_lipi.shape[1], "cols.")
print("Lipi data table filename:", lipid_dataframe)

logfile = (logfile + "Lipi data table: "
           + str(df_lipi.shape[0])
           + " rows, "
           + str(df_lipi.shape[1])
           + " cols.\n"
           )


# ===============
# setting runmode
# ===============

start_seconds = time.time()

# Various runmodes have been kept separate to avoid to check at every
# iteration all conditions. Not a groundbreaking speed bump, but a good
# compromise to keep things as simple as possible

if runmode == "strict" and statistics_test == "pearson":
    print("Running in strict mode (tests are run on all x-to-y permutations).")
    execute_pearsons_with_permutations(
        output_string = output_string
        )
elif runmode == "strict" and statistics_test == "spearman":
    print("Running in strict mode (tests are run on all x-to-y permutations).")
    execute_spearman_with_permutations(
        output_string = output_string
        )
elif runmode == "average" and statistics_test == "pearson":
    print("Running in average mode (tests are run on value averages).")
    execute_pearson_with_averages(
        output_string = output_string
        )
elif runmode == "average" and statistics_test == "spearman":
    print("Running in average mode (tests are run on value averages).")
    execute_spearman_with_averages(
        output_string = output_string
        )
elif runmode == "statfriendly":
    if argv.merge_diets:
        print("Doing magic in statfriendly mode.")
        execute_statfriendly_mergediets(
            output_string = output_string,
            statistics_test = statistics_test
            )
    else:
        print("Doing magic in statfriendly mode. Separate diets.")
        execute_statfriendly_singlediets(
            output_string = output_string,
            statistics_test = statistics_test
            )

benchmark = round(stop_seconds - start_seconds, 2)
print("Analysis took {} seconds to perform.".format(benchmark))
print("(Import times and dataframe tidying are not taken into account.)")
