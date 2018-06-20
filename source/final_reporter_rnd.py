#!/Applications/Anaconda3/anaconda/bin/python3
# Author: Manzini Stefano; stefano.manzini@gmail.com

# This script is derived from final_reporter.py, and it should be used to
# get a final report of random reconciler.py reports because final_reporter.py
# also categorize lipids by classes, but random lipids have no class.

__version__ = "0.8.130217"

"""
This program is meant to read tables produced by Reconciler.
It gives back enrichment tables that count:

- how many times a miRNA is associated with a passed test
- how many times a lipid is associated with a passed test
- detailed stats and enrichment info about miRNAs and lipids

/* known bugs!!
 * when writing the final detailed lipid report, there are chances of
 * IndexError if some tests are fewer than expected. To be fixed!!
 */


################
# DEPENDENCIES #
################

To get Numpy, Scipy, etc out of the box: install Anaconda:

https://www.continuum.io/downloads

Almost all other dependencies can be installed via

pip install <module>

Google for further info!
"""


# ==================================================
# ============     arguments parsing    ============
# ==================================================

import argparse
import sys
parser = argparse.ArgumentParser(
    description = ("Final Reporter produces "
                   "enrichment reports for given Reconciler-produced table.")
    )

parser.add_argument(
        "dataframe",
        help = "Filename of Reconciler DataFrame.")

parser.add_argument(
        "-s", 
        "--sep",
        help = ("Sets the preferred csv separator. Defaults to <tab>."),
        metavar = "",
        dest = "outfile_separator")
        
parser.add_argument(
        "-m", 
        "--mirna-header",
        help = ("Sets a custom mirna header ID. Defaults to \"mirna\"."),
        metavar = "",
        dest = "mirna_header")
        
parser.add_argument(
        "-l", 
        "--lipid-header",
        help = ("Sets a custom lipid header ID. Defaults to \"lipid\"."),
        metavar = "",
        dest = "lipid_header")

parser.add_argument(
        "-t", 
        "--shared-mirnas-threshold",
        help = ("Sets the number of miRNA names, shared by a whole lipid "
                "class, to include in the final report. Defaults to 20."),
        metavar = "",
        dest = "shared_mirnas_number_threshold")

parser.add_argument(
        "-c", 
        "--report-most-common-mirnas",
        help = ("Sets the number of most shared miRNAs to detail in the final "
                "report. For each, the name and the number of occurrence is "
                "detailed. Defaults to 10."),
        metavar = "",
        dest = "most_common_mirnas_to_report")

# new in 0.8!
parser.add_argument(
    "-M",
    "--tidy-for-mirna-report",
    help = """
        This option is meant for working with sliced Reconciler Tables. If
        selected, it moves all sliced tables to /sliced/, to the exception
        of *_mir_finalreport.csv, kept in the working directory for later
        processing by do_tidy_reported.py (which is also invoked).
        """,
    action = "store_true",
)

parser.add_argument(
    "-v",
    "--version",
    action = "version",
    version = __version__)

argv = parser.parse_args()

# =========================================================
# ==========//     end of arguments parsing    //==========
# =========================================================


# ================================
# =====     FUNCTIONS     ========
# ================================

# =========================================
def gen_class_heads(class_list, sep = "\t"):
# =========================================
    
    """Writes class-specific headers.
    
    Expects a <list> of strings: class_list
    
    Returns a <string>.
    """
    
    returnlist = []
    
    common_heads = [    # list of parameters specific for *each* lipid class
        "_ave_corr",    # average of all correlation values
        "_ave_abs_corr",    # average of modulus of all correlation values
        "_ave_pos_corr", # average of all positive correlation values
        "_ave_neg_corr", # average of all negative correlation values
        "_curr_elems",  # number of lipids in this class that passed test
        "_total",   # total number of lipids in input dataframe
        "_passed",  # total number of lipids that passed *any* test
        "_total%",  # % passed lipids for this mirna, over total lipids
        "_passed%", # % passed lipids for this mirna, over *any* passed lipid
        "_pos", # number of lipids with positive correlation
        "_neg", # number of lipids with negative correlation
        "_pos%",    # % positive passed lipids, over total passed lipids
        "_neg%" # % negative passed lipids, over total passed lipids
        ]
    
    for c in class_list:
        for h in common_heads:
            returnlist.append(c + h)
    
    return sep.join(returnlist)



# =====================
def m_average(listlike):
# =====================
    
    """Special wrapper over numpy average function that deals with
    empty list/vectors.
    """
    
    if len(listlike) > 0:
        return average(listlike)
    else:
        return 0
    

# =======================================
# ===//     END OF FUNCTIONS     //======
# =======================================


# ==============================
# =====     imports     ========
# ==============================
print("Importing modules..")

from numpy import mean, std
                  
from scipy.stats import (pearsonr, # (<pearson’s coefficient>, <2-tail-p-value)
                         spearmanr)
# the following is safe, no namespace will be harmed.
# required manzutils 0.17.300316 at least
from manzutils import *

import minibar
import pickle
from io import StringIO
import os
import subprocess

# these are imported with manzutils
#import regex as re    # for regular expressions
#import pandas as pd
#import time
# from numpy import (average,
#                   polyfit,
#                   poly1d,
#                   linspace)

#import plotly.tools as tls
#tls.set_credentials_file(username = "Manz", api_key = "2gmj2v0hxj")
#import plotly.plotly as py
#import plotly.graph_objs as go

# =====================================
# ===//     end of imports     //======
# =====================================

# ==============================
# =====     GLOBALS     ========
# ==============================

dataframe = argv.dataframe

# input/output file separator
if argv.outfile_separator is None:
    sep = "\t"
else:
    sep = argv.outfile_separator

if argv.mirna_header:
    mirna_header = argv.mirna_header
else:
    mirna_header = "mirna"
    
if argv.lipid_header:
    lipid_header = argv.lipid_header
else:
    lipid_header = "lipid"

if argv.shared_mirnas_number_threshold:
    shared_mirnas_number_threshold = int(argv.shared_mirnas_number_threshold)
else:
    shared_mirnas_number_threshold = 20

if argv.most_common_mirnas_to_report:
    if int(argv.most_common_mirnas_to_report) > 0:
        most_common_mirnas_to_report = int(argv.most_common_mirnas_to_report)
    else:
        most_common_mirnas_to_report = 1
else:
    most_common_mirnas_to_report = 10


# hard coded lipid species. These are all lipid species present in
# megatable_lipi.csv, as extracted by get_lipid_class2()
# (manzutils 0.24+)
all_lipid_classes = ['CE', 'Cer', 'DAG', 'FC', 'Gb3', 'Glc', 'LPC', 'LPE', 
                     'LPI', 'Lac', 'PA', 'PC', 'PC O', 'PC P', 'PE', 'PE O',
                     'PE P', 'PG', 'PI', 'PS', 'SM', 'TAG']

# ===============     defining output filenames      ===============

out_filename_prefix = ""
if dataframe.endswith(".csv"):
    out_filename_prefix = dataframe.replace(".csv", "")
elif dataframe.endswith(".txt"):
    out_filename_prefix = dataframe.replace(".txt", "")
else:
    out_filename_prefix = dataframe

out_filename_mir = (
        out_filename_prefix
        + "_miR_enrich.csv"
        )

out_filename_lip = (
        out_filename_prefix
        + "_lip_enrich.csv"
        )

out_filename_finalreport_mir = (
        out_filename_prefix
        + "_mir_finalreport.csv"
        )

out_filename_finalreport_lip = (
        out_filename_prefix
        + "_lip_finalreport.csv"
        )
# =============//     end of output filenames      //=============

# ========================================
# =====//    END OF GLOBALS     //========
# ========================================


print("Working on:\n" + argv.dataframe)

if dataframe.endswith(".xlsx"):
    df = pd.read_excel(dataframe)
else:
    df = pd.read_csv(dataframe, sep = sep)

print("Gathering miRNAs enrichment..")
mir_enrichment_dict = count_repetitive_column_elements(df, mirna_header)
write_dict_to_csv(mir_enrichment_dict, out_filename_mir)

print("Gathering lipid enrichment..")
lip_enrichment_dict = count_repetitive_column_elements(df, lipid_header)
write_dict_to_csv(lip_enrichment_dict, out_filename_lip)


# now getting serious #

mm = {}


print("Gathering detailed miRNAs enrichment..")

mm["mirna"] = {}

# creating a handle
mmirna = mm["mirna"]

i = 0   # goodbye, enumerate(). Needed by minibar
for mirna in minibar.bar(df["mirna"]):

    mmirna.setdefault(mirna, {})

    # setting current mirna mm dict handle
    mmc = mmirna[mirna]

    # populating current miRNA dict, if not already done
    # same key as table headers
    mmc.setdefault("lipid", [])
    mmc.setdefault("pvalue", [])
    mmc.setdefault("stats_corr", [])
    mmc.setdefault("class", [])
    # new elements
    mmc.setdefault("lipids_classes", set())

    # taking notes #
    
    mmc["lipid"].append(df["lipid"][i])
    mmc["pvalue"].append(df["pvalue"][i])
    mmc["stats_corr"].append(df["stats_corr"][i])
    #mmc["class"].append(get_lipid_class2(df["lipid"][i]))
    mmc["class"].append("NO_CLASS")

    #mmc["lipids_classes"].add(get_lipid_class2(df["lipid"][i]))
    mmc["lipids_classes"].add("NO_CLASS")

    i += 1


print("\nGathering detailed lipid enrichment..")

mm["lipid"] = {}

# creating a handle
mmlip = mm["lipid"]

i = 0   # goodbye, enumerate(). Needed by minibar
for lipid in minibar.bar(df["lipid"]):
    
    # trying to fix indexes
    
    mmlip.setdefault(lipid, {})
    
    # setting current mirna mm dict handle
    mml = mmlip[lipid]
    
    # populating current lipid dict, if not already done
    # same key as table headers
    mml.setdefault("mirna", [])
    mml.setdefault("pvalue", [])
    mml.setdefault("stats_corr", [])
    # taking notes #
    
    mml["mirna"].append(df["mirna"][i])
    mml["pvalue"].append(df["pvalue"][i])
    mml["stats_corr"].append(df["stats_corr"][i])
    
    i += 1

print("\nDone indexing.")


# to report, we need to know stuff about the original lipid database, and
# also stuff of the Reconciler-produced table.

try:
    all_lipids_and_classes_info = pickle.load(
        open("all_lipids_and_classes_info.p", "rb")
    )
except FileNotFoundError:
    all_lipids_and_classes_info = pickle.load(
        open("./source/all_lipids_and_classes_info.p", "rb")
    )

#reconciled_lipids_and_classes_info = index_lipid_species(df, "lipid")


# ============================================================================#
# =======================   now dealing with miRNAs    =======================#
# ============================================================================#

# headers! #
headers = (
    "mirna" + sep +
    "n_lipids" + sep +
    "n_classes" + sep +
    "ave_pval" + sep
    )

headers = headers + gen_class_heads(all_lipid_classes)

headers = headers + sep + "delete_me" # supernumerary column to be deleted at the end

# starting to write output table as string#

returnstring = ""
returnstring = returnstring + headers + "\n"

print("\nProcessing detailed miRNA enrichment..")
for m in minibar.bar(mmirna.keys()):    # for every miRNA key
    
    returnstring = (returnstring +
                   m + sep +
                   str(len(mmirna[m]["lipid"])) + sep +
                   str(len(mmirna[m]["lipids_classes"])) + sep +
                   str(m_average(mmirna[m]["pvalue"])) + sep
                   )
    
    mmc = mmirna[m]  # current mirna dict handle
    
    for c in all_lipid_classes: # for every possible lipid class
        
        if c in mmc["class"]:
            
            # for the current class, we need to pull and process info regarding
            # that class only.
            
            # this must be done element-wise; can't do that with comprehension.
        
            lipids = []
            pvalues = []
            stats_corrs = []
        
            positive_lipids = []
            negative_lipids = []
        
            positive_pvalues = []
            negative_pvalues = []
        
            positive_stats_corrs = []
            negative_stats_corrs = []
        
        
            for elem in enumerate(mmc["class"]):
                if elem[1] == c:
                    lipids.append(mmc["lipid"][elem[0]])
                    pvalues.append(mmc["pvalue"][elem[0]])
                    stats_corrs.append(mmc["stats_corr"][elem[0]])
                
                    if mmc["stats_corr"][elem[0]] >= 0:
                        positive_lipids.append(mmc["lipid"][elem[0]])
                        positive_pvalues.append(mmc["pvalue"][elem[0]])
                        positive_stats_corrs.append(mmc["stats_corr"][elem[0]])
                    else:
                        negative_lipids.append(mmc["lipid"][elem[0]])
                        negative_pvalues.append(mmc["pvalue"][elem[0]])
                        negative_stats_corrs.append(mmc["stats_corr"][elem[0]])
            
            # DONE populating current class lists for writing to table
            
            all_lipids_number = len(all_lipids_and_classes_info[c])
            passed_lipids_number = len(
                set(reconciled_lipids_and_classes_info[c])
                )
            
            # PLEASE NOTE
            # Reconciler-produced tables ma contain any number of duplicate
            # lipids. This can result if multiple organ/tissue intersections
            # are represented at once. Thus, for a given miRNA, a lipid
            # can appear multiple times if it passed the test with that miRNA
            # in more than one organ/tissue combination.
            # To avoid confusions with total possible lipid species, the
            # numbers are calculated *at worst*, which means that for a given
            # miRNA multiple positivities are *not* considered.
            # i.e. If there are 30 individual ceramides, it may be that a
            # given miRNA has 60 passed ceramides, but they are not unique,
            # coming from pooled organ/tissue combinations.
            # These are reduced to unique elements prior to evaluating.
            #
            # *Vice versa*, positive and negative correlation scores are
            # *not* rendered unique, because in different organ/tissue
            # combinations it may happen that the very same lipid has a
            # different correlation; to avoid confusion, every hit
            # is regarded as an independent event.
            
            returnstring = (
                returnstring +
                str(m_average(stats_corrs)) + sep +
                str(m_average(absolute(stats_corrs))) + sep +
                str(m_average(positive_stats_corrs)) + sep +
                str(m_average(negative_stats_corrs)) + sep +
                
                #str(len(lipids)) + sep +
                str(len(set(lipids))) + sep +
                
                str(all_lipids_number) + sep +    # from GLOBALS
                str(passed_lipids_number) + sep +
                
                #str(round((len(lipids) / all_lipids_number) * 100, 1)) + sep +
                str(round((len(set(lipids)) / all_lipids_number) * 100, 1)) + sep +
                
                #str(round((len(lipids) / passed_lipids_number) * 100, 1)) + sep +
                str(round((len(set(lipids)) / passed_lipids_number) * 100, 1)) + sep +
                
                #str(len(positive_lipids)) + sep +
                str(len(set(positive_lipids))) + sep +
                
                #str(len(negative_lipids)) + sep +
                str(len(set(negative_lipids))) + sep +
                
                str(round((len(positive_lipids) / len(lipids)) * 100, 1)) + sep +
                # str(round((len(set(positive_lipids)) / len(set(lipids))) * 100, 1)) + sep +
                
                str(round((len(negative_lipids) / len(lipids)) * 100, 1)) + sep
                # str(round((len(set(negative_lipids)) / len(set(lipids))) * 100, 1)) + sep
                )
        
        else:   # if for this miRNA there's no lipid for this class
                # manually put values (all zeroes?)
            
            all_lipids_number = len(all_lipids_and_classes_info[c])
            
            try:
                passed_lipids_number = len(
                    set(reconciled_lipids_and_classes_info[c])
                    )
            except:
                passed_lipids_number = "n.a."
            
            returnstring = (
                returnstring +
                "n.d." + sep +
                "n.d." + sep +
                "n.d." + sep +
                "n.d." + sep +
                "0" + sep +
                str(all_lipids_number) + sep +    # from GLOBALS
                str(passed_lipids_number) + sep +
                "0" + sep +
                "0" + sep +
                "0" + sep +
                "0" + sep +
                "0" + sep +
                "0" + sep
                )
    
    # JOB done for all lipid classes of current miRNA
    returnstring = returnstring + "\n"

# JOB done for all miRNAs in the dictionary

# time to write out the final table

out_mirna_df = pd.read_csv(StringIO(returnstring), sep = sep)
del out_mirna_df["delete_me"]

# now ordering out_mirna_df by header: "n_lipids", descending
out_mirna_df = out_mirna_df.sort_values(by = "n_lipids", ascending = False)

out_mirna_df.to_csv(out_filename_finalreport_mir, sep = sep, index = False)
print("\nFinal detailed miRNA report written to:\n", out_filename_finalreport_mir)

#                         // JOB DONE with miRNAS //

print("Job done with miRNAs. Stopping here. Fake lipids are not divided into classes.")
sys.exit()


# ============================================================================#
# =======================   now dealing with lipids    =======================#
# ============================================================================#

# gathering together, class-wise, miRNAs that are shared by *all* lipids
# of that class
mm["mirnas_shared"] = {}
mmsham = mm["mirnas_shared"]    # handle

for l in mmlip.keys():
    
    lipid_class = get_lipid_class2(l)
    temp_mirna_set = set(mmlip[l]["mirna"])
    
    mmsham.setdefault(lipid_class, temp_mirna_set)
    mmsham[lipid_class] = mmsham[lipid_class] & temp_mirna_set

# ==============================================================
# now counting the number of times miRNAs hit per class elements
mm["mirnas_per_lipid_class"] = {}
mmclam = mm["mirnas_per_lipid_class"]   # handle

for l in mmlip.keys():
    
    lipid_class = get_lipid_class2(l)
    mmclam.setdefault(lipid_class, {})
    
    temp_mirna_list = mmlip[l]["mirna"]
    
    for m in temp_mirna_list:
        
        mmclam[lipid_class].setdefault(m, 0)
        mmclam[lipid_class][m] += 1
# =========// end of making mmclam // ==========================


# making string containing lipid headers #
lheaders = (
    "lipid" + sep +
    "n_mirnas" + sep +  # n. of mirnas correlating with current lipid
    "ave_pval" + sep +  # average of all pvalues
    "ave_corr" + sep +  # average of all correlation values (pos and neg)
    "ave_abs_corr" + sep +  # average of absolute correlation values
    "ave_pos_corr" + sep +  # average of positive correlation values
    "ave_neg_corr" + sep +  # average of negative correlation values
    "tot_lipids_in_class" + sep +   # total lipids in this class
    "tot_passed_lipids_in_class" + sep +    # total passed lipids in this class
    "common_mirnas_n" + sep +   # number miRNAs shared by all elements of class
    "common_mirnas" + sep   # shared miRNAs, if fewer than threshold
    )

# adding headers specific to the intended number of most common miRNAs, for
# that class, to report. For each one, it is reported:
#
# - the name of the miRNA
# - the number of occurrence within lipids that passed test in that class
#
# * always * bear in mind that the number of passed tests can exceed the number
# of actual lipids, since Reconciler tables can contain multiple organ/sample
# comparisons.

temphead = ""
for e in range(most_common_mirnas_to_report):
    currmir = "mir_" + str(e + 1)
    temphead = (temphead + currmir + sep + (currmir + "_hits") + sep)
temphead = temphead + "delete_me"

lheaders = lheaders + temphead
# lipid headers are now complete

# starting to write output table as string#

returnstring = ""   # come on in, garbage collector
returnstring = returnstring + lheaders + "\n"

# this try block, inserted in 0.8, addresses the issue of having
# the "bug with few lipids" that still *needs* t be addressed, and
# having final_reporter go on when the -M flag is selected.
# FIX THIS BUG!
try:
    print("\nProcessing detailed lipid enrichment..")
    for l in minibar.bar(mmlip.keys()):    # for every lipid key
    
        mml = mmlip[l]  # current lipid dict handle
        current_class = get_lipid_class2(l)
    
        returnstring = (returnstring +
                        l + sep +
                        str(len(mml["mirna"])) + sep +
                        str(m_average(mml["pvalue"])) + sep +
                        str(m_average(mml["stats_corr"])) + sep +
                        str(m_average(absolute(mml["stats_corr"]))) + sep
                        )
    
        positive_stats_corr = [x for x in mml["stats_corr"] if x >= 0]
        negative_stats_corr = [x for x in mml["stats_corr"] if x < 0]
    
        returnstring = (
            returnstring +
            str(m_average(positive_stats_corr)) + sep + 
            str(m_average(negative_stats_corr)) + sep +
            str(len(all_lipids_and_classes_info[current_class])) + sep +
            str(len(reconciled_lipids_and_classes_info[current_class])) + sep +
            str(len(mmsham[current_class])) + sep
            )
    
        if len(mmsham[current_class]) <= shared_mirnas_number_threshold:
            separator = ";"
            shared_mirnas = separator.join(mmsham[current_class])
        else:
            shared_mirnas = "more than " + str(shared_mirnas_number_threshold)
    
        returnstring = returnstring + shared_mirnas + sep
    
    
        # now populating the most common miRNAs
        sorted_keys = sorted(
            mmclam[current_class],
            key = lambda k: mmclam[current_class][k],
            reverse = True
            )

        for cm in range(most_common_mirnas_to_report):
            returnstring = (returnstring +
                            sorted_keys[cm] + sep +
                            str(mmclam[current_class][sorted_keys[cm]]) + sep
                            )
    
        returnstring = returnstring + "\n"
    
    out_lipid_df = pd.read_csv(StringIO(returnstring), sep = sep)
    del out_lipid_df["delete_me"]
    out_lipid_df.to_csv(out_filename_finalreport_lip, sep = sep, index = False)
    print("\nFinal detailed lipid report written to:\n", out_filename_finalreport_mir)

    #                         // JOB DONE with lipids //
except IndexError:
    print("\nBUG: fewer lipids thatn expected, this needs debugging")

mm_dump_filename = out_filename_prefix + ".p"
pickle.dump(mm, open(mm_dump_filename, "wb"))
print(
    "\n"
    "Gathered and ordered data has been pickled to a Python dict object"
    " for further inspection; data has been dumped to:\n"
    + mm_dump_filename)

# new in 0.8!
#out_filename_mir
#out_filename_lip
#out_filename_finalreport_mir # the one we don't want to move
#out_filename_finalreport_lip
if argv.tidy_for_mirna_report:
    print("Tidying all but *mir_finalreport files. Putting everything in /sliced ..")

    if not os.path.exists("sliced"):
        os.mkdir("sliced")

    if os.path.isdir("sliced"):
        path = "./sliced/"
        for file in (
            out_filename_mir,
            out_filename_lip,
            out_filename_finalreport_lip,
            mm_dump_filename,
        ):
            destpath = path + file
            try:
                os.rename(file, destpath)
            except FileNotFoundError:
                print("Expected file {} not found, skipping..".format(file))

        # now calling ./do_tidy_reported.py
        print("Calling do_tidy_reported.py..")
        command = "python do_tidy_reported.py"
        subprocess.call(command, shell = True)
    else:   # "sliced" path exists but is not a directory
        print("Path /sliced exists, but is not a directory. Stop tidying now.")
        sys.exit()
