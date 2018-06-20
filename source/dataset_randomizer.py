#!/Applications/Anaconda3/anaconda/bin/python3
# Author: Manzini Stefano; stefano.manzini@gmail.com


__version__ = "1.2.131016"

"""This program:

- is a modified version of Reconciler's engine aimed at getting info about
    variability of data in the two datasets.
    Of course, this is done to randomize new datasets that closely resemble
    the parent ones, but random.
    
*** LICENSE ***
This program comes under:
Creative Commons Attribution-ShareAlike 4.0 International License
Full license here: http://creativecommons.org/licenses/by-sa/4.0/legalcode

Following here is a human-readable summary of (and NOT a substitute for)
the License:

Briefly, it grants you the rights to copy and distribute the work if you agree
to acknowledge to always cite the author and to always make clear to others the 
terms of this license.

You are free to make derivative works and use both derivative
and the original work for commercial purposes as long
as you distribute all that under under the same license
as the original.
"""

# ==================================================
# ============     arguments parsing    ============
# ==================================================

import argparse
import sys
parser = argparse.ArgumentParser(
    description = ("Dataset Randomizer takes two Cariplo 2.1-compliant miRNA "
            "and lipid dataframes, and produces new ones with the same data "
            "shape as the parent dataframes, but random."
            )
    )

parser.add_argument(
        "-m",
        "--mirna-dataframe",
        help = ("Sets the parent miRNA dataframe"),
        dest = "mirna_dataframe",
        metavar = "",
        )

parser.add_argument(
        "-l",
        "--lipid-dataframe",
        help = ("Sets the parent lipid dataframe"),
        dest = "lipid_dataframe",
        metavar = "",
        )

parser.add_argument(
        "-n",
        "--number-to-make",
        help = ("Tells Dataset Randomizer how many dataframes to make (both)"
                ". Defaults to 1"),
        dest = "dataframes_to_yeld",
        metavar = "",
        )

parser.add_argument(
        "-s", 
        "--sep",
        help = ("Sets the preferred csv separator. Defaults to <tab>."),
        metavar = "",
        dest = "outfile_separator")

parser.add_argument(
    "-v",
    "--version",
    action = "version",
    version = __version__)

argv = parser.parse_args()

# =========================================================
# =========///     end of arguments parsing    ///=========
# =========================================================



# ===============================
# ===         imports         ===
# ===============================

print("Importing modules..")
import sys

from numpy import mean, std
from io import StringIO

# the following is safe, no namespace will be harmed.
from manzutils import *


# ===============================
# ===         globals         ===
# ===============================

if argv.dataframes_to_yeld:
    dataframes_to_yeld = int(argv.dataframes_to_yeld)
else:
    dataframes_to_yeld = 1

if argv.mirna_dataframe:
    mirna_dataframe = argv.mirna_dataframe
else:
    mirna_dataframe = "miRNAs_dataset.csv"

if argv.lipid_dataframe:
    lipid_dataframe = argv.lipid_dataframe
else:
    lipid_dataframe = "lipidomics_dataset.csv"

if argv.outfile_separator:
    sep = argv.outfile_separator
else:
    sep = "\t"

# value to change NaN values in lipi dataframe to. Use a NUMBER! Not string.
nan_value = 0  # lipid dataframe only!

# rounding of the decimals of % SD deviation (to avoid dict overcrowding)
#sd_rounding = 0

allowed_decimals = 4    # this is for not getting dicts too sparse

# ======================================
# ===         end of globals         ===
# ======================================


# =====================================================================
# ===                           FUNCTIONS                           ===
# =====================================================================

# == DEBUG ==
def nan_test(something):
# == DEBUG ==
    """Debug function. Quits if "nan" found.
    """
    
    if something == "nan":
        sys.exit()


# ====================
def dict_to_list(dict):
# ====================
    
    """ This function is designed to work with this script.
    Takes appropriate dicts (where each key is a str(<number>) and each
    value is the number of times that was found, and returns a list
    populated with that many times those values.
    
    This also handles nasty "nan" values.
    """
    
    working_dict = dict.copy()
    result = []
    
    try:
        for i in range(working_dict["nan"]):
            result.append(0.0)
            
        del working_dict["nan"]
    except:
        pass
        
    for k, v in working_dict.items():
        for i in range(working_dict[k]):
            result.append(float(k))
            
    return result



# ==================================================
def extract_mir_values(df_mirna,
                   diet = None,
                   organ = None
                   ):
# ==================================================
    
    """ This is a *modified function* from manzutils's extract_values()
    modified to just work with one special dataframe at a time.
    Following doc is related to the original function:
    
    ***
    
    It digs into mirnomics and lipidomics dataframes, and yields, at
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
    
    extracted_values = {
        "mirna_expression_levels_b_desV2" : [],
        "mirna_expression_levels_p_desV2" : [],
        "mirna_expression_levels_l_desV2" : [],
        "mirna_expression_levels_b_UQ" : [],
        "mirna_expression_levels_p_UQ" : [],
        "mirna_expression_levels_l_UQ" : [],
        "curr_mirna_ID" : "",
        "curr_diet" : "",
        "curr_organ" : "",
        "curr_algorithm" : ""
    }
    
    c = (organ, None, diet)    # current c (condition) parameters
    # written like this to match indexes from mirna_ids() used later on
    
    # c[0] = organ
    # c[2] = diet
    
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
        
        df_mirna_slice = df_mirna[i : i + 1]  # still has headers and all!
        
        # =============================================================
        # appending current miRNA ID name
        extracted_values["curr_mirna_ID"] = strip_mirname(df_mirna_slice["mirna"][i])
        # fetching expression levels and putting them into extracted_values dict
        
        for mh_hit in mh_hits:
            if re.match(r".*_desV2$", mh_hit):
                if mirna_ids(mh_hit)[1] == "b":
                    extracted_values["mirna_expression_levels_b_desV2"].append(df_mirna_slice[mh_hit][i])
                elif mirna_ids(mh_hit)[1] == "p":
                    extracted_values["mirna_expression_levels_p_desV2"].append(df_mirna_slice[mh_hit][i])
                elif mirna_ids(mh_hit)[1] == "l":
                    extracted_values["mirna_expression_levels_l_desV2"].append(df_mirna_slice[mh_hit][i])
            elif re.match(r".*_UQ$", mh_hit):
                if mirna_ids(mh_hit)[1] == "b":
                    extracted_values["mirna_expression_levels_b_UQ"].append(df_mirna_slice[mh_hit][i])
                elif mirna_ids(mh_hit)[1] == "p":
                    extracted_values["mirna_expression_levels_p_UQ"].append(df_mirna_slice[mh_hit][i])
                elif mirna_ids(mh_hit)[1] == "l":
                    extracted_values["mirna_expression_levels_l_UQ"].append(df_mirna_slice[mh_hit][i])
        # ======== done fetching expression levels for current miRNA ==========
        
        yield extracted_values



# ==================================================
def extract_lip_values(
                   df_lipi,
                   diet = None,
                   organ = None
                   ):
# ==================================================
    
    """ This is a *modified function* from manzutils's extract_values()
    modified to just work with one special dataframe at a time.
    Following doc is related to the original function:
    
    ***
    
    It digs into mirnomics and lipidomics dataframes, and yields, at
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
    
    extracted_values = {
        "lipid_amount_b_pla" : [],
        "lipid_amount_p_pla" : [],
        "lipid_amount_l_pla" : [],
        "lipid_amount_b_aor" : [],
        "lipid_amount_p_aor" : [],
        "lipid_amount_l_aor" : [],
        "lipid_amount_b_liv" : [],
        "lipid_amount_p_liv" : [],
        "lipid_amount_l_liv" : [],
        "curr_lipid_ID" : "",
        "curr_diet" : "",
        "curr_organ" : ""
    }
    
    extracted_values["curr_diet"] = diet
    extracted_values["curr_organ"] = organ
    
    c = (organ, None, diet) # useful to lipi_ids()

    lh_hits = []    # container for lipid headers
    for lh in df_lipi.columns[1:]:
        lh_ids = lipi_ids(lh) # (sample, genotype, diet, animal)
        
        # only diet matters now. We're going to pick aorta, liver and
        # plasma lipid levels regardless of the organ we're coming from
        # the miRNA table.
        # That would be 6(values per condition)*3(genotypes)*3(samples)=
        # = 54 total values (out of 108 total headers).
        
        if c[2] == lh_ids[2]:   # diet
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
            

# ============================================================================
# ===                         end of   FUNCTIONS                           ===
# ============================================================================


        




# SDs will be counted in here
mastermind = {
    "mir_bpl_var_desV2" : {},
    "mir_biorep_var_desV2" : {},
    "mir_bpl_var_UQ" : {},
    "mir_biorep_var_UQ" : {},
    "lip_bpl_var" : {},
    "lip_biorep_var" : {},
    "mir_all_averages_desV2" : {},
    "mir_all_averages_UQ" : {},
    "lip_all_averages" : {}
    }


print("Hello, I'm Dataset Randomizer, version ",
         __version__, ".", sep = "")


# printing stuff, logging stuff
print("Starting analysis on", get_current_time().lower())
print("\nLoading lipi and mirna tables..")

df_mirna = pd.read_csv(mirna_dataframe, sep = "\t").replace(to_replace = "NaN", value = 0)
print("Mirna data table:", df_mirna.shape[0], "rows,", df_mirna.shape[1], "cols.")
print("Mirna data table filename:", mirna_dataframe)

df_lipi = pd.read_csv(lipid_dataframe, sep = "\t").replace(to_replace = "NaN", value = nan_value)
# nan_value was defined in globals
print("Lipi data table:", df_lipi.shape[0], "rows,", df_lipi.shape[1], "cols.")
print("Lipi data table filename:", lipid_dataframe)




print("\nGetting info on data variability..")

# =============================================================================
# starting with miRNA DataFrame
# =============================================================================
counter = 0 # needed for the progress counter
for organ, diet in gen_cond():

    print("\nWorking on:", organ, diet)
    counter += 1    # needed for the progress counter
    # beginning to yield values cycling through miRNAs and lipids
    
    for e in extract_mir_values(df_mirna,
                            diet = diet,
                            organ = organ,
                            ):
        counter += 1    # needed for the progress counter
    
        print(counter, "\r", end = "")
                
        for a in ("desV2", "UQ"):   # cycling algorithms
            minimind = {
                "b_temp" : 0,
                "p_temp" : 0,
                "l_temp" : 0
                }
            for g in ("b_", "p_", "l_"):    # cycling genotypes
                
                # getting info inter biological replicates
                
                # all biological replicates averages
                global_average_key = "mir_all_averages_" + a
                
                # inter-genotype dict handle:
                biorep_key = "mir_biorep_var_" + a
                
                # needed later for intra-genotype averages:
                minimind_key = g + "temp"  
                
                ave = mean(e["mirna_expression_levels_" + g + a])
                ave_perc_diffs = [(x / ave) for x in e["mirna_expression_levels_" + g + a]]
                
                # applying the decimals threshold, also adjusting type
                ave = str(round(float(ave), allowed_decimals))
                ave_perc_diffs = [str(round(float(x), allowed_decimals)) for x in ave_perc_diffs]
                
                mastermind[global_average_key].setdefault(ave, 0)
                mastermind[global_average_key][ave] += 1
                
                for i in ave_perc_diffs:
                    mastermind[biorep_key].setdefault(i, 0)
                    mastermind[biorep_key][i] += 1
                
                minimind[minimind_key] = ave

            
            # getting info on the differences among b,p,l averages of same
            # organ, diet and algorithm    
            
            # dict handle for % variability within genotypes
            bpl_key = "mir_bpl_var_" + a
            
            ave = average([float(x) for x in minimind.values()])
            ave_perc_diffs = [(float(x)/ave) for x in minimind.values()]
            
            # applying the decimals threshold, also adjusting type
            
            ave_perc_diffs = [str(round(float(x), allowed_decimals)) for x in ave_perc_diffs]
            
            for i in ave_perc_diffs:
                mastermind[bpl_key].setdefault(i, 0)
                mastermind[bpl_key][i] += 1

# =============================================================================
#starting again, with lipi DataFrame
# =============================================================================

for organ, diet in gen_cond():

    print("\nWorking on:", organ, diet)
    
    for e in extract_lip_values(df_lipi,
                            diet = diet,
                            organ = organ,
                            ):
        counter += 1    # needed for the progress counter
        
        print(counter, "\r", end = "")
        
        # lipi
        
        # considering each sample separately but putting all results together
        for s in ("pla", "aor", "liv"):

            minimind = {
                "b_temp" : 0,
                "p_temp" : 0,
                "l_temp" : 0
                }
                
            for g in ("b_", "p_", "l_"):
                
                # =============================================================
                # getting info inter biological replicates
                # =============================================================
                
                # invariant key: "lip_all_averages"
                # invariant key: "lip_biorep_var"
                
                # needed later for intra-genotype averages:
                minimind_key = g + "temp"  
                
                ave = mean(e["lipid_amount_" + g + s])
                ave_perc_diffs = [(x /ave) for x in e["lipid_amount_" + g + s]]
                
                # applying the decimals threshold, also adjusting type
                ave = str(round(float(ave), allowed_decimals))
                ave_perc_diffs = [str(round(float(x), allowed_decimals)) for x in ave_perc_diffs]
                
                mastermind["lip_all_averages"].setdefault(ave, 0)
                mastermind["lip_all_averages"][ave] += 1
                
                for i in ave_perc_diffs:
                    mastermind["lip_biorep_var"].setdefault(i, 0)
                    mastermind["lip_biorep_var"][i] += 1
                
                minimind[minimind_key] = ave
            
            # =================================================================
            # getting info on the differences among b,p,l averages of same
            # organ, diet and algorithm
            # =================================================================
            
            # dict handle for % variability within genotypes
            # invariant key: bpl_key = "lip_bpl_var"
            
            ave = average([float(x) for x in minimind.values()])
            ave_perc_diffs = [(float(x) /ave) for x in minimind.values()]
            
            # applying the decimals threshold, also adjusting type
            
            ave_perc_diffs = [str(round(float(x), allowed_decimals)) for x in ave_perc_diffs]
            
            for i in ave_perc_diffs:
                mastermind["lip_bpl_var"].setdefault(i, 0)
                mastermind["lip_bpl_var"][i] += 1

print("Populating values lists..")
masterlist = {}

for k in mastermind.keys():
    masterlist[k] = dict_to_list(mastermind[k])

#=======================
#making a miRNA DataFrame
#=======================

for index in range(dataframes_to_yeld):
    
    print("Building a random mirna DataFrame:", index)
    
    fake_headers_list = ["mirna"]
    for f in gen_fake_mirna_cond_full_sex():
        fake_headers_list.append(f)
    
    fake_headers_list.append("remove_me")
    fake_mir_df = sep.join(fake_headers_list) + "\n"
    
    mirlist = [("mmu-miR-" + str(n)
        + "-0p_ID=MIM0000000_MI0000000_mmu-mir-"
        + str(n)) for n in range(df_mirna.shape[0])]
    
    algorithm_switcher = 0  # trick to avoid getting mad with headers
    for f in mirlist:
        fake_mir_df = fake_mir_df + f + sep
        for n in range(
                0,
                (len(fake_headers_list) - 2),
                9   # amplitude of miRNA bio replicates * number of genotypes
                ):
            if algorithm_switcher < 14:
                a = "desV2"
                ave_key = "mir_all_averages_" + a
                biorep_key = "mir_biorep_var_" + a
                bpl_key = "mir_bpl_var_" + a
                
                seed = choice(masterlist[ave_key])
                fake_mir_df = fake_mir_df + str(seed) + sep
                modifier_1 = choice(masterlist[biorep_key])
                fake_mir_df = fake_mir_df + str(seed * modifier_1) + sep
                modifier_2 = choice(masterlist[biorep_key])
                fake_mir_df = fake_mir_df + str(seed * modifier_2) + sep
                
                # new seed for new genotype
                seed = seed * choice(masterlist[bpl_key])
                fake_mir_df = fake_mir_df + str(seed) + sep
                modifier_1 = choice(masterlist[biorep_key])
                fake_mir_df = fake_mir_df + str(seed * modifier_1) + sep
                modifier_2 = choice(masterlist[biorep_key])
                fake_mir_df = fake_mir_df + str(seed * modifier_2) + sep
                
                # new seed for last genotype
                seed = seed * choice(masterlist[bpl_key])
                fake_mir_df = fake_mir_df + str(seed) + sep
                modifier_1 = choice(masterlist[biorep_key])
                fake_mir_df = fake_mir_df + str(seed * modifier_1) + sep
                modifier_2 = choice(masterlist[biorep_key])
                fake_mir_df = fake_mir_df + str(seed * modifier_2) + sep
                            
                algorithm_switcher += 1
            else:
                a = "UQ"
                ave_key = "mir_all_averages_" + a
                biorep_key = "mir_biorep_var_" + a
                bpl_key = "mir_bpl_var_" + a
                
                seed = choice(masterlist[ave_key])
                fake_mir_df = fake_mir_df + str(seed) + sep
                modifier_1 = choice(masterlist[biorep_key])
                fake_mir_df = fake_mir_df + str(seed * modifier_1) + sep
                modifier_2 = choice(masterlist[biorep_key])
                fake_mir_df = fake_mir_df + str(seed * modifier_2) + sep
                
                # new seed for new genotype
                seed = seed * choice(masterlist[bpl_key])
                fake_mir_df = fake_mir_df + str(seed) + sep
                modifier_1 = choice(masterlist[biorep_key])
                fake_mir_df = fake_mir_df + str(seed * modifier_1) + sep
                modifier_2 = choice(masterlist[biorep_key])
                fake_mir_df = fake_mir_df + str(seed * modifier_2) + sep
                
                # new seed for last genotype
                seed = seed * choice(masterlist[bpl_key])
                fake_mir_df = fake_mir_df + str(seed) + sep
                modifier_1 = choice(masterlist[biorep_key])
                fake_mir_df = fake_mir_df + str(seed * modifier_1) + sep
                modifier_2 = choice(masterlist[biorep_key])
                fake_mir_df = fake_mir_df + str(seed * modifier_2) + sep
                
        fake_mir_df = fake_mir_df + "\n"
    
    fmir_df = pd.read_csv(StringIO(fake_mir_df), sep = "\t")
    del fmir_df["remove_me"]    # this eliminates the supernumerary column
    del fake_headers_list
    del fake_mir_df
    
    out_mirna_file_name = "rnd_mirna_" + str(index) + ".csv"
    fmir_df.to_csv(out_mirna_file_name, index = False, sep = sep)

# =======================
# making a lipid DataFrame
# =======================

for index in range(dataframes_to_yeld):
    print("Building a random lipid DataFrame:", index)
    
    fake_headers_list = ["LIPID_NAME"]
    for f in gen_fake_lipid_cond_full_sex():
        fake_headers_list.append(f)
    
    fake_headers_list.append("remove_me")
    fake_lip_df = sep.join(fake_headers_list) + "\n"
    
    liplist = [("fakenoic acid " + str(n)) for n in range(df_lipi.shape[0])]
    
    # so that I can recycle code later on:
    ave_key = "lip_all_averages"
    biorep_key = "lip_biorep_var"
    bpl_key = "lip_bpl_var"
    
    bio_reps = 6
            
    for f in liplist:
        fake_lip_df = fake_lip_df + f + sep
        for n in range(
                0,
                (len(fake_headers_list) - 2),
                bio_reps   # amplitude of lipid biological replicates
                ):
            seed = choice(masterlist[ave_key])
            fake_lip_df = fake_lip_df + str(seed) + sep
            for t in range(bio_reps - 1):
                modifier = choice(masterlist[biorep_key])
                fake_lip_df = fake_lip_df + str(seed * modifier) + sep
        
        fake_lip_df = fake_lip_df + "\n"
    
    flip_df = pd.read_csv(StringIO(fake_lip_df), sep = "\t")
    del flip_df["remove_me"]    # this eliminates the supernumerary column
    del fake_headers_list
    del fake_lip_df
    
    out_lipi_file_name = "rnd_lipi_" + str(index) + ".csv"
    flip_df.to_csv(out_lipi_file_name, index = False, sep = sep)
