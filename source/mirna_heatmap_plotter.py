#!/Applications/Anaconda3/anaconda/bin/python3
# Author: Manzini Stefano; stefano.manzini@gmail.com

__version__ = "0.5.120217"

# TODO:
# write another table, reduced with the same data that is plotted
# fix the issue with names in heatmaps (also possibly replacing seaborn..)
# add --limit-to-mirnas functionality (or rather rows?)


"""
This program is meant to read tables produced by Reconciler.

It counts stuff, and produces an heatmap of the miRNAs
passing the tests, with some inclusion criteria.

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
    description = """miRNA heatmap plotter reads Reconciler tables, counts
                  passed tests and draws heatmaps accordingly.
                  """
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
        "-t", 
        "--passed-tests",
        help = """Discards the data of a given miRNA, if it fails to
                pass that many tests at least.
                """,
        metavar = "",
        dest = "tests_threshold")

parser.add_argument(
        "-O",
        "--limit-to-organs",
        help = ("""
            Limits the plot to selected miRNA organs. 
            Selectable organs: aor, liv, duo, jej, ile, wat, bra
            . Quits if none is given. Implicitly selects all organs if
            this flag is not given.
            """
        ),
        metavar = "",
        dest = "organs",
        nargs = "*")

parser.add_argument(
        "-S",
        "--limit-to-samples",
        help = ("Limits the plot to selected lipid samples. "
                "Selectable samples: aor, liv, pla"
                ". Quits if none is given."
                ),
        metavar = "",
        dest = "samples",
        nargs = "*")

parser.add_argument(
        "-p",
        "--preserve-filenames",
        help = """Filenames for statfriendly, mergediets analyses are shortened 
               to keep filenames at decent length. This specifies to refrain 
               from it.""",
        action = "store_true")

parser.add_argument(
        "-a", 
        "--annotate",
        help = ("""
            Writes values into heatmap cells. Off by default. For reasons 
            unknown to me, it requires -x and -y flags, otherwise it won't 
            annotate anything even if selected.
            """
        ),
        action = "store_true")

parser.add_argument(
        "-f", 
        "--fontsize",
        help = """Sets the font size. It affects both axes and annotated 
            values. If not set, Seaborn attempts at choosing one.
            """,
        metavar = "",
        dest = "fontsize")

parser.add_argument(
        "-x", 
        "--x-labels",
        help = ("Shows the x labels names. Off by default."),
        action = "store_true")

parser.add_argument(
        "-y", 
        "--y-labels",
        help = ("Shows the y labels names. Off by default."),
        action = "store_true")

parser.add_argument(
        "-d",
        "--decent-names",
        help = """
            Annotates miRNA with decent names (like miR-111-3p) instead of the
            long blob. This is useful for making readable heatmaps, but should
            be avoided in order to keep track of the accurate miRNA IDs.
        """,
        action = "store_true")

parser.add_argument(
        "-r",
        "--order-names",
        help = """
            Orders miRNA names alphabetically.
        """,
        action = "store_true")

parser.add_argument(
        "-m", 
        "--mode",
        help = """Sets the preferred behavior. Modes: "save" or "plot"
                . "save" saves the plot; "plot" directly 
                outputs to video. Defaults to save. The filename 
                is auto-generated and probably kept at reasonable 
                length.""",
        metavar = "",
        dest = "mode")

parser.add_argument(
        "-o", 
        "--picture-format",
        help = """
            Sets the output picture format. Defaults to "pdf". Accepts 
            every format matplotlib.pyplot is comfortable with (without 
            checking). Just type the extension, like: -o png
            """,
        metavar = "",
        dest = "picformat")

parser.add_argument(
    "-k",
    "--axes-aspect",
    help = """
        Sets another geometry for the heatmap. Best determined by trial-and-error,
        try 0.5 to 1.5 to begin.
    """,
    metavar = "",
    dest = "aspect",
    type = float
)

parser.add_argument(
    "-X",
    "--x-rotation",
    help = """
        Sets the x tick labels rotation, in integer degrees.
        """,
    metavar = "",
    type = int,
    dest = "x_rotation"
)

parser.add_argument(
    "-Y",
    "--y-rotation",
    help = """
        Sets the y tick labels rotation, in integer degrees.
        """,
    metavar = "",
    type = int,
    dest = "y_rotation"
)

parser.add_argument(
    "-T",
    "--transpose-axes",
    help = """
        Switches x and y axes.
        """,
    action = "store_true"
)

parser.add_argument(
    "-R",
    "--robust",
    help = """
        Draws the colorbar by considering quantiles, and not just the linearity
        from minimum to maximum value.
    """,
    action = "store_true"
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


# ==============================
def humanize_mirname(stringlike):
# ==============================
    """
    Takes something like:
    mmu-miR-434-5p_ID=MIM0001421_MI0001526_mmu-mir-434

    Spits something like:
    miR-434-5p
    """

    match = re.search(r"(.+?)_", stringlike)
    return match.group(1)[4:]



# ============================
def add_limited_organs(string):
# ============================
    
    """This adds, at the end of the filename, the limited organ(s).
    """
    
    extension = string[-4:]
    core = string[:-4] + "_O"
    
    for o in organs:
        core = core + o[0]
    
    return core + extension


# ============================
def add_limited_samples(string):
# ============================
    
    """This adds, at the end of the filename, the limited samples(s).
    """
    
    extension = string[-4:]
    core = string[:-4] + "_S"
    
    for s in samples:
        core = core + s[0]
    
    return core + extension

        
# =======================================
# ===//     END OF FUNCTIONS     //======
# =======================================


# ==============================
# =====     imports     ========
# ==============================
print("Importing modules..")

import os

import minibar
from io import StringIO

# this fixes a bug in the MacOSx backend
import matplotlib
matplotlib.use("Agg")

import regex as re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from manzutils import generate_random_key

#import numpy as np
import pickle

# =====================================
# ===//     end of imports     //======
# =====================================

# ==============================
# =====     GLOBALS     ========
# ==============================

#new in 0.5!
if argv.robust:
    robust = True
else:
    robust = False

dataframe = argv.dataframe

# input/output file separator
if argv.outfile_separator is None:
    sep = "\t"
else:
    sep = argv.outfile_separator

if argv.tests_threshold:
    tests_threshold = int(argv.tests_threshold)
# if not set, the test won't be run (automatically included)

legit_organs = ("aor", "liv", "duo", "jej", "ile", "wat", "bra")
if argv.organs:
    organs = argv.organs
    for o in organs:
        if o not in legit_organs:
            print(o + ": Weird organ. "
                  "Selectable organs: aor, liv, duo, jej, ile, wat, bra.")
            sys.exit()
    if len(organs) < 1:
        print("Please specify at least one organ of the miRNA dataset you "
              "want to keep.")
        sys.exit()
else:
    # all of them
    organs = ("aor", "liv", "duo", "jej", "ile", "wat", "bra")



legit_samples = ("aor", "liv", "pla")
if argv.samples:
    samples = argv.samples
    for s in samples:
        if s not in legit_samples:
            print(s + ": Weird sample. "
                  "Selectable samples: aor, liv, pla.")
            sys.exit()
    if len(samples) < 1:
        print("Please specify at least one sample of the lipid dataset you "
              "want to keep.")
else:
    # all of them
    samples = ("aor", "liv", "pla")

if argv.mode:
    if argv.mode in ("save", "plot"):
        mode = argv.mode
    else:
        print("Unknown mode:", argv.mode)
        sys.exit()
else:
    mode = "save"

# ===============     defining output filenames      ===============

if argv.picformat:
    if len(argv.picformat) > 3:
        sys.exit("Are you sure about picture format?")
    else:
        picformat = argv.picformat
else:
    picformat = "pdf"

# core #
out_filename_prefix = ""
if dataframe.endswith(".csv"):
    out_filename_prefix = dataframe.replace(".csv", "")
elif dataframe.endswith(".txt"):
    out_filename_prefix = dataframe.replace(".txt", "")
else:
    out_filename_prefix = dataframe

# unless specified, we're (probably) trimming the filename
if not argv.preserve_filenames:
    out_filename_prefix = (
        out_filename_prefix
        .replace("statfriendly_mergediets_spearman_", "")
        .replace("allowed0aves", "a0a")
        )

# we're not making a pickled dict with info for every modification
# of the core filename (result of different slices of data)
# we take just this one
pickled_filename = out_filename_prefix
mirna_heatmap_data = (pickled_filename + "_mirheatdata.p")


if argv.organs:
    out_filename_prefix = add_limited_organs(out_filename_prefix)

if argv.samples:
    out_filename_prefix = add_limited_samples(out_filename_prefix)
# // core // #



if argv.tests_threshold:
    out_filename_mirna_heatmap = (
            out_filename_prefix
            + "_mir_heatmap_t"
            + argv.tests_threshold
            + ".csv"
            )
else:
    out_filename_mirna_heatmap = (
            out_filename_prefix
            + "_mir_heatmap.csv"
            )

# the following does not care if the threshold actually exists,
# if it does not, it will not be in the filename
out_thresholded_data_table_filename = out_filename_mirna_heatmap

if argv.organs:
    out_filename_mirna_heatmap = add_limited_organs(out_filename_mirna_heatmap)

if argv.samples:
    out_filename_mirna_heatmap = add_limited_samples(out_filename_mirna_heatmap)


if argv.tests_threshold:
    out_filename_mirna_heatmap_picture = (
            out_filename_prefix
            + "_mir_heatmap_t"
            + argv.tests_threshold
            + "." + picformat
            )
else:
    out_filename_mirna_heatmap_picture = (
            out_filename_prefix
            + "_mir_heatmap"
            + "." + picformat
            )
if argv.organs:
    out_filename_mirna_heatmap_picture = add_limited_organs(
                                            out_filename_mirna_heatmap_picture
                                            )

if argv.samples:
    out_filename_mirna_heatmap_picture = add_limited_samples(
                                            out_filename_mirna_heatmap_picture
                                            )

# =============//     end of output filenames      //=============

# ==============================
# ==== // END OF GLOBALS  // ===
# ==============================

# /*
# Making data for a miRNA-centric heatmap of passed tests #
# */

# making a table, one miRNA per row, counting the occurrences of passed tests
# organs - to - samples; we still have to go through the entire input table

df = pd.read_csv(dataframe, sep)

# loading/creating info #
if os.path.exists(mirna_heatmap_data):
    print("Using previously pickled data for this Reconciler table.")
    mheat = pickle.load(open(mirna_heatmap_data, "rb"))
else:
    mheat = {}
    
    # reading all mirnas from the dataframe, and making an ordered, unique list
    all_mirnas = list(sorted(set(df["mirna"])))
    all_headers = []    # headers of the future table
    for o in legit_organs:  # all organs
        for s in legit_samples: # all samples
            all_headers.append(o + "→" + s)
    
    all_headers = sorted(all_headers)
    
    # populating the mheat dict with the miRNA keys, and every header each key
    for m in all_mirnas:
        mheat.setdefault(m, {})
        for h in all_headers:
            mheat[m].setdefault(h, 0)
    
    print("\nGathering info for miRNA-centric passed tests heatmap..")
    
    for i, mirna in enumerate(minibar.bar(df["mirna"])):
    
        sample_organ_key = df["organ"][i] + "→" + df["sample"][i]
        mheat[mirna][sample_organ_key] += 1
    
    pickle.dump(mheat, open(mirna_heatmap_data, "wb"))
# // done loading/creating info // #


# now we wipe away miRNAs passing fewer tests than wanted:
mirna_keys_to_delete = []
if argv.tests_threshold:
    for m in mheat.keys():
        passed_tests = 0
        for k,v in mheat[m].items():
            passed_tests += v
        if passed_tests < tests_threshold:
            mirna_keys_to_delete.append(m)
    
for d in mirna_keys_to_delete:
    del mheat[d]
# // threshold test done // #  

# now making the table
all_headers = []
for o in legit_organs:
    for s in legit_samples:
        all_headers.append(o + "→" + s)

returnstring = "mirna" + sep

returnstring = returnstring + sep.join(all_headers)
returnstring = returnstring + sep + "delete_me\n"

for m in mheat.keys():
    returnstring = returnstring + m + sep
    for h in all_headers:
        returnstring = returnstring + str(mheat[m][h]) + sep
    returnstring = returnstring + "\n"

# from string to Pandas dataframe.
out_mirna_heatmap = pd.read_csv(
    StringIO(returnstring),
    sep = sep,
    )
del out_mirna_heatmap["delete_me"]

# saving the .csv file | whole (thresholded)  data
out_mirna_heatmap.to_csv(out_thresholded_data_table_filename,
                        sep = sep, 
                        index = False
                        )
print("\n\nData of all mirnas and tests has been saved to:\n{}".format(
    out_thresholded_data_table_filename
    ))


# adjusting the dataframe for plotting #

# we need to pass column 0 as rows index, sns.heatmap will complain otherwise
out_mirna_heatmap = pd.read_csv(
    StringIO(returnstring),
    sep = sep,
    index_col = 0
    )
del out_mirna_heatmap["delete_me"]

# wiping unwanted columns
headers_to_preserve = []
for o in organs:
    for s in samples:
        headers_to_preserve.append(o + "→" + s)

headers_to_delete = [x for x in all_headers if x not in headers_to_preserve]

for h in headers_to_delete:
    del out_mirna_heatmap[h]

# new in 0.5!
# indexes are immutable, but we require editing them if -d flag is turned on
# some trick is required
if argv.decent_names:

    tempfile = "__temp" + generate_random_key(10)

    out_mirna_heatmap.to_csv(tempfile, sep = sep)
    tempdf = pd.read_csv(tempfile, sep = sep)
    tempdf["mirna"] = list(map(humanize_mirname, tempdf["mirna"]))
    tempdf.to_csv(tempfile, sep = sep, index = False)
    out_mirna_heatmap = pd.read_csv(tempfile, sep = sep, index_col = 0)
    os.remove(tempfile)

# new in 0.5!
# indexes are immutable, but we require editing them if -r flag is turned on
# some trick is required
if argv.order_names:

    tempfile = "__temp" + generate_random_key(10)

    out_mirna_heatmap.to_csv(tempfile, sep = sep)
    tempdf = pd.read_csv(tempfile, sep = sep)
    tempdf.sort_values(by = "mirna", inplace = True, axis = 0)
    tempdf.to_csv(tempfile, sep = sep, index = False)
    out_mirna_heatmap = pd.read_csv(tempfile, sep = sep, index_col = 0)
    os.remove(tempfile)


# saving the .csv file | what we actually plot
out_mirna_heatmap.to_csv(out_filename_mirna_heatmap, sep = sep)
print("\n\nData that has been plotted is saved to:\n{}".format(
    out_filename_mirna_heatmap
    ))

# new in 0.4!
if argv.fontsize:
    sns.set(font_scale = float(argv.fontsize))

# new in 0.5!
if argv.transpose_axes:
    out_mirna_heatmap = out_mirna_heatmap.T

# we can now finally plot
heatmap = sns.heatmap(
    out_mirna_heatmap,
    linewidth = 0.5,
    #vmax = 1,  # for lipids - to actually plot something
    xticklabels = argv.x_labels,
    yticklabels  = argv.y_labels,
    annot = argv.annotate,
    fmt = "g",    # with annotate
    cbar = True, # draw the color bar
    robust = robust # new in 0.5!
    )

# new in 0.5!
if argv.x_rotation:
    #heatmap.set_xticklabels(rotation = argv.x_rotation)
    plt.xticks(rotation = argv.x_rotation)

# new in 0.5!
if argv.y_rotation:
    #heatmap.set_yticklabels(rotation = argv.y_rotation)
    plt.yticks(rotation = argv.y_rotation)

# new in 0.5!
if argv.aspect:
    #plt.Axes.set_aspect(out_mirna_heatmap, aspect = argv.aspect)
    axes = plt.gca()
    axes.set_aspect(argv.aspect)

if mode == "plot":
    print("To make this work on OSX, I had to switch backends, since "
          "MacOSX backend has known bugs with the orientation of y-axis "
          "labels. Unfortunately, the working backend, \"agg\", does "
          "not support on-the-fly plots. Writing instead..")
    #sns.plt.savefig(out_filename_mirna_heatmap_picture)
    plt.savefig(out_filename_mirna_heatmap_picture)
    print("\nHeatmap drawn:\n{}".format(out_filename_mirna_heatmap_picture))  
    #plt.show()
else:
    #sns.plt.savefig(out_filename_mirna_heatmap_picture)
    plt.savefig(out_filename_mirna_heatmap_picture)
    print("\nHeatmap drawn:\n{}".format(out_filename_mirna_heatmap_picture))