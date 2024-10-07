#!/bin/bash

echo
echo "========================================================================="
echo "This simple script performs all steps needed to perform all correlations"
echo "   in both miRNA and lipid datasets and reproduce the results of the"
echo "   paper. Please note that, since the picking of values in the lipid"
echo "  dataframe and their pairing with miRNA values can vary from analysis"
echo "    to analysis (see the paper), actual results might slightly vary."
echo "========================================================================="
echo
echo "REQUIREMENTS"
echo "============"
echo
echo "Please note that Python 3.6+ is required to run supplied Python programs."
echo "There are also required dependencies to be installed, see each"
echo "file requirements and/or see ModuleNotFoundError complaints."
echo
echo "** If in trouble with dependencies, try: less ./lib/manzutils.py **"
echo
echo "MAKING & ANALYZING RANDOMIZED DATASETS"
echo "======================================="
echo
echo "In the paper, one thousand datasets (both miRNAs and lipids) have been"
echo "randomized (see paper), and tested for correlation. This script"
echo "produces one random dataset for both miRNAs and lipids, then performs"
echo "the merged diets analysis on it."
echo "The script can be easily modified to produce as many random dataset"
echo "as required, as well as to perform as many analyses as needed."
echo "Please consider that one full run produces 15-20 MB of data, and takes"
echo "a few hours to complete."
echo
echo "Now running tests. Please note that this might take a few hours."

export PYTHONPATH=$PYTHONPATH:./lib/

# making the random datasets to analyze
python -W ignore ./source/dataset_randomizer.py -m data/miRNAs_dataset.csv -l data/lipidomics_dataset.csv

# correlating miRNAs with lipids (default options, but merging diets)
python ./reconciler/reconciler.py rnd_mirna_0.csv rnd_lipi_0.csv -d

# changing the name of resulting csv table into something more apt
mv statfriendly_mergediets_spearman_rnd_mirna_0_rnd_lipi_0_0.01_0.7_SD_150p_allowed0aves_1.csv merged_diets_random_correlations.csv
echo
echo "** The correlations can be found in: merged_diets_random_correlations.csv **"
echo

# after calculating all correlations, we gather data in a final report.
python ./source/final_reporter_rnd.py merged_diets_random_correlations.csv
echo
echo "**Open correlations_miR_enrich.csv to see the number of miRNAs passing tests.**"
echo

python ./source/mirna_heatmap_plotter.py merged_diets_random_correlations.csv -x -y -r -t 25 -f 0.7 -R
echo
echo "See how tests cluster by opening merged_diets_random_correlations_mir_heatmap_t25.pdf"
echo
echo Done!
