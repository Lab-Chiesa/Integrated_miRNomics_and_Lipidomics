#!bin/bash

echo
echo =========================================================================
echo This simple script performs all steps needed to perform all correlations
echo    in both miRNA and lipid datasets and reproduce the results of the
echo    paper. Please note that, since the picking of values in the lipid
echo   dataframe and their pairing with miRNA values can vary from analysis
echo     to analysis (see the paper), actual results might slightly vary.
echo =========================================================================
echo
echo Please note that Python 3.6+ is required to run Python programs.
echo There are also required dependencies to be installed, see each
echo file requirements and/or see ModuleNotFoundError complaints.
echo If in trouble, try: less /lib/manzutils.py
echo

export PYTHONPATH=$PYTHONPATH:./lib/

# correlating miRNAs with lipids (default options, but merging diets)
python ./reconciler/reconciler.py data/miRNAs_dataset.csv data/lipidomics_dataset.csv -d

# changing the name of resulting csv table into something more apt
mv statfriendly_mergediets_spearman_mirnas_anon_lipids_anon_0.01_0.7_SD_150p_allowed0aves_1.csv correlations.csv
echo The correlations can be found in: correlations.csv

# after calculating all correlations, we gather data in a final report.
python final_reporter.py correlations.csv
echo "Open correlations_miR_enrich.csv to see the number of miRNAs passing tests."

python mirna_heatmap_plotter.py correlations.csv -x -y -a -r -t 350 -f 0.7 -R
echo See how tests cluster by opening correlations_mir_heatmap_t350.pdf