#!bin/bash

# This script 

# correlating miRNAs with lipids (default options, but merging diets)
python reconciler.py mirnas_anon.csv lipids_anon.csv -d

# changing the name of resulting csv table into something more apt
mv statfriendly_mergediets_spearman_mirnas_anon_lipids_anon_0.01_0.7_SD_150p_allowed0aves_1.csv correlations.csv
echo The correlations can be found in: correlations.csv

# after calculating all correlations, we gather data in a final report.
python final_reporter.py correlations.csv
echo "Open correlations_miR_enrich.csv to see the number of miRNAs passing tests."

python mirna_heatmap_plotter.py correlations.csv -x -y -a -r -t 350 -f 0.7 -R
echo See how tests cluster by opening correlations_mir_heatmap_t350.pdf
