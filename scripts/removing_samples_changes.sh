# This file manually re-runs all code to re-generate main figures, and supporting files to edit main text,
#   when needing to do this to remove a handful of subjects

# Note, this is not written all that well, and uses hard coded links within scripts 
# To change, either modify clozuk_overlap external file
# or grep "clozuk_overlap" as the R object with the list of samples to be removed
# and add additional links nearby as appropriate

# run with bash "/Users/robert davies/proj/IBBC_Aim2_22Q11DS/removing_samples_changes.sh"
set -e 
# If removing samples (e.g. 3 clozuk samples)
# can execute / work interactively
OUT_DIR="/Users/robert davies/IBBC/plots_2020_03_10/"
mkdir -p "${OUT_DIR}"



# Items:  Main figures, main tables
# How:    Automatic
# What:   Many
# File:   R/anlyze_prs.R
# manual switch on
export MANUAL_ANALYSIS_SWITCH=TRUE # to avoid running whole makefile
R -f "/Users/robert davies/proj/IBBC_Aim2_22Q11DS/R/analyze_prs.R"
# Notes:
# For Table 1
rsync -av "/Users/robert davies/IBBC/2018_11_28/sample_characteristics.csv" "${OUT_DIR}/Table 1.csv"
# that can be swapped into Table 1, propagating values into the text appropriately
# For Table 2
# Not in a single file, go through below for figure 1 and figure 2, but with "withbeta" versions
# Not clear there is a single file? Read off of output plots
# Figure 1
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2B.main.old.pdf" "${OUT_DIR}/Figure1.pdf"
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2B.main.old.withbeta.pdf" "${OUT_DIR}/for_table2.old.pdf"
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2B.main.scz.withbeta.pdf.definite.withbeta.pdf" "${OUT_DIR}/for_text_definite_see_bottom_left.pdf"
# Figure 2
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2A.main.pdf" "${OUT_DIR}/Figure2.pdf"
# Figure 3
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2B.main.new.pdf" "${OUT_DIR}/Figure3.pdf"
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2B.main.new.withbeta.pdf" "${OUT_DIR}/for_table2.new.pdf"
# Figure 4
rsync -av "/Users/robert davies/IBBC/2018_11_28/cutoffPSile_both.pdf" "${OUT_DIR}/Figure4.pdf"
# Supplemental Table 1
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2X.tableSCZ.csv" "${OUT_DIR}/Supplemental Table 1.csv"
# Supplemental Table 2
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2X.tableIQ.csv" "${OUT_DIR}/Supplemental Table 2.csv"
# Supplemental Table 3 - mediation analysis
rsync -av "/Users/robert davies/IBBC/2018_11_28/mediation.csv" "${OUT_DIR}/Supplemental Table 3.csv"
# Supplemental Figure 2
rsync -av "/Users/robert davies/IBBC/2018_11_28/cont.frac.ps_scz.pdf" "${OUT_DIR}/Supplemental Figure 2.pdf"


# Item:   Supplementary Figure 4, 5
# How:    Manual
# What:   Flowchart of study participants
# Notes:  On Jacob


# Items:  Quant subthreshold 
# How:    Semi-automatic
R -f "/Users/robert davies/proj/IBBC_Aim2_22Q11DS/R/examine_quant_subthreshold.R"
rsync -av "/Users/robert davies/IBBC/2018_11_28/sips1.ps_scz.pdf" "${OUT_DIR}/Supplemental Figure 3.pdf"
rsync -av "/Users/robert davies/IBBC/2018_11_28/quant_sub.ori_only.pdf" "${OUT_DIR}/Supplemental Figure 8.pdf"




# Items:  PCA
# How:    Semi-automatic
# File:   R/summarize_pca.R
export MANUAL_ANALYSIS_SWITCH=TRUE # to avoid running whole makefile
R -f "/Users/robert davies/proj/IBBC_Aim2_22Q11DS/R/summarize_pca.R"
rsync -av "/Users/robert davies/IBBC/2018_11_28/pc1.fancy.pdf" "${OUT_DIR}/Supplemental Figure 1.pdf"






# Items:  Power
# How:    Semi-automati
# File:   R/simulate_develop.R
# run this file semi-interactively
export MANUAL_ANALYSIS_SWITCH=TRUE # to avoid running whole makefile
# R -f "/Users/robert davies/proj/IBBC_Aim2_22Q11DS/R/simulate_develop.R"
# Supplementary Table 4
rsync -av "/Users/robert davies/IBBC/2018_11_28/parameter.estimates.csv" "${OUT_DIR}/Supplemental Table 4.csv"
# Rest below on server
# Supplementary Table 5 - based on Supp Table 6, 7
# Requires running on smew, after copying over model parameters from laptop
# Supplemental Figure 9 also on smew
if [ 1 == 0 ]
then
    # Supp tables 6, 7
    # manually specify user name
    FROMDIR="rescompNew:/well/davies/users/${USER_NAME}/22Qresults/"
    rsync -av ${FROMDIR}aim2a.power.bothPS.pdf "${OUT_DIR}/Supplemental Figure 6.pdf"
    rsync -av ${FROMDIR}aim2b.power.only4.pdf "${OUT_DIR}/Supplemental Figure 7.pdf"
    rsync -av ${FROMDIR}aim2ab.quantSIPS.pdf "${OUT_DIR}/Supplemental Figure 9.pdf"
    rsync -av ${FROMDIR}power.simple.csv "${OUT_DIR}/Supplemental Table 5.csv"
fi

exit
