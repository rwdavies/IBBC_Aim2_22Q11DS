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
OUT_DIR="/Users/robert davies/IBBC/plots_2020_09_02/"
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
# (supplementary figure, formerly Figure 1)
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2B.main.old.pdf" "${OUT_DIR}/prev_Figure1.pdf"
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2B.main.old.withbeta.pdf" "${OUT_DIR}/for_table2.old.pdf"
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2B.main.scz.withbeta.pdf.definite.withbeta.pdf" "${OUT_DIR}/for_text_definite_see_bottom_left.pdf"
# Figure 1 (former Figure 2)
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2A.main.pdf" "${OUT_DIR}/Figure 1.pdf"
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2A.main.eps" "${OUT_DIR}/Figure 1.eps"
# Figure 2 (new Figure 2)
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2B.main.new.special.pdf" "${OUT_DIR}/Figure 2.pdf"
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2B.main.new.special.eps" "${OUT_DIR}/Figure 2.eps"
# EDF 1 (formerly a main text Figure, then at one point supp figure 2)
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2B.main.old.special.pdf" "${OUT_DIR}/EDF_1.pdf"
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2B.main.old.special.eps" "${OUT_DIR}/EDF_1.eps"
# Formerly Figure 3, now related to Figure 2, but new formatting version
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2B.main.new.pdf" "${OUT_DIR}/prev_Figure3.pdf"
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2B.main.new.withbeta.pdf" "${OUT_DIR}/for_table2.new.pdf"
# Figure 3
rsync -av "/Users/robert davies/IBBC/2018_11_28/cutoffPSile_both.pdf" "${OUT_DIR}/Figure 3.pdf"
rsync -av "/Users/robert davies/IBBC/2018_11_28/cutoffPSile_both.eps" "${OUT_DIR}/Figure 3.eps"
# Supplemental Table 1
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2X.tableSCZ.csv" "${OUT_DIR}/Supplemental Table 1.csv"
# Supplemental Table 2
rsync -av "/Users/robert davies/IBBC/2018_11_28/aim2X.tableIQ.csv" "${OUT_DIR}/Supplemental Table 2.csv"
# Supplemental Table 3 - mediation analysis
rsync -av "/Users/robert davies/IBBC/2018_11_28/mediation.csv" "${OUT_DIR}/Supplemental Table 3.csv"

# Now EDF2, Formally Supplemental Figure 2
rsync -av "/Users/robert davies/IBBC/2018_11_28/cont.frac.ps_scz.pdf" "${OUT_DIR}/EDF_2.pdf"
rsync -av "/Users/robert davies/IBBC/2018_11_28/cont.frac.ps_scz.eps" "${OUT_DIR}/EDF_2.eps"

exit

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
    echo get files from rescomp
    # manually specify user name
    USER_NAME=$1
    PASSWORD=$2 ## temp
    FROMDIR="rescompNew2:/well/davies/users/${USER_NAME}/22Qresults/"
    SSHPASS="/Users/robert davies/Downloads/sshpass/sshpass"
    "${SSHPASS}" -p ${PASSWORD} rsync -av ${FROMDIR}aim2a.power.bothPS.pdf "${OUT_DIR}/Supplemental Figure 6.pdf"
    "${SSHPASS}" -p ${PASSWORD} rsync -av ${FROMDIR}aim2b.power.only4.pdf "${OUT_DIR}/Supplemental Figure 7.pdf"
    "${SSHPASS}" -p ${PASSWORD} rsync -av ${FROMDIR}aim2ab.quantSIPS.pdf "${OUT_DIR}/Supplemental Figure 9.pdf"
    "${SSHPASS}" -p ${PASSWORD} rsync -av ${FROMDIR}power.simple.csv "${OUT_DIR}/Supplemental Table 5.csv"
    "${SSHPASS}" -p ${PASSWORD} rsync -av ${FROMDIR}q*csv "${OUT_DIR}/"
fi



exit
