IBBC Aim 2 22Q11DS
==================

Bespoke analysis code to support "Aim 2" of the IBBC i.e. investigating the relationship between schizophrenia and schizophrenia related phenotypes of subthreshold psychosis, low baseline intellectual function and cognitive decline. Additionally supports making and investigating relationships between polygenic risk scores and both schizophrenia and intellectual ability.

Ideally everything is well contained in this repository and can be run using code like the following. Note that very specific input files are required. Please email Robert Davies for details or concerns robertwilliamdavies@gmail.com.

```
export ANALYSIS_VERSION_DATE="2018_11_28"
export PROJ_DIR="~/Dropbox/22Q11/" ## where this repo lives
export STAGE1_DIR=~/IBBC/
export STAGE2_DIR=~/Dropbox/IBBC\ Aim\ II\ group,\ Toronto\ division/analyses/

## generate PRS, PCA (Robbie laptop)
(set -e && . activate && make --file "${PROJ_DIR}/Makefile" -C "${STAGE1_DIR}" all1)

## setup joint Dropbox directory
(set -e && . activate && mkdir -p "${STAGE2_DIR}/${ANALYSIS_VERSION_DATE}" && rsync -av "${STAGE1_DIR}/${PHENO_FILE_WITH_PCA_AND_PRS}" "${STAGE2_DIR}/${ANALYSIS_VERSION_DATE}")
(set -e && . activate && make --file "${PROJ_DIR}/Makefile" -C "${STAGE2_DIR}" "${PHENO_FILE_WITH_PCA_AND_PRS}" -t && touch "${STAGE2_DIR}/${PHENO_FILE_WITH_PCA_AND_PRS}")

## run analysis
(set -e && . activate && make --file "${PROJ_DIR}/Makefile" -C "${STAGE2_DIR}" all2  --dry-run)
```

math.pdf does not require -C

This is partly active as this is a work in progress. Notably external data acquisition is not in Makefile and is in `scripts/external.sh` or is not fully documented yet. 

```
To copy to shared Dropbox for Toronto analysis run
rsync -av /Users/robert\ davies/Dropbox/22Q11/* ~/Dropbox/IBBC\ Aim\ II\ group,\ Toronto\ division/robbie/code/
mkdir ~/Dropbox/IBBC\ Aim\ II\ group,\ Toronto\ division/robbie/2018_07_05
rsync -av ~/IBBC/2018_07_05/iBBC_AIMIIdata_14June2018.withPRS.csv ~/Dropbox/IBBC\ Aim\ II\ group,\ Toronto\ division/robbie/2018_07_05
rsync -av ~/IBBC/2018_07_05/*pdf ~/Dropbox/IBBC\ Aim\ II\ group,\ Toronto\ division/robbie/2018_07_05/
```