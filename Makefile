R_pheno_info=${R_DIR}/pheno_info.R
R_functions=${R_DIR}/functions.R
R_generate_prs=${R_DIR}/generate_prs.R
R_generate_prs_manually=${R_DIR}/generate_prs_manually.R
R_analyze_prs=${R_DIR}/analyze_prs.R
R_intersect_snps=${R_DIR}/intersect_snps.R
R_summarize_pca=${R_DIR}/summarize_pca.R
R_simulate=${R_DIR}/simulate_develop.R
R_reformat_scz_2018=${R_DIR}/reformat_scz_2018.R

PRSice=${BIN_DIR}/PRSice_mac
PRSice_R=${BIN_DIR}/PRSice.R
PLINK=$(BIN_DIR)/plink

IBBC_FREQ=${RESULTS_DIR}/${IBBC_PREFIX}.frqx
SNPLIST_INTERSECT=${RESULTS_DIR}/snplist.intersected.txt
SNPLIST_IBBC=${RESULTS_DIR}/snplist.ibbc.txt
LDPRUNED_SNPLIST_PREFIX=${RESULTS_DIR}/snplist.intersected.LDpruned

IBBC_GOOD_SNPS_PREFIX=${RESULTS_DIR}/${IBBC_PREFIX}.good_snps

PCA_PREFIX=${RESULTS_DIR}/pca
PCA1_EIGENVEC=$(PCA_PREFIX)1.eigenvec
PCA1_OUTLIERS=$(PCA_PREFIX)1.outliers
PCA2_EIGENVEC=$(PCA_PREFIX)2.eigenvec

PCA_ALLSNPS_PREFIX=${RESULTS_DIR}/pca.allsnps
PCA_ALLSNPS_EIGENVEC=$(PCA_ALLSNPS_PREFIX).eigenvec
PCA_ALLSNPS_OUTLIERS=$(PCA_ALLSNPS_PREFIX).outliers

AIM2A_HIST_FILENAME=${RESULTS_DIR}/aim2A.hist.pdf
AIM2A_MAIN_FILENAME=${RESULTS_DIR}/aim2A.main.pdf
AIM2A_R2_FILENAME=${RESULTS_DIR}/aim2A.r2.pdf
AIM2B_HIST_FILENAME=${RESULTS_DIR}/aim2B.hist.pdf
AIM2B_MAIN_FILENAME=${RESULTS_DIR}/aim2B.main.pdf
AIM2A_SIM_FILENAME=${RESULTS_DIR}/aimab.power.pdf
AIM2B_SIM_FILENAME=${RESULTS_DIR}/aim2b.power.pdf
AIM2X_TABLE_SCZ_FILENAME=${RESULTS_DIR}/aim2x.tableSCZ.txt
AIM2X_TABLE_IQ_FILENAME=${RESULTS_DIR}/aim2x.tableIQ.txt

all1:install $(PHENO_FILE_WITH_PCA_AND_PRS) $(PCA_ALLSNPS_OUTLIERS)

all2:$(AIM2A_HIST_FILENAME) $(AIM2A_MAIN_FILENAME) $(AIM2B_HIST_FILENAME) $(AIM2B_MAIN_FILENAME) $(AIM2A_SIM_FILENAME) $(AIM2B_SIM_FILENAME)

math.pdf:math.tex
	pdflatex math.tex && pdflatex math.tex

clean:
	set -e && \
	if [ "$(BIN_DIR)" != "" ]; then cd $(BIN_DIR); echo cleaning bin; rm *; fi && \
	if [ "$(RESULTS_DIR)" != "" ]; then cd $(RESULTS_DIR); echo cleaning results; rm *; fi && \
	echo done clean

install:$(PRSice) $(PLINK)

$(PRSice):
	mkdir -p "$(BIN_DIR)" && cd "$(BIN_DIR)" && \
	curl -L -O https://github.com/choishingwan/PRSice/releases/download/2.1.2.beta/PRSice_mac.zip && \
	unzip PRSice_mac.zip

$(PLINK):
	mkdir -p $(BIN_DIR) && cd $(BIN_DIR) && \
	curl -L -O https://www.cog-genomics.org/static/bin/plink180612/plink_mac.zip && \
	unzip plink_mac.zip

## see ${R_dir}/pheno_info.R
NEALE_FLUID_INTELLIGENCE_SUMSTATS=${EXTERNAL_DIR}/fluid_intelligence.neale.tsv.gz
$(NEALE_FLUID_INTELLIGENCE_SUMSTATS):
	R -f ${R_DIR}/download_and_reformat_neale_fluid_intelligence.R

## must be downloaded through a browser
SCZ_SUMSTATS=${EXTERNAL_DIR}/ckqny.scz2snpres.gz

SCZ_2018_SUMSTATS_ORIGINAL=${EXTERNAL_DIR}/clozuk_pgc2.meta.sumstats.txt.gz

SCZ_2018_SUMSTATS=${EXTERNAL_DIR}/clozuk_pgc2.meta.sumstats.reformatted.txt.gz

## second one is to get variants.tsv with rsids in them
$(SCZ_2018_SUMSTATS):$(SCZ_2018_SUMSTATS_ORIGINAL) $(NEALE_FLUID_INTELLIGENCE_SUMSTATS)
	R -f $(R_reformat_scz_2018) --args \
	$(SCZ_2018_SUMSTATS_ORIGINAL) \
	$(SCZ_2018_SUMSTATS)

## https://www.dropbox.com/s/ibjoh0g5e3sdd8t/GWAS_CP_all.txt
LEE_CP_SUMSTATS=${EXTERNAL_DIR}/GWAS_CP_ALL.txt.gz

GWAS_DATA_DOWNLOADED=${EXTERNAL_DIR}/gwas_data_downloaded.txt

$(GWAS_DATA_DOWNLOADED):$(SCZ_SUMSTATS) $(SCZ_2018_SUMSTATS) $(NEALE_FLUID_INTELLIGENCE_SUMSTATS)
	touch $(GWAS_DATA_DOWNLOADED)

$(IBBC_FREQ):$(PLINK)
	mkdir -p $(RESULTS_DIR) && \
	echo "calculate allele frequency" && \
	$(PLINK) \
	--bfile ${IBBC_PLINK} \
	--freqx \
	--out ${RESULTS_DIR}/${IBBC_PREFIX}

$(SNPLIST_IBBC):$(IBBC_FREQ) $(GWAS_DATA_DOWNLOADED)
	R -f $(R_intersect_snps) \
	--args \
	$(IBBC_FREQ) \
	$(SNPLIST_INTERSECT) \
	$(SNPLIST_IBBC)

$(SNPLIST_INTERSECT):$(SNPLIST_IBBC)

## LD-prune input dataset
$(LDPRUNED_SNPLIST_PREFIX).prune.in:$(PLINK) $(SNPLIST_INTERSECT)
	mkdir -p $(RESULTS_DIR) && \
	$(PLINK) \
	--extract $(SNPLIST_INTERSECT) \
	--bfile ${IBBC_PLINK} \
	--indep-pairwise 250kb 1 0.1 \
	--out $(LDPRUNED_SNPLIST_PREFIX)

## PCA1
$(PCA1_EIGENVEC):$(PLINK) $(LDPRUNED_SNPLIST_PREFIX).prune.in
	echo "run PCA" && \
	$(PLINK) \
	--extract $(LDPRUNED_SNPLIST_PREFIX).prune.in \
	--bfile ${IBBC_PLINK} \
	--out $(PCA_PREFIX)1 \
	--pca header

## honestly, these look fine?
$(PCA1_OUTLIERS):$(PCA1_EIGENVEC)
	R -f $(R_summarize_pca) \
	--args \
	$(PCA1_EIGENVEC) \
	$(PCA1_OUTLIERS)

## all-SNPs from IBBC that are MAF > 0.10 and not in MHC and not in 22Q
$(PCA_ALLSNPS_EIGENVEC):$(PLINK) $(SNPLIST_IBBC)
	$(PLINK) \
	--extract $(SNPLIST_IBBC) \
	--bfile ${IBBC_PLINK} \
	--out $(PCA_ALLSNPS_PREFIX) \
	--pca header

temp:$(PCA_ALLSNPS_EIGENVEC)

## honestly, these look fine?
$(PCA_ALLSNPS_OUTLIERS):$(PCA_ALLSNPS_EIGENVEC)
	R -f $(R_summarize_pca) \
	--args \
	$(PCA_ALLSNPS_EIGENVEC) \
	$(PCA_ALLSNPS_OUTLIERS)

$(IBBC_GOOD_SNPS_PREFIX).bed:$(PLINK) $(SNPLIST_IBBC)
	mkdir -p $(RESULTS_DIR) && \
	$(PLINK) \
	--extract $(SNPLIST_IBBC) \
	--bfile ${IBBC_PLINK} \
	--out ${RESULTS_DIR}/${IBBC_PREFIX}.good_snps --make-bed 

## now, build!
## $(R_generate_prs) $(R_functions) 
$(PHENO_FILE_WITH_PCA_AND_PRS):$(PCA_ALLSNPS_EIGENVEC) $(SNPLIST_IBBC) $(IBBC_GOOD_SNPS_PREFIX).bed $(PCA_ALLSNPS_OUTLIERS)
	R -f $(R_generate_prs) \
	--args \
	$(PRSice_R) \
	$(PRSice) \
	$(PCA_ALLSNPS_EIGENVEC) \
	$(SNPLIST_IBBC) \
	$(IBBC_GOOD_SNPS_PREFIX)


## second part of analysis
## simulations and analysis of PRS

## simulations
SIMULATION_COMPLETE=${RESULTS_DIR}/simulations_complete.txt
$(SIMULATION_COMPLETE):
	R -f $(R_simulate) \
	&& touch $@

$(AIM2A_SIM_FILENAME):$(SIMULATION_COMPLETE)

$(AIM2B_SIM_FILENAME):$(SIMULATION_COMPLETE)


PRS_ANALYSIS_COMPLETE=${RESULTS_DIR}/prs_analysis_complete.txt
$(PRS_ANALYSIS_COMPLETE):$(PHENO_FILE_WITH_PCA_AND_PRS)
	R -f $(R_analyze_prs) \
	--args \
	$(PHENO_FILE_WITH_PCA_AND_PRS) \
	$(AIM2A_HIST_FILENAME) \
	$(AIM2A_MAIN_FILENAME) \
	$(AIM2A_R2_FILENAME) \
	$(AIM2B_HIST_FILENAME) \
	$(AIM2B_MAIN_FILENAME) \
	$(AIM2X_TABLE_SCZ_FILENAME) \
	$(AIM2X_TABLE_IQ_FILENAME) \
	${EXTERNAL_DIR}/List_samples_Overlapping_with_CLOZUK.txt
	&& touch $@


all2:$(AIM2A_HIST_FILENAME) $(AIM2A_MAIN_FILENAME) $(AIM2B_HIST_FILENAME) $(AIM2B_MAIN_FILENAME) $(AIM2A_SIM_FILENAME) $(AIM2B_SIM_FILENAME)

## not quite right
$(AIM2A_HIST_FILENAME):$(PRS_ANALYSIS_COMPLETE)

$(AIM2A_MAIN_FILENAME):$(PRS_ANALYSIS_COMPLETE)

$(AIM2A_R2_FILENAME):$(PRS_ANALYSIS_COMPLETE)

$(AIM2B_HIST_FILENAME):$(PRS_ANALYSIS_COMPLETE)

$(AIM2B_MAIN_FILENAME):$(PRS_ANALYSIS_COMPLETE)

