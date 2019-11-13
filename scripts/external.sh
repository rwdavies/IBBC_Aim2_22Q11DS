ANALYSIS_DIR="/hpf/largeprojects/tcagstor/scratch/rwdavies/iBBC_22Q11/"

## setup using external data
mkdir -p ${ANALYSIS_DIR}
cd ${ANALYSIS_DIR}
mkdir external
cd external

## manually rsync from client (laptop) to ${ANALYSIS_DIR}/external/ the file
## 22q_IBBC_forRobbie_Imputed_QCed.zip
## check md5sum here
md5sum 22q_IBBC_forRobbie_Imputed_QCed.zip ## should be 5dbac05506681303c5c20117ececea23
unzip 22q_IBBC_forRobbie_Imputed_QCed.zip

## PGC - transferred over - md5sum
390004e6bd7f89300e3bbacbe5807dcd  iPSYCH-PGC_ASD_Nov2017 

## 
