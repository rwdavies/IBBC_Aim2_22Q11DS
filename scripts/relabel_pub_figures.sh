set -e

## Power figures
rsync -av florence:~/22Q/* ~/IBBC/2018_11_28/

home="/Users/rwdavies/"
home="/Users/robert davies/"
from=${home}"IBBC/2018_11_28/"
dest=${home}"Dropbox/IBBC Aim II group/Aim II paper/AIM II figures/1. Current AIM II figures/"

## PCA plot
cp "${from}pc1.fancy.pdf" "${dest}Suppl Figure 1.pdf"
exit

## Power is, currently, supp figures 3, 4, 7, supp table 2
cp "${from}aim2a.power.bothPS.pdf" "${dest}Suppl Figure 3.pdf"
cp "${from}aim2b.power.only4.pdf" "${dest}Suppl Figure 4.pdf"
cp "${from}aim2ab.quantSIPS.pdf" "${dest}Suppl Figure 7.pdf"
cp "${from}power.simple.csv" "${dest}Suppl Table 2.csv"

## Figure 1
cp "${from}aim2A.main.pdf"  "${dest}Main Figure 1.pdf"

## Table 1
cp "${from}sample_characteristics.csv"  "${dest}Table 1.csv"


## 
