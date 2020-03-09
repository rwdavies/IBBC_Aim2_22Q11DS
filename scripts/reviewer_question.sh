## look at the sensitivity of the simulations with respect to
## h2_g for subthreshold psychosis
## h2_g for viq decline

set -e
R -f ~/proj/IBBC_Aim2_22Q11DS/R/simulate_develop.R

for sub in 0.33 0.38 0.43
do
    for viq in 0.15 0.20 0.25
    do
	R -f ~/proj/IBBC_Aim2_22Q11DS/R/simulate_develop.R --args ${sub} ${viq}
    done
done

## can now group together!
