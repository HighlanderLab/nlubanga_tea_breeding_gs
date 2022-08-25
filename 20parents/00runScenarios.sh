#!/bin/bash
## This bash script runs all scenarios at once

# no spaces when giving commands
SCENARIOS='02_PS 03_PedigreeBlup_Seed 04_Seed_GSconst_pheno_noECT 05_Seed_GSunconst_pheno_noECT 06_ECT_GS'

# run scenarios 
for i in $SCENARIOS;
do
  	qsub ${i}.R
done

echo Jobs submitted
