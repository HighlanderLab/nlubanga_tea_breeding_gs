#!/bin/bash
## This bash script runs all scenarios at once

# no spaces when giving commands
SCENARIOS='02_PS5 03_PedigreeBlup_Seed5 04_Seed_GSconst_pheno_noECT5 05_Seed_GSunconst_pheno_noECT5 06_ECT_GS5'

# run scenarios 
for i in $SCENARIOS;
do
  	qsub ${i}.R
done

echo Jobs submitted
