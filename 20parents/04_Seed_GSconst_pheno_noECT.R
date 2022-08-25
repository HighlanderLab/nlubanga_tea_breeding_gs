#!/exports/cmvm/eddie/eb/groups/HighlanderLab/communal/R-4.1.3/R-4.1.3/bin/Rscript
#$ -N SdCon_phEC
#$ -cwd
#$ -R y
#$ -pe sharedmem 1
#$ -l h_vmem=50G
#$ -l h_rt=010:00:00
#$ -j y
#$ -V
#$ -P roslin_HighlanderLab
#$ -M jbancic@exseed.ed.ac.uk

rm(list = ls())
library(AlphaSimR)

## Cluster parameters
REP = Sys.getenv("SGE_TASK_ID")

# Load data
load(paste0("BURNIN_", REP, ".RData"))

# Set alternative scenario name
scenarioName = "Seed_GSconst_pheno_noECT"

# Change scenario
output$scenario = scenarioName

# Run 20 cycle of breeding program
for(year in (nBurnIn + 1):nCycles){ 
  cat(" Working on year ", year, "\n", sep = "")
  # Run GS model
  gsmodel = RRBLUP(pop=trainPop, useReps = TRUE)
  
  # Select new parents in current year before advancing the material
  ParentPool = Seedlings   
  ParentPool = setEBV(ParentPool,gsmodel,value = "bv")
  Parents = selectInd(ParentPool, nInd = 20, use = "ebv")
  
  # Year 13
  ECT6 = setPheno(ECT5, varE = VarE, reps = repECT, p = P[year])
  
  # Year 12
  ECT5 = ECT4
  
  # Year 11
  ECT4 = ECT3
  
  # Year 10
  ECT3 = ECT2
  
  # Year 9
  ECT2 = ECT1
  
  # Year 8
  output$accACT[year] = cor(gv(ACT5), pheno(ACT5)) # accuracy based on 500 inds
  ECT1 = selectInd(ACT5, nInd = nClonesPhenoECT, use = "pheno")
  
  # Year 7
  ACT5 = setPheno(ACT4, varE = VarE, reps = repACT, p = P[year])
  
  # Year 6
  ACT4 = ACT3
  
  # Year 5
  ACT3 = ACT2
  
  # Year 4
  ACT2 = ACT1
  
  # Year 3
  ACT1 = setEBV(Seedlings, gsmodel, value = "bv")
  output$accSeed[year] = cor(gv(ACT1), ebv(ACT1)) # accuracy based on 800 inds
  ACT1 = selectInd(ACT1, nInd = ncloneSeedGScostACT, use = "ebv")
  
  # Year 2
  Seedlings = setPheno(F1, reps = repHPT, p = P[year])
  
  # Year 1
  # Crossing block
  F1 = randCross(Parents, nCrosses=ncross, nProgeny = nProgSeed_GSconst)
  
  # Update training population
  ACT5@fixEff <- as.integer(rep(year,nInd(ACT5))) #specify year
  trainPop = c(trainPop,ACT5)
  
  # Report mean and variance
  output$meanParent[year] = meanG(Parents)
  output$varParent[year] = varG(Parents)
  
  output$meanSeed[year] = meanG(Seedlings)
  output$meanACT[year] = meanG(ACT5)
  output$meanECT[year] = meanG(ECT6)
  
  output$varSeed[year] = varG(Seedlings)
  output$varACT[year] = varG(ACT5)
  output$varECT[year] = varG(ECT6)
  
  gp = genParam(Seedlings)
  output$GenicVar[year] = genicVarG(Seedlings)  
  output$CovG_HW[year] = gp$covG_HW        #Genetic Covariance nonHWE
  
}

# Save output
cat(paste0("saving REP ", REP, "\n"))
saveRDS(output, paste0(scenarioName, "_", REP, ".rds"))  

