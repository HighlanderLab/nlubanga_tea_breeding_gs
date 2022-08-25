#!/exports/cmvm/eddie/eb/groups/HighlanderLab/communal/R-4.1.3/R-4.1.3/bin/Rscript
#$ -N ECT_GS
#$ -cwd
#$ -R y
#$ -pe sharedmem 1
#$ -l h_vmem=50G
#$ -l h_rt=020:00:00
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
scenarioName = "ECT_GS"

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
  
  # Year 9
  ECT7 = setPheno(ECT6, varE = VarE, reps = repECT, p = P[year])
  
  # Year 8
  ECT6 = ECT5
  
  # Year 7
  ECT5 = ECT4
  
  # Year 6
  ECT4 = ECT3
  
  # Year 5
  ECT3 = ECT2
  
  # Year 4
  ECT2 = ECT1
  
  # Year 3
  ECT1 = setEBV(Seedlings, gsmodel, value = "bv")
  output$accSeed[year] = cor(gv(ECT1), ebv(ECT1)) # accuracy based on 800 inds
  ECT1 = selectInd(ECT1, nInd = ncloneECTGS, use = "ebv")
  
  # Year 2
  Seedlings = setPheno(F1, reps = repHPT, p = P[year])
  
  # Year 1
  # Crossing block
  F1 = randCross(Parents, nCrosses=ncross, nProgeny = nPrognProgECT_GS)
  
  # Update training population
  if (year < 46) {
    ECT6 = setPheno(ECT6, varE = VarE, reps = repECT, p = P[year])
    ECT6@fixEff <- as.integer(rep(year,nInd(ECT6)))
    trainPop = c(trainPop, ECT6) 
  } else if (year == 46) {
    #phenos from inds from year 41 will be added to trainpop one year earlier but with 
    #lower heritability (0.61)
    ECT6 = setPheno(ECT6, varE = VarE, reps = repECT, p = P[year]) #higher VarE
    ECT6@fixEff <- as.integer(rep(year,nInd(ECT6))) #specify year
    trainPop = c(trainPop, ECT6)
  } else if (year == 47) {
    #remove records from previous year and replace them with the right ones
    trainPop <- trainPop[-c((nInd(trainPop)+1-nInd(ECT6)):nInd(trainPop))]
    ECT7@fixEff <- as.integer(rep(year,nInd(ECT7)))
    trainPop = c(trainPop, ECT7) 
  } else {
    ECT7@fixEff <- as.integer(rep(year,nInd(ECT7)))
    trainPop = c(trainPop, ECT7) 
  }
  
  # Report mean and variance
  if (year <= 47) {
    output$meanParent[year] = meanG(Parents)
    output$meanSeed[year] = meanG(Seedlings)
    output$meanECT[year] = meanG(ECT6)
    output$varParent[year] = varG(Parents)
    output$varSeed[year] = varG(Seedlings)
    output$varECT[year] = varG(ECT6)
    gp = genParam(Seedlings)
    output$GenicVar[year] = genicVarG(Seedlings)  
    output$CovG_HW[year] = gp$covG_HW        #Genetic Covariance nonHWE
    
  } else {
    output$meanParent[year] = meanG(Parents)
    output$meanSeed[year] = meanG(Seedlings)
    output$meanECT[year] = meanG(ECT7)
    output$varParent[year] = varG(Parents)
    output$varSeed[year] = varG(Seedlings)
    output$varECT[year] = varG(ECT7)
    gp = genParam(Seedlings)
    output$GenicVar[year] = genicVarG(Seedlings)  
    output$CovG_HW[year] = gp$covG_HW        #Genetic Covariance nonHWE
    
  }
  
}

# Save output
cat(paste0("saving REP ", REP, "\n"))
saveRDS(output, paste0(scenarioName, "_", REP, ".rds"))  

