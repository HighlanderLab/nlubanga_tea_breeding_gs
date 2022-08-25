rm(list=ls())
library(tidyverse)
library(dplyr)
library(ggplot2)
library(coda)
library(gridExtra)
library(reshape2)
library(ggpubr)
library(wesanderson)

nReps = 30
scenarios = c("PS","PedigreeBlup_Seed",
              "Seed_GSconst_pheno_noECT",
              "Seed_GSunconst_pheno_noECT",
              "ECT_GS")
title = "Replace 100% parents per year" # plot title

# Check if any files missing
temp <- paste(rep(scenarios,each = nReps),"_", rep(1:nReps, times = length(scenarios)), ".rds",sep = "")
temp[!file.exists(temp)]

# Read in output data
rawData = vector("list",nReps*length(scenarios))
i = 0L
for(SCENARIO in scenarios){
  for(REP in 1:nReps){
    i = i+1L
    FILE = paste0(SCENARIO,"_",REP,".rds")
    temp = readRDS(FILE)
    
    # Reset year to years since burn-in
    temp$year = temp$year-40
    temp$scenario = SCENARIO
    rawData[[i]] = temp
    
  }
}
rawData = bind_rows(rawData)

# Replace Scenario names with prettier names
rawData$Scenario = rawData$scenario
rawData$scenario = sub("^PS$","PS",rawData$scenario)
rawData$scenario = sub("^PedigreeBlup_Seed$","Pedigree",rawData$scenario)
# rawData$scenario = sub("^PedigreeBlup$","Pedigree",rawData$scenario)
rawData$scenario = sub("^Seed_GSconst_pheno_noECT$","Seed-GSc",rawData$scenario)
rawData$scenario = sub("^Seed_GSunconst_pheno_noECT$","Seed-GSunc",rawData$scenario)
rawData$scenario = sub("^ECT_GS$","ECT-GS",rawData$scenario)

# Re-order the scenarios for the legend
rawData$scenario <- factor(rawData$scenario, 
                           levels=c("PS", 
                                    "Pedigree",
                                    "Seed-GSc",
                                    "Seed-GSunc",
                                    "ECT-GS"))
levels(as.factor(rawData$scenario))

dat = rawData

initial = dat %>%
  group_by(rep) %>% # each replicate separately
  filter(year == min(year), scenario == "PS") %>% # get the initial values
  dplyr::select(rep, meanSeed, varSeed, GenicVar) # these are the vars we need
colnames(initial) = c("rep", "initialMean", "initialVarG", "initialVarGenic")
dat = full_join(dat, initial, by="rep")

# Genetic mean set to zero in first generation and in units of first generation genetic variance
dat$zMean      = (dat$meanSeed - dat$initialMean) / sqrt(dat$initialVarG)
dat$zMeanGenic = (dat$meanSeed - dat$initialMean) / sqrt(dat$initialVarGenic)

# Variance set relative to initial variance
dat$zVarG     = dat$varSeed     / dat$initialVarG
dat$zVarGenic = dat$GenicVar / dat$initialVarGenic

# Standard deviations
dat$sdG      = sqrt(dat$varSeed)
dat$sdGenic  = sqrt(dat$GenicVar)
dat$zSdG     = sqrt(dat$zVarG)
dat$zSdGenic = sqrt(dat$zVarGenic)


dataArrow = dat %>%
  group_by(scenario, rep) %>%
  dplyr::mutate(zSdGeneticMin = min(zSdGenic), zSdGeneticMax = max(zSdGenic)) %>% 
  #dplyr::mutate(zMeanSlopeSdGenic = mean(coef(lm(zMeanT1 ~ 1 + I(1-zSdGenic_T1) ))[2])) %>%
  #dplyr::mutate(zMeanInterceptSdGenic = mean(coef(lm(zMeanT1 ~ 1 + I(1-zSdGenic_T1) ))[1])) %>%
  dplyr::mutate(zMeanSlopeSdGenic = mean(coef(lm(zMean ~ 1 + zSdGenic ))[2])) %>%
  dplyr::mutate(zMeanInterceptSdGenic = mean(coef(lm(zMean ~ 1 + zSdGenic ))[1])) %>%
  summarize(zMeanInterceptSdGenic = mean(zMeanInterceptSdGenic, na.rm = TRUE),
            zMeanSlopeSdGenic     = mean(zMeanSlopeSdGenic,     na.rm = TRUE),
            zSdGeneticMin = mean(zSdGeneticMin, na.rm = TRUE),
            zSdGeneticMax = mean(zSdGeneticMax, na.rm = TRUE)) %>%
  dplyr::mutate(zMeanStartSdGenetic = zMeanInterceptSdGenic + zMeanSlopeSdGenic*zSdGeneticMax, 
                zMeanEndSdGenetic   = zMeanInterceptSdGenic + zMeanSlopeSdGenic*zSdGeneticMin) %>%
  group_by(scenario) %>%
  summarise_all(mean, na.rm = TRUE) 

color_roslin=c("black",wes_palette("FantasticFox1")[2:5])

optns <- theme_bw(base_size = 16, base_family = "sans") +
  theme(panel.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 15,colour = "black"),
        axis.title = element_text(size = 16),
        strip.text = element_text(face = "bold", size = 20, colour = 'blue'))

Plot1 = ggplot(dat, aes(x = zSdGenic, y = zMean, colour = scenario, group = interaction(scenario))) + 
  theme_classic() +  geom_path(size = 2, alpha = 0.1) + 
  scale_x_reverse(sec.axis = sec_axis(trans=~1-., name = "Converted/Lost genic standard deviation")) +  
  xlab("Genic standard deviation") + ylab("Genetic mean") +
  geom_segment(data = dataArrow, mapping = aes(x = zSdGeneticMax, xend = zSdGeneticMin,
                                               y = zMeanStartSdGenetic, yend = zMeanEndSdGenetic,
                                             colour = scenario), arrow = arrow(), show.legend = TRUE) +
  scale_color_manual(values = c(color_roslin, color_roslin, color_roslin))+optns

ggsave(filename = paste("PlotProgress_finalPlots/Efficiency_plot.png"), plot = Plot1, width = 5, height = 5, scale = 1.1)

