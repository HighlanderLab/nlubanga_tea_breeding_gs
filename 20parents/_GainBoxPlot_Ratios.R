# Run after PlotProgress_finalPlots

#########################################################################
# Perform pairwise stat tests using Tukey Test and make a boxplot 
#########################################################################
library(multcompView)
test <- rawData[rawData$year == 40,]
str(test)

# Analysis of variance
summary(anova <- aov(meanSeed ~ Scenario, data = test))

# Tukey's test
(tukey <- TukeyHSD(anova,conf.level = 0.05))

# Compact letter display
(cld <- multcompLetters4(anova, tukey))

# table with factors and 3rd quantile
Tk <- ddply(test,c("year","scenario","Scenario"), summarize,
      mean = mean(meanSeed),
      quant = quantile(meanSeed, probs = 0.75))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$Scenario)
Tk$cld <- cld$Letters[match(Tk$Scenario,rownames(cld))]
Tk$cld <- c("a","b","c","c","c")

a1 <- ggplot(rawData[rawData$year == 40,], aes(x=scenario, y=meanSeed,color=scenario)) +
  geom_boxplot() +
  stat_summary(fun ="mean", color="gray", shape=1) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  geom_hline(yintercept = mean(rawData$meanSeed[rawData$year == 40 & rawData$scenario == "Seed-GSc"]), linetype = 2) +
  scale_colour_manual(values=c("black",wes_palette("FantasticFox1")[2:5])) +
  geom_text(data = Tk, aes(x = scenario, y = quant, label = cld),color="black", size = 5, vjust=-1, hjust =-2) +
  # ggtitle(title) +
  scale_y_continuous("Genetic mean",limits = c(4000,18000)) +
  scale_x_discrete("Program") +
  optns +
  theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1))
ggsave(filename = paste("PlotProgress_finalPlots/Gain_boxplot_Seed.png"), plot = a1, width = 6, height = 5, scale = 1.1)
 
#########################################################################
# Calculate ratios between PS and othe programs
#########################################################################
# Calculates the ratio between two variables
# Returns a 95% confidence-interval for the ratio
# Use paired=TRUE if values come from a common burn-in
calcRatio = function(x,y,paired=FALSE,silent=FALSE){
  if(any(x<=0) | any(y<=0)){
    stop("x and y values must be greater than 0")
  }
  x = log(x)
  y = log(y)
  model = t.test(x,y,paired=paired)
  if(paired){
    ratio = exp(model$estimate)
  }else{
    ratio = exp(model$estimate[1]-model$estimate[2])
  }
  ratio = round(unname(ratio),2)
  l95 = round(exp(model$conf.int[1]),2)
  u95 = round(exp(model$conf.int[2]),2)
  if(!silent){
    cat(round(ratio,2)," (95% CI [",l95,", ",u95,"])", sep = "")
  }
  output = list(ratio=ratio,CI=c(l95,u95))
  invisible(output)
}

# Comparison with PS program
for (i in levels(test$scenario)[-1]) {
  cat("PS vs.",i,"   ")
  calcRatio(x = as.matrix(test[test$scenario == i,"meanSeed"]), 
            y = as.matrix(test[test$scenario == "PS","meanSeed"]), 
            paired = T) 
  cat("\n")
}
# PS vs. Seed-Ped    1.25 (95% CI [1.21, 1.29])
# PS vs. Seed-GSc    1.71 (95% CI [1.66, 1.76])
# PS vs. Seed-GSunc  1.75 (95% CI [1.7, 1.81])
# PS vs. ECT-GS      1.64 (95% CI [1.59, 1.68])
# (1.71+1.75+1.64)/3

# Comparison with Seed-Ped program
for (i in levels(test$scenario)[3:5]) {
  cat("Seed-Ped vs.",i,"   ")
  calcRatio(x = as.matrix(test[test$scenario == i,"meanSeed"]), 
            y = as.matrix(test[test$scenario == "Seed-Ped","meanSeed"]), 
            paired = T) 
  cat("\n")
}
# Seed-Ped vs. Seed-GSc     1.37 (95% CI [1.34, 1.4])
# Seed-Ped vs. Seed-GSunc   1.40 (95% CI [1.37, 1.44])
# Seed-Ped vs. ECT-GS       1.31 (95% CI [1.28, 1.34])
# (1.37+1.40+1.31)/3
 
save.image("PlotProgress_finalPlots/Data100.RData")
