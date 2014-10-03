# Anne S. Warlaumont

# Test to see if there are statistically significant differences between Realistic and Abstract simulations at generation 100

# If you run with the following command, a file containing the output will be saved to the directory containing this script:
# R CMD BATCH AnalyzeFitness.r

rm(list=ls())

library(lme4)
library(lmerTest)

# Change the directory as needed:
#simdata = read.csv('~/Downloads/MultipleRunsData.csv',header=T)
simdata = read.csv('/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/MultipleRunsData.csv',header=T)

simdata$useVocalTract = as.factor(simdata$useVocalTract);
simdata$medProFitnessScaled = scale(simdata$medProFitness)
simdata$medPerFitnessScaled = scale(simdata$medPerFitness)
simdata$generationScaled = scale(simdata$generation)

t.test(subset(subset(simdata,useVocalTract==0),generation==100)$medProFitness,subset(subset(simdata,useVocalTract==1),generation==100)$medProFitness)
t.test(subset(subset(simdata,useVocalTract==0),generation==100)$medPerFitness,subset(subset(simdata,useVocalTract==1),generation==100)$medPerFitness)

t.test(subset(subset(simdata,useVocalTract==0),generation==500)$medProFitness,subset(subset(simdata,useVocalTract==1),generation==500)$medProFitness)
t.test(subset(subset(simdata,useVocalTract==0),generation==500)$medPerFitness,subset(subset(simdata,useVocalTract==1),generation==500)$medPerFitness)

proFitModel = lmer(medProFitnessScaled~(1|runfolder)+useVocalTract+generationScaled+useVocalTract*generationScaled,data=subset(simdata,generation<=100))
summary(proFitModel)

perFitModel = lmer(medPerFitnessScaled~(1|runfolder)+useVocalTract+generationScaled+useVocalTract*generationScaled,data=subset(simdata,generation<=100))
summary(perFitModel)

var.test(subset(subset(simdata,useVocalTract==0),generation==100)$medProFitness,subset(subset(simdata,useVocalTract==1),generation==100)$medProFitness)
var.test(subset(subset(simdata,useVocalTract==0),generation==100)$medPerFitness,subset(subset(simdata,useVocalTract==1),generation==100)$medPerFitness)

var.test(subset(subset(simdata,useVocalTract==0),generation==500)$medProFitness,subset(subset(simdata,useVocalTract==1),generation==500)$medProFitness)
var.test(subset(subset(simdata,useVocalTract==0),generation==500)$medPerFitness,subset(subset(simdata,useVocalTract==1),generation==500)$medPerFitness)