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

aggdataPro = aggregate(medProFitness~runfolder+useVocalTract,subset(simdata,generation>25 & generation < 150),mean)
t.test(subset(aggdataPro,useVocalTract==0)$medProFitness,subset(aggdataPro,useVocalTract==1)$medProFitness)

aggdataPer = aggregate(medPerFitness~runfolder+useVocalTract,subset(simdata,generation>25 & generation < 150),mean)
t.test(subset(aggdataPer,useVocalTract==0)$medPerFitness,subset(aggdataPer,useVocalTract==1)$medPerFitness)

proFitModel = lmer(medProFitnessScaled~(1|runfolder)+useVocalTract+generationScaled+useVocalTract*generationScaled,data=subset(simdata,generation<=100 & generation > 50))
summary(proFitModel)
quartz(); qqnorm(resid(proFitModel))

perFitModel = lmer(medPerFitnessScaled~(1|runfolder)+useVocalTract+generationScaled+useVocalTract*generationScaled,data=subset(simdata,generation<=100 & generation > 50))
summary(perFitModel)
quartz(); qqnorm(resid(perFitModel))

var.test(subset(subset(simdata,useVocalTract==0),generation==100)$medProFitness,subset(subset(simdata,useVocalTract==1),generation==100)$medProFitness)
var.test(subset(subset(simdata,useVocalTract==0),generation==100)$medPerFitness,subset(subset(simdata,useVocalTract==1),generation==100)$medPerFitness)

var.test(subset(subset(simdata,useVocalTract==0),generation==500)$medProFitness,subset(subset(simdata,useVocalTract==1),generation==500)$medProFitness)
var.test(subset(subset(simdata,useVocalTract==0),generation==500)$medPerFitness,subset(subset(simdata,useVocalTract==1),generation==500)$medPerFitness)
