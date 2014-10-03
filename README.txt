This project focuses on the simulation of the evolution of neural networks that produce and receive vocal signals. The aim is to model the evolution of reflexive signals, such as screams and cries, in primates. It compares models of vocalization that treat signals as abstract entities to models that consider vocal tract physiology and acoustics.

Most of the project is in MATLAB. Some statistical analyses are done in R.

Main simulation code:
NeuralNetVocalControlEvolution.m, which depends on getProducerParentSounds.m

Plotting and analysis code:
NeuralNetVocalControlEvolutionMultipleRunsGetDataFromWorkspaces.m
genesOverTime.m
landscape.m
wavefigs.m
AnalyzeFitness.r
AnalyzeSignalDistances.r

Author: Anne S. Warlaumont
