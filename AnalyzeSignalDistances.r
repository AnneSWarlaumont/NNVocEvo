# Anne S. Warlaumont

# Determine how far apart the signal types are in 6-dimensional space as well as how far apart the exemplars are within a signal type
# Also run k means clustering on the signal exemplars

# To run:
# R CMD BATCH AnalyzeSignalDistances.r

rm(list=ls())

library(mclust)

rundirs = c(
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20131126T160313/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20131126T160318/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20131126T160323/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20131126T160328/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20131126T160333/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20131126T160338/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20131126T160348/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20140327T154338/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20140327T154343/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20140327T154348/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20140527T171008/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20140527T171011/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20140527T171016/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20140527T171022/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20140527T171027/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20140620T141715/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20140620T141720/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20140620T141725/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20140620T141730/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20140620T141736/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160354/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160359/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160404/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160409/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160414/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160425/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160430/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20140327T154353/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20140327T154358/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20140327T154403/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20140612T105820/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20140612T105825/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20140612T105831/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20140612T105836/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20140612T105841/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20140620T141341/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20140620T141346/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20140620T141351/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20140620T141356/",
"/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20140620T141401/")

versions = c("Abstract","Abstract","Abstract","Abstract","Abstract","Abstract","Abstract","Abstract","Abstract","Abstract","Abstract","Abstract","Abstract","Abstract","Abstract","Abstract","Abstract","Abstract","Abstract","Abstract","Realistic","Realistic","Realistic","Realistic","Realistic","Realistic","Realistic","Realistic","Realistic","Realistic","Realistic","Realistic","Realistic","Realistic","Realistic","Realistic","Realistic","Realistic","Realistic","Realistic")

allsims = data.frame(rundir=rundirs,version=versions)

for (count in 1:length(rundirs)) {
	
	rundir = rundirs[count]
	
	sig1_full = read.csv(paste(rundir,"perceiverParentInputsDiary_gen500_sig1.txt",sep=""),"header"=FALSE,"sep"="\t")
	sig1 = sig1_full[seq(1,10000,100),1:6]
	
	sig2_full = read.csv(paste(rundir,"perceiverParentInputsDiary_gen500_sig2.txt",sep=""),"header"=FALSE,"sep"="\t")
	sig2 = sig2_full[seq(1,10000,100),1:6]
	
	sig3_full = read.csv(paste(rundir,"perceiverParentInputsDiary_gen500_sig3.txt",sep=""),"header"=FALSE,"sep"="\t")
	sig3 = sig3_full[seq(1,10000,100),1:6]
	
	sig1_dist = dist(sig1)
	sig1_dist_mean = mean(sig1_dist)
	
	sig2_dist = dist(sig2)
	sig2_dist_mean = mean(sig2_dist)
	
	sig3_dist = dist(sig3)
	sig3_dist_mean = mean(sig3_dist)
	
	sig12_dist = dist(rbind(sig1,sig2))
	sig12_distbetween_mean = (sum(sig12_dist)-sum(sig1_dist)-sum(sig2_dist))/10000
	
	sig13_dist = dist(rbind(sig1,sig3))
	sig13_distbetween_mean = (sum(sig13_dist)-sum(sig1_dist)-sum(sig3_dist))/10000
	
	sig23_dist = dist(rbind(sig2,sig3))
	sig23_distbetween_mean = (sum(sig23_dist)-sum(sig2_dist)-sum(sig3_dist))/10000
	
	within_sig_dist_mean = mean(c(sig1_dist_mean,sig2_dist_mean,sig3_dist_mean))
	between_sig_dist_mean = mean(c(sig12_distbetween_mean,sig13_distbetween_mean,sig23_distbetween_mean))
	within_between_sig_dist_mean_ratio = within_sig_dist_mean/between_sig_dist_mean
	
	print(rundir)
	print(within_sig_dist_mean)
	print(between_sig_dist_mean)
	print(within_sig_dist_mean/between_sig_dist_mean)
	
	allsims$within_sig_dist_mean[count] = within_sig_dist_mean
	allsims$between_sig_dist_mean[count] = between_sig_dist_mean
	allsims$within_between_sig_dist_mean_ratio[count] = within_between_sig_dist_mean_ratio
	
	# Plot the variance explained when different numbers of clustered are used with k means clustering
	# Code adapted from http://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
	allsigs = rbind(sig1,sig2,sig3)
	wss = (nrow(allsigs)-1)*sum(apply(allsigs,2,var))
	for (k in 2:15){
		wss[k] = sum(kmeans(allsigs,k)$withinss)
	}
	plot(1:15,wss,type="b",xlab="Number of Clusters",ylab="Within groups sum of squares")
	
	# Get the optimal number of clusters using the Bayesian Information Criterion (BIC)
	# Code adapted from http://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
	# and //www.jstatsoft.org/v18/i06/paper
	# http://cran.r-project.org/web/packages/mclust/mclust.pdf
	# Run the function to see how many clusters
	# it finds to be optimal, set it to search for
	# at least 1 model and up 20.
	allsigs_clust = Mclust(as.matrix(allsigs), G=1:20)
	m.best = dim(allsigs_clust$z)[2]
	cat("Model-based optimal number of clusters:",m.best,"\n")
	allsims$numclust[count] = m.best
	#plot(allsigs_clust)
	
}

mean(subset(allsims,version=="Abstract")$numclust)
mean(subset(allsims,version=="Realistic")$numclust)
t.test(subset(allsims,version=="Abstract")$numclust,subset(allsims,version=="Realistic")$numclust)

mean(subset(allsims,version=="Abstract")$within_sig_dist_mean)
mean(subset(allsims,version=="Abstract")$between_sig_dist_mean)
mean(subset(allsims,version=="Abstract")$within_between_sig_dist_mean_ratio)
sd(subset(allsims,version=="Abstract")$within_between_sig_dist_mean_ratio)

mean(subset(allsims,version=="Realistic")$within_sig_dist_mean)
mean(subset(allsims,version=="Realistic")$between_sig_dist_mean)
mean(subset(allsims,version=="Realistic")$within_between_sig_dist_mean_ratio)
sd(subset(allsims,version=="Realistic")$within_between_sig_dist_mean_ratio)

t.test(subset(allsims,version=="Abstract")$within_sig_dist_mean,subset(allsims,version=="Realistic")$within_sig_dist_mean)
t.test(subset(allsims,version=="Abstract")$between_sig_dist_mean,subset(allsims,version=="Realistic")$between_sig_dist_mean)
t.test(subset(allsims,version=="Abstract")$within_between_sig_dist_mean_ratio,subset(allsims,version=="Realistic")$within_between_sig_dist_mean_ratio)

