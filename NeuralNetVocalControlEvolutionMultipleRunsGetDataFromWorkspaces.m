function [] = NeuralNetVocalControlEvolutionMultipleRunsGetDataFromWorkspaces(runsFolder,csvFile,longVer)

% Anne S. Warlaumont
%
% Example usage:
% NeuralNetVocalControlEvolutionMultipleRunsGetDataFromWorkspaces('/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/','/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/MultipleRunsData.csv',1)
% Set longVer to plot detailed muscle and signal vectors, save the sounds of those signals, and plot genes as they change over time. This will take a lot longer to run.

% Specify where the simulation data are located:
workspaceFolders = {};
dirs = dir([runsFolder,'AbstractReady/run*']);
for run = 1:size(dirs,1)
    workspaceFolders = [workspaceFolders;['AbstractReady/',dirs(run).name]];
end
dirs = dir([runsFolder,'RealisticReady/run*']);
for run = 1:size(dirs,1)
    workspaceFolders = [workspaceFolders;['RealisticReady/',dirs(run).name]];
end

% Specify the number of generations to analyze:
numgens = 500;

% Open the file simulation summary information file for writing, for later analysis in R:
fid2 = fopen(csvFile,'w');
% Write column headers into that summary information file:
fprintf(fid2,'runfolder,useVocalTract,mutationProbability,numIndividuals,generation,medProFitness,medPerFitness\n');

% Initialize the figure where median fitness values across generation will be plotted:
medianTrajFig = figure('visible','off');
set(0,'defaultaxesfontsize',12);

% Initialize counters for the number of Abstract (numUseVocalTract0) and Realistic (numUseVocalTract1) simulations
numUseVocalTract0 = 0;
numUseVocalTract1 = 0;

% Initialize matrices that will store the median producer and perceiver fitnesses within each simulation at each generation, as well as a vector that will store whether each simulation used Abstract or Realistic signals
medianProdFitnesses=NaN(size(workspaceFolders,1),numgens);
medianPercFitnesses=NaN(size(workspaceFolders,1),numgens);
useVocalTracts = NaN(size(workspaceFolders,1),1);

% Iterate through all the simulations
for run=1:size(workspaceFolders,1)
	
	display(['run ',num2str(run),' of ',num2str(size(workspaceFolders,1))]);
    
    % Load the simulation's MATLAB workspace
    load([runsFolder,char(workspaceFolders(run,1)),'/NeuralNetVocalControlEvolutionWorkspace.mat']);
    display(['Loaded ',runsFolder,char(workspaceFolders(run,1)),'/NeuralNetVocalControlEvolutionWorkspace.mat']);
    
    medianProdFitness = median(producerParentFitnessesDiary,1);
    medianProdFitnesses(run,:) = medianProdFitness(1,1:numgens);
    medianPercFitness = median(perceiverParentFitnessesDiary,1);
    medianPercFitnesses(run,:) = medianPercFitness(1,1:numgens);
    useVocalTracts(run,1) = useVocalTract;
    
    if useVocalTract == 0
        numUseVocalTract0 = numUseVocalTract0 + 1;
        subplot(2,2,1); hLine1 = plot(median(producerParentFitnessesDiary(:,1:numgens)),'Color',[.5,.5,.5],'LineWidth',1); xlabel('Generation'); ylabel('Median producer fitness score'); ylim([0,4]); hold on;
        subplot(2,2,3); hLine2 = plot(median(perceiverParentFitnessesDiary(:,1:numgens)),'Color',[.5,.5,.5],'LineWidth',1); xlabel('Generation'); ylabel('Median perceiver fitness score'); ylim([0,4]); hold on;
    elseif useVocalTract == 1
        numUseVocalTract1 = numUseVocalTract1 + 1;
        figure(medianTrajFig);
        subplot(2,2,2); hLine1 = plot(median(producerParentFitnessesDiary(:,1:numgens)),'Color','black','LineWidth',1); xlabel('Generation'); ylabel('Median producer fitness score'); ylim([0,4]); hold on;
        subplot(2,2,4); hLine2 = plot(median(perceiverParentFitnessesDiary(:,1:numgens)),'Color','black','LineWidth',1); xlabel('Generation'); ylabel('Median perceiver fitness score'); ylim([0,4]); hold on;
    end
    
    % Save the individuals' spectrograms and waveforms, Praat scripts, and muscle outputs (where applicable).
    
    for signalNum=1:3
        temp = perceiverParentInputsDiary{numgens,signalNum};
        save([runsFolder,char(workspaceFolders(run,1)),'/perceiverParentInputsDiary_gen',num2str(numgens),'_sig',num2str(signalNum),'.txt'],'temp','-ascii','-tabs');
    end
        
    for gencount=1:numgens
        fprintf(fid2,[[runsFolder,char(workspaceFolders(run,1))],',',num2str(useVocalTract),',',num2str(mutationProbability),',',num2str(numIndividuals),',',num2str(gencount),',',num2str(median(producerParentFitnessesDiary(:,gencount))),',',num2str(median(perceiverParentFitnessesDiary(:,gencount))),'\n']);
    end
    
    if longVer
        
        [~,sortIndGen1]=sort(producerParentFitnessesDiary(:,1),'descend');
        [~,sortInd]=sort(producerParentFitnessesDiary(:,numgens),'descend');
        
        % Generate all signals and save their images and, if applicable, their sounds, at generation 1 and generation numgens:
        for savesiggen = [1,numgens]
            for producerParent=1:numIndividuals
                if savesiggen==1
                    fitRank = find(sortIndGen1==producerParent);
                elseif savesiggen==numgens
                    fitRank = find(sortInd==producerParent);
                end
                for signalNum=1:3
                    muscleOutputsReshaped = reshape((producerParentOutputsDiary{savesiggen,1}(2*((producerParent-1)*3+signalNum),:)+1)/2,3,2);
                    caxisBounds = [0,1];
                    muscleOutputsFilename = [runsFolder,char(workspaceFolders(run,1)),'/generation',num2str(savesiggen),'_producerParent',num2str(producerParent),'_fitRank',num2str(fitRank),'_vocalization',num2str(signalNum),'_muscleOutputs'];
                    muscleOutputs_fig = figure('visible','off'); colormap(flipud(gray)); image(flipud(muscleOutputsReshaped),'CDataMapping','scaled'); caxis(caxisBounds); set (gca,'XTick',[]); set(gca,'YTick',[]); print(muscleOutputs_fig,'-dtiff',muscleOutputsFilename); close(muscleOutputs_fig);
                end
                saveSounds = 1;
                [producerParentOutputsDiary,~,~,~,~] = getProducerParentSounds(producerParent,producerParentGenesDiary{savesiggen,1},numProducerNetInputs,numProducerNetHidden,numProducerNetOutputs,numSignals,producerInputs,[runsFolder,char(workspaceFolders(run,1)),'/'],savesiggen,useVocalTract,timestep,duration,producerParentOutputsDiary,melfcc_wintime,melfcc_nbands,melfcc_hoptime,numIndividuals,perceiverParentGenesDiary{savesiggen,1},numPerceiverNetInputs,numPerceiverNetHidden,numPerceiverNetOutputs,perceiverTargets,perceiverParentOutputsDiary,perceiverParentInputsDiary,perceiverParentCorrectness,saveSounds);
            end
        end
        
        % Plot the genes of all individuals having a given rank (constantRank) as the population evolves:
        % TODO: replace the code below with a call to genesOverTime.m
        constantRankProGenes=[];
        constantRankPerGenes=[];
        constantRank = 50;
        constantFitProGenesFilename = [runsFolder,char(workspaceFolders(run,1)),'/fitRank',num2str(constantRank),'_producergenes_gen1to',num2str(numgens),'.tif'];
        constantFitPerGenesFilename = [runsFolder,char(workspaceFolders(run,1)),'/fitRank',num2str(constantRank),'_perceivergenes_gen1to',num2str(numgens),'.tif'];
        % Producers:
        for generationcount=1:numgens
            [~,sortInd]=sort(producerParentFitnessesDiary(:,generationcount),'descend');
            constantRankInd(generationcount) = sortInd(constantRank,1);
            constantRankProGenes=[constantRankProGenes,producerParentGenesDiary{generationcount,1}(constantRankInd(generationcount),:)'];
        end
        constantFitProGenes_fig = figure('visible','off');
        imagesc(constantRankProGenes);xlabel('Generation'),ylabel('Gene number'); colormap(flipud(gray)); colorbar();
        print(constantFitProGenes_fig,'-dtiff',constantFitProGenesFilename); close(constantFitProGenes_fig);
        % Perceivers:
        for generationcount=1:numgens
            [~,sortInd]=sort(perceiverParentFitnessesDiary(:,generationcount),'descend');
            constantRankInd(generationcount) = sortInd(constantRank,1);
            constantRankPerGenes=[constantRankPerGenes,perceiverParentGenesDiary{generationcount,1}(constantRankInd(generationcount),:)'];
        end
        constantFitPerGenes_fig = figure('visible','off');
        imagesc(constantRankPerGenes);xlabel('Generation'),ylabel('Gene number'); colormap(flipud(gray)); colorbar();
        print(constantFitPerGenes_fig,'-dtiff',constantFitPerGenesFilename); close(constantFitPerGenes_fig);
        
    end
    
    clearvars -except runsFolder fid2 workspaceFolders mutationProbabilities mutationSDs medianFitnessAfters numgens medianTrajFig numUseVocalTract0 numUseVocalTract1 medianProdFitnesses medianPercFitnesses useVocalTracts numgens longVer
    
end

fclose(fid2);

% Add titles and legend, where appropriate, to the plots of individual runs, then save the plots
figure(medianTrajFig);
subplot(2,2,1);
title('Abstract');
subplot(2,2,2);
title('Realistic');
print(medianTrajFig,'-depsc',[runsFolder,'fig_Medians.eps']); close(medianTrajFig);

save([runsFolder,'NeuralNetVocalControlEvolutionMultipleRunsGetDataFromWorkspaces.mat']);

% Make a figure that shows mean of the median fitness across different runs
% TODO: Add confidence intervals (might be easier in R)
fig_MeanOfMedian = figure;
subplot(2,1,1);
plot(mean(medianProdFitnesses(find(useVocalTracts==0),:)),'Color',[.5,.5,.5]); xlabel('Generation'); ylabel('Producers'); title('Mean median fitness score'); hold on;
plot(mean(medianProdFitnesses(find(useVocalTracts==1),:)),'Color','black'); hold on;
legend('Abstract','Embodied','Location','SouthEast'); legend('boxoff');
subplot(2,1,2);
plot(mean(medianPercFitnesses(find(useVocalTracts==0),:)),'Color',[.5,.5,.5]); xlabel('Generation'); ylabel('Perceivers'); hold on;
plot(mean(medianPercFitnesses(find(useVocalTracts==1),:)),'Color','black'); hold on;
print('-depsc',[runsFolder,'fig_MeanOfMedian.eps']); close(fig_MeanOfMedian);

% Make a figure that shows SD of the median fitness across different runs
fig_SDOfMedian = figure;
subplot(2,1,1);
plot(std(medianProdFitnesses(find(useVocalTracts==0),:)),'Color',[.5,.5,.5]); xlabel('Generation'); ylabel('Producers'); title('Standard deviation of median fitness score'); hold on;
plot(std(medianProdFitnesses(find(useVocalTracts==1),:)),'Color','black'); hold on;
legend('Abstract','Embodied','Location','NorthEast'); legend('boxoff');
subplot(2,1,2);
plot(std(medianPercFitnesses(find(useVocalTracts==0),:)),'Color',[.5,.5,.5]); xlabel('Generation'); ylabel('Receivers'); hold on;
plot(std(medianPercFitnesses(find(useVocalTracts==1),:)),'Color','black'); hold on;
print('-depsc',[runsFolder,'fig_SDOfMedian.eps']); close(fig_SDOfMedian);

% Make a figure that shows coefficient of variation of the median fitness across different runs
fig_CoVOfMedian = figure;
subplot(2,1,1);
plot(std(medianProdFitnesses(find(useVocalTracts==0),:))./mean(medianProdFitnesses(find(useVocalTracts==0),:)),'Color',[.5,.5,.5]); xlabel('Generation'); ylabel('Producers'); title('Coefficient of variation of median fitness score'); hold on;
plot(std(medianProdFitnesses(find(useVocalTracts==1),:))./mean(medianProdFitnesses(find(useVocalTracts==1),:)),'Color','black'); hold on;
legend('Abstract','Embodied','Location','NorthEast'); legend('boxoff');
subplot(2,1,2);
plot(std(medianPercFitnesses(find(useVocalTracts==0),:))./mean(medianPercFitnesses(find(useVocalTracts==0),:)),'Color',[.5,.5,.5]); xlabel('Generation'); ylabel('Perceivers'); hold on;
plot(std(medianPercFitnesses(find(useVocalTracts==1),:))./mean(medianPercFitnesses(find(useVocalTracts==1),:)),'Color','black'); hold on;
print('-depsc',[runsFolder,'fig_CoVOfMedian.eps']); close(fig_CoVOfMedian);

