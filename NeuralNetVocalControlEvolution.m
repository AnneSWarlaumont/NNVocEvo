function [] = NeuralNetVocalControlEvolution(PraatDir,numIndividuals,numGenerations,mutationProbability,useVocalTract,producerInputs,perceiverTargets,varargin)

% Anne S. Warlaumont
% Collaboration with Andrew Olney
% Evolving a simple recurrent network for human vocal motor control
% "Muscles" used: Lungs, Interarytenoid, Cricothyroid, Vocalis, Thyroarytenoid, PosteriorCricoarytenoid, LateralCricoarytenoid
% Upper vocal tract muscles are ignored.
% PraatDir sets the directory in which Praat scripts and wav files will be written
%
% Example usage:
% NeuralNetVocalControlEvolution('/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/Realistic/run20140527T170143/',100,500,.05,1,[1,0,0;0,1,0;0,0,1],[1,0,0;0,1,0;0,0,1])
% NeuralNetVocalControlEvolution('/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/ExtremeWeightsExp/run20140812_capWeights_1/',100,500,.05,0,[1,0,0;0,1,0;0,0,1],[1,0,0;0,1,0;0,0,1],'capWeights',1)
% NeuralNetVocalControlEvolution('/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/ExtremeWeightsExp/run20140812_arbitraryFitness_1/',100,500,.05,0,[1,0,0;0,1,0;0,0,1],[1,0,0;0,1,0;0,0,1],'arbitraryFitness',1)

% NeuralNetVocalControlEvolution('/Users/awarlau/Downloads/noiseAmnt_1_run_1/',100,500,.05,1,[1,0,0;0,1,0;0,0,1],[1,0,0;0,1,0;0,0,1],'noiseAmnt',1)

% NeuralNetVocalControlEvolution('/Users/awarlau/Downloads/noiseAmnt_1_maxFreq_500_noNN_1_run_1/',50,500,.05,1,[1,0,0;0,1,0;0,0,1],[1,0,0;0,1,0;0,0,1],'noiseAmnt',1,'maxFreq',500,'noNN',1)

% NeuralNetVocalControlEvolution('/Users/awarlau/Downloads/noiseAmnt_Pt001_nInd_50/',50,500,.05,1,[1,0,0;0,1,0;0,0,1],[1,0,0;0,1,0;0,0,1],'noiseAmnt',.001)

p = inputParser;
addRequired(p,'PraatDir');
addRequired(p,'numIndividuals');
addRequired(p,'numGenerations');
addRequired(p,'mutationProbability');
addRequired(p,'useVocalTract');
addRequired(p,'producerInputs');
addRequired(p,'perceiverTargets');
addParameter(p,'capWeights',0);
addParameter(p,'arbitraryFitness',0);
addParameter(p,'crossover',0);
addParameter(p,'noiseAmnt',0);
addParameter(p,'maxFreq',2000);
addParameter(p,'noNN',0)
parse(p,PraatDir,numIndividuals,numGenerations,mutationProbability,useVocalTract,producerInputs,perceiverTargets,varargin{:});

capWeights = p.Results.capWeights;
arbitraryFitness = p.Results.arbitraryFitness;
crossover = p.Results.crossover;
noiseAmnt = p.Results.noiseAmnt;
maxFreq = p.Results.maxFreq;
noNN = p.Results.noNN;

if ~exist(PraatDir, 'dir')
    mkdir(PraatDir);
end

addpath(genpath('rastamat'));
addpath(genpath('randp'));

if exist([PraatDir,'NeuralNetVocalControlEvolutionWorkspace.mat']) == 0
    
    % Set simulation time and time step size paramenters:
    duration = .5;
    timestep = .5;
    melfcc_wintime = .3;
    melfcc_hoptime = .2;
    melfcc_nbands = 3;
    
    % Set the number of fixed signals that the population should produce:
    numSignals = size(producerInputs,1);
    
    numProducerNetOutputs = 6; % This is the number of muscles
    numProducerNetInputs = size(producerInputs,2);
	if ~noNN
		numProducerNetHidden = 3; % This is chosen rather arbitrarily
	else
		numProducerNetHidden = 0;
	end
    
    numPerceiverNetInputs = melfcc_nbands*2; % This will have to be changed to the correct value, which is the number of pixels in each spectrogram.
    numPerceiverNetOutputs = size(producerInputs,2);
	if ~noNN
 		numPerceiverNetHidden = 3; % This is chosen rather arbitrarily
	else
		numPerceiverNetHidden = 0;
	   
    end
    
    % numGenes = # of input to hidden weights, # of context to hidden weights,
    % # of hidden to output weights, # of bias to hidden weights, # of bias to
    % output weights
    if ~noNN
		numProducerGenes = (numProducerNetInputs*numProducerNetHidden)+(numProducerNetHidden*numProducerNetOutputs);
	    numPerceiverGenes = (numPerceiverNetInputs*numPerceiverNetHidden)+(numPerceiverNetHidden*numPerceiverNetOutputs);
	else
		numProducerGenes = numProducerNetInputs*numProducerNetOutputs;
		numPerceiverGenes = numPerceiverNetInputs*numPerceiverNetOutputs;
	end
    
    % Randomly initialize parent genes
    producerParentGenes = rand(numIndividuals,numProducerGenes)*2-1;
    perceiverParentGenes = rand(numIndividuals,numPerceiverGenes)*2-1;
    
    producerParentFitnessesDiary = zeros(numIndividuals,numGenerations);
    perceiverParentFitnessesDiary = zeros(numIndividuals,numGenerations);
    parentSig1LoudnessDiary = zeros(numIndividuals,numGenerations);
    parentSig2LoudnessDiary = zeros(numIndividuals,numGenerations);
    parentSig3LoudnessDiary = zeros(numIndividuals,numGenerations);
    parentSigDiffDiary = zeros(numIndividuals,numGenerations);
    producerParentGenesDiary = cell(numGenerations,1);
    perceiverParentGenesDiary = cell(numGenerations,1);
    producerParentOutputsDiary = cell(numGenerations,numSignals);
    perceiverParentOutputsDiary = cell(numGenerations,numSignals);
    perceiverParentInputsDiary = cell(numGenerations,numSignals);
    generation = 1;
else
    load([PraatDir,'NeuralNetVocalControlEvolutionWorkspace.mat'],'-regexp','^(?!numGenerations$).')
	
    producerParentFitnessesDiary = [producerParentFitnessesDiary,zeros(numIndividuals,numGenerations-generation)];
    perceiverParentFitnessesDiary = [perceiverParentFitnessesDiary,zeros(numIndividuals,numGenerations-generation)];
    parentSig1LoudnessDiary = [parentSig1LoudnessDiary,zeros(numIndividuals,numGenerations-generation)];
    parentSig2LoudnessDiary = [parentSig2LoudnessDiary,zeros(numIndividuals,numGenerations-generation)];
	parentSig3LoudnessDiary = [parentSig3LoudnessDiary,zeros(numIndividuals,numGenerations-generation)];
    parentSigDiffDiary = [parentSigDiffDiary,zeros(numIndividuals,numGenerations-generation)];
    producerParentGenesDiary{numGenerations,1} = [];
    perceiverParentGenesDiary{numGenerations,1} = [];
    producerParentOutputsDiary{numGenerations,numSignals} = [];
    perceiverParentOutputsDiary{numGenerations,numSignals} = [];
    perceiverParentInputsDiary{numGenerations,numSignals} = [];
	
	generation = generation + 1;
end

for generation = generation:numGenerations
    
    producerParentGenesDiary{generation} = producerParentGenes;
    perceiverParentGenesDiary{generation} = perceiverParentGenes;
    
    perceiverParentCorrectness = zeros(numIndividuals,numIndividuals);
    producerParentFitnesses = zeros(numIndividuals,1);
    perceiverParentFitnesses = zeros(numIndividuals,1);
    
    sig1Loudnesses = zeros(numIndividuals,1);
    sig2Loudnesses = zeros(numIndividuals,1);
	sig3Loudnesses = zeros(numIndividuals,1);
    
    for producerParent=1:numIndividuals
        
        saveSounds = 1;
        [producerParentOutputsDiary,aspectra,perceiverParentOutputsDiary,perceiverParentInputsDiary,perceiverParentCorrectness] = getProducerParentSounds(producerParent,producerParentGenes,numProducerNetInputs,numProducerNetHidden,numProducerNetOutputs,numSignals,producerInputs,PraatDir,generation,useVocalTract,timestep,duration,producerParentOutputsDiary,melfcc_wintime,melfcc_nbands,melfcc_hoptime,numIndividuals,perceiverParentGenes,numPerceiverNetInputs,numPerceiverNetHidden,numPerceiverNetOutputs,perceiverTargets,perceiverParentOutputsDiary,perceiverParentInputsDiary,perceiverParentCorrectness,saveSounds,noiseAmnt,maxFreq,noNN);
        
        sig1Loudnesses(producerParent,1) = sum(sum(aspectra{1}));
        sig2Loudnesses(producerParent,1) = sum(sum(aspectra{2}));
		sig3Loudnesses(producerParent,1) = sum(sum(aspectra{3}));
        
    end
    
    for producerParent = 1:numIndividuals
        producerParentFitnesses(producerParent,1) = mean(perceiverParentCorrectness(:,producerParent))+.01;
    end
    for perceiverParent = 1:numIndividuals
        perceiverParentFitnesses(perceiverParent,1) = mean(perceiverParentCorrectness(perceiverParent,:))'+.01;
    end
    
    parentSig1LoudnessDiary(:,generation) = sig1Loudnesses;
    parentSig2LoudnessDiary(:,generation) = sig2Loudnesses;
	parentSig3LoudnessDiary(:,generation) = sig3Loudnesses;
    producerParentFitnessesDiary(:,generation) = producerParentFitnesses;
    perceiverParentFitnessesDiary(:,generation) = perceiverParentFitnesses;
    
    [~,ind] = sort(producerParentFitnesses,'ascend');
    for producerParent = 1:numIndividuals
        producerParentFitnessesOrdinal(producerParent) = find(ind==producerParent);
    end
    
    [~,ind] = sort(perceiverParentFitnesses,'ascend');
    for perceiverParent = 1:numIndividuals
        perceiverParentFitnessesOrdinal(perceiverParent) = find(ind==perceiverParent);
    end
    
    % Select pairs of producer parents to reproduce:
    producerParentPairs = zeros(numIndividuals,2);
    for pairNum = 1:numIndividuals;
        pairFound=0;
        while(~pairFound)
            if arbitraryFitness
                producerParentPairs(pairNum,:)=randp(ones(size(producerParentFitnesses)),[1,2]);
            else
                producerParentPairs(pairNum,:)=randp(producerParentFitnesses,[1,2]);
            end
            if producerParentPairs(pairNum,1)~=producerParentPairs(pairNum,2)
                pairFound=1;
            end
        end
    end
    
    % Create the producer children:
    producerChildGenes=zeros(numIndividuals,numProducerGenes);
    for child=1:numIndividuals;
        if crossover
            for gene=1:numProducerGenes
                producerChildGenes(child,gene)=producerParentGenes(producerParentPairs(child,round(rand(1))+1),gene);
            end
            producerChildGenes(child,:) = mean(producerParentGenes(producerParentPairs(child,:),:));
        else
            producerChildGenes(child,:) = producerParentGenes(producerParentPairs(child,1),:);
        end
    end
    
    % Mutate the producer children:
    mutationLoci=randp([1-mutationProbability,mutationProbability],size(producerChildGenes))-1;
    mutationAmounts=(mutationLoci.*((rand(size(producerChildGenes))-.5)/.7));
    producerChildGenes(1:end)=producerChildGenes(1:end)+mutationAmounts(1:end);
    if capWeights
        producerChildGenes(producerChildGenes<-1)=-1;
        producerChildGenes(producerChildGenes>1)=1;
    end
    
    % Make the producer children the producer parents:
    producerParentGenes=producerChildGenes;
    
    % Select pairs of perceiver parents to reproduce:
    perceiverParentPairs = zeros(numIndividuals,2);
    for pairNum = 1:numIndividuals;
        pairFound=0;
        while(~pairFound) 
            if arbitraryFitness
                perceiverParentPairs(pairNum,:)=randp(ones(size(perceiverParentFitnesses)),[1,2]);
            else
                perceiverParentPairs(pairNum,:)=randp(perceiverParentFitnesses,[1,2]);
            end
            if perceiverParentPairs(pairNum,1)~=perceiverParentPairs(pairNum,2)
                pairFound=1;
            end
        end
    end
    
    % Create the perceiver children:
    perceiverChildGenes=zeros(numIndividuals,numPerceiverGenes);
    for child=1:numIndividuals;
        if crossover
            for gene=numPerceiverGenes
                perceiverChildGenes(child,gene)=perceiverParentGenes(perceiverParentPairs(child,round(rand(1))+1),gene);
            end
            perceiverChildGenes(child,:) = mean(perceiverParentGenes(perceiverParentPairs(child,:),:));
        else
            perceiverChildGenes(child,:) = perceiverParentGenes(perceiverParentPairs(child,1),:);
        end
    end
    
    % Mutate the perceiver children:
    mutationLoci=randp([1-mutationProbability,mutationProbability],size(perceiverChildGenes))-1;
    mutationAmounts=(mutationLoci.*((rand(size(perceiverChildGenes))-.5)/.7));
    perceiverChildGenes(1:end)=perceiverChildGenes(1:end)+mutationAmounts(1:end);
    if capWeights
        perceiverChildGenes(perceiverChildGenes<-1)=-1;
        perceiverChildGenes(perceiverChildGenes>1)=1;
    end
    
    % Make the perceiver children the perceiver parents:
    perceiverParentGenes=perceiverChildGenes;
    
    maxGeneration=generation;
    
    myfig = figure;
    subplot(1,5,1); plot(median(producerParentFitnessesDiary(:,1:maxGeneration))); xlabel('Generation');ylabel('Median producer fitness score');
    subplot(1,5,2); plot(median(perceiverParentFitnessesDiary(:,1:maxGeneration))); xlabel('Generation');ylabel('Median perceiver fitness score');
    subplot(1,5,3); plot(median(parentSig1LoudnessDiary(:,1:maxGeneration))); xlabel('Generation'); ylabel('Median signal 1 loudness');
    subplot(1,5,4); plot(median(parentSig2LoudnessDiary(:,1:maxGeneration))); xlabel('Generation'); ylabel('Median signal 2 loudness');
    subplot(1,5,5); plot(median(parentSig3LoudnessDiary(:,1:maxGeneration))); xlabel('Generation'); ylabel('Median signal 3 loudness');
    %set(myfig,'Units','inches');
    %set(myfig,'Position',[1,1,18,3]);
    set(myfig,'PaperPositionMode','manual');
    set(myfig,'PaperUnits','inches');
    set(myfig,'PaperSize',[16 4]);
    set(myfig,'PaperPosition',[.5,.5,15,3]);
    print(myfig,'-dpdf',[PraatDir,'MedianPerformanceAcrossGenerations.pdf']);
    close(myfig);
    
    myfig = figure;
    subplot(1,5,1); plot(max(producerParentFitnessesDiary(:,1:maxGeneration))); xlabel('Generation');ylabel('Max producer fitness score');
    subplot(1,5,2); plot(max(perceiverParentFitnessesDiary(:,1:maxGeneration))); xlabel('Generation');ylabel('Max perceiver fitness score');
    subplot(1,5,3); plot(max(parentSig1LoudnessDiary(:,1:maxGeneration))); xlabel('Generation'); ylabel('Max signal 1 loudness');
    subplot(1,5,4); plot(max(parentSig2LoudnessDiary(:,1:maxGeneration))); xlabel('Generation'); ylabel('Max signal 2 loudness');
    subplot(1,5,5); plot(max(parentSig3LoudnessDiary(:,1:maxGeneration))); xlabel('Generation'); ylabel('Max signal 3 loudness');
    %set(myfig,'Units','inches');
    %set(myfig,'Position',[1,1,18,3]);
    set(myfig,'PaperPositionMode','manual');
    set(myfig,'PaperUnits','inches');
    set(myfig,'PaperSize',[16 4]);
    set(myfig,'PaperPosition',[.5,.5,15,3]);
    print(myfig,'-dpdf',[PraatDir,'MaxPerformanceAcrossGenerations.pdf']);
    close(myfig);
    
    myfig = figure;
    subplot(1,5,1); plot(min(producerParentFitnessesDiary(:,1:maxGeneration))); xlabel('Generation');ylabel('Min producer fitness score');
    subplot(1,5,2); plot(min(producerParentFitnessesDiary(:,1:maxGeneration))); xlabel('Generation');ylabel('Min perceiver fitness score');
    subplot(1,5,3); plot(min(parentSig1LoudnessDiary(:,1:maxGeneration))); xlabel('Generation'); ylabel('Min signal 1 loudness');
    subplot(1,5,4); plot(min(parentSig2LoudnessDiary(:,1:maxGeneration))); xlabel('Generation'); ylabel('Min signal 2 loudness');
    subplot(1,5,5); plot(min(parentSig3LoudnessDiary(:,1:maxGeneration))); xlabel('Generation'); ylabel('Min signal 3 loudness');
    %set(myfig,'Units','inches');
    %set(myfig,'Position',[1,1,18,3]);
    set(myfig,'PaperPositionMode','manual');
    set(myfig,'PaperUnits','inches');
    set(myfig,'PaperSize',[15 4]);
    set(myfig,'PaperPosition',[.25,.25,14.5,3.5]);
    print(myfig,'-dpdf',[PraatDir,'MinPerformanceAcrossGenerations.pdf']);
    close(myfig);
    
    save([PraatDir,'NeuralNetVocalControlEvolutionWorkspace.mat']);
    
	if (generation ~= 1) && (generation ~= numGenerations)
    	delete([PraatDir,'generation',num2str(generation),'_producerParent*']);
	end
	
end

% Plot the genes of all individuals having a given rank (constantRank) as the population evolves:
if 0
    maxGeneration=100;
    constantRankGenes=[];
    constantRank = 50;
    for generation=1:maxGeneration
        [~,sortInd]=sort(producerParentFitnessesDiary(:,generation),'descend');
        constantRankInd(generation) = sortInd(50,1);
        constantRankGenes=[constantRankGenes,producerParentGenesDiary{generation,1}(constantRankInd(generation),:)'];
    end
    figure; imagesc(constantRankGenes);xlabel('Generation'),ylabel('Gene number'); colormap('gray'); colorbar(); caxis([-5,5]);
end







