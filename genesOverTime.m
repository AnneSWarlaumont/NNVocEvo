function [] = genesOverTime(runsFolder,workspaceFolder,numgens,varargin)

% Plot the Nth most fit (i.e the constantRank) individual's genes over the course of numgens generations of evolution.
% Example: genesOverTime('/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/Abstract_CapWeights/','run1/',400,'constantRank',50);
%
% Anne S. Warlaumont

p = inputParser;
addRequired(p,'runsFolder');
addRequired(p,'workspaceFolder');
addRequired(p,'numgens');
addParameter(p,'constantRank',50);

parse(p,runsFolder,workspaceFolder,numgens,varargin{:});

constantRank = p.Results.constantRank;

load([runsFolder,workspaceFolder,'NeuralNetVocalControlEvolutionWorkspace.mat'],...
    'producerParentFitnessesDiary','producerParentGenesDiary',...
    'perceiverParentFitnessesDiary','perceiverParentGenesDiary');

% Initialize the matrix that will hold the genes for the constantRank fitness individual over time:
constantRankProGenes=[];
constantRankPerGenes=[];

% Set the file names that the genes over time figures will be written to:
constantFitProGenesFilename = [runsFolder,workspaceFolder,'fitRank',num2str(constantRank),'_producergenes_gen1to',num2str(numgens),'.tif'];
constantFitPerGenesFilename = [runsFolder,workspaceFolder,'fitRank',num2str(constantRank),'_perceivergenes_gen1to',num2str(numgens),'.tif'];

% Get the genes over time for the producers:
for generationcount=1:numgens % Loop over each generation
    [~,sortInd]=sort(producerParentFitnessesDiary(:,generationcount),'descend'); % find out the relative fitnesses of each individual
    constantRankInd(generationcount) = sortInd(constantRank,1); % find the index of the individual with the desired relative fitness
    constantRankProGenes=[constantRankProGenes,producerParentGenesDiary{generationcount,1}(constantRankInd(generationcount),:)']; % Append the individuals' genes to the record of genes across generations
end

% Plot the producer genes over time:
constantFitProGenes_fig = figure('visible','off');
imagesc(constantRankProGenes);xlabel('Generation'),ylabel('Gene number'); colormap(flipud(gray)); colorbar();
print(constantFitProGenes_fig,'-dtiff',constantFitProGenesFilename); close(constantFitProGenes_fig);

% Get the genes over time for the perceivers:
for generationcount=1:numgens
    [~,sortInd]=sort(perceiverParentFitnessesDiary(:,generationcount),'descend'); % find out the relative fitnesses of each individual
    constantRankInd(generationcount) = sortInd(50,1); % find the index of the individual with the desired relative fitness
    constantRankPerGenes=[constantRankPerGenes,perceiverParentGenesDiary{generationcount,1}(constantRankInd(generationcount),:)']; % Append the individuals' genes to the record of genes across generations
end

% Plot the perceiver genes over time:
constantFitPerGenes_fig = figure('visible','off');
imagesc(constantRankPerGenes);xlabel('Generation'),ylabel('Gene number'); colormap(flipud(gray)); colorbar();
print(constantFitPerGenes_fig,'-dtiff',constantFitPerGenesFilename); close(constantFitPerGenes_fig);
