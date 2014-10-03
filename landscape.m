function [fits] = landscape(runDir,gen,PraatDir)
%landscape Plot a fitness landscape for a NeuralNetVocalControlEvolution
%run
%
%   [] = landscape(runDir,gen,saveDir) plots a fitness landscape for the 
%   simulation within runDir. The landscape is defined by taking each gene 
%   from the 50th most fit producer at generation gen and systematically 
%   varying that gene while holding all other genes constant. The fitness 
%   of the individual with those different gene values constitutes the 
%   approximation of the fitness landscape. In fact, the landscape would 
%   fall within a gen-dimensional space, but that would be too hard to 
%   visualize and would take a very long time to run. The image of the
%   fitness landscape is saved within runDir under the name 
%   "fitnessLandscape.eps".
%
%   Examples:
%   
%   fits = landscape('/Volumes/Storage/NeuralNetVocalTractEvolutionRuns/AbstractReady/run20131126T160313/',500,'/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/Temp/');
%
%   /Applications/MATLAB_R2014a.app/bin/matlab nohup -nodisplay -r "landscape('/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/AbstractReady/run20131126T160313/',500,'/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/Temp/')" > "temp.txt" &
%
%   Anne S. Warlaumont

display(['run: ',runDir]);

% Load the necessary variables from the simulation workspace:
load([runDir,'NeuralNetVocalControlEvolutionWorkspace.mat'],'duration','melfcc_hoptime',...
    'melfcc_nbands','melfcc_wintime','numIndividuals',...
    'numPerceiverNetHidden','numPerceiverNetInputs','numPerceiverNetOutputs',...
    'numProducerNetHidden','numProducerNetInputs','numProducerNetOutputs',...
    'numSignals','perceiverParentGenesDiary','perceiverTargets',...
    'producerInputs','producerParentFitnessesDiary','producerParentGenesDiary','timestep','useVocalTract');

% Find the median performing producer:
[~,sortInd]=sort(producerParentFitnessesDiary(:,gen),'descend');
producerParent = sortInd(50,1);
origGenes = producerParentGenesDiary{gen,1}(producerParent,:);

% Get the perceiver genes (needed to evaluate producer fitness):
perGenes = perceiverParentGenesDiary{gen,1};

numGenes = size(origGenes,2); % Get the number of genes.
changes = [-5/7,-2.5/7,0,2.5/7,5/7]; % Set the manipulations that will be made to each gene.
fits = NaN(numGenes,length(changes)); % Initialize a matrix that will hold all the fitnesses of the individual with each gene modified by each amount.

for gene = 1:numGenes % Cycle through the genes and populate the fits matrix.
    for changeind = 1:length(changes)
        change = changes(changeind);
        if change == 0 % if change is 0, the fitness should already be available in producerParentFitnessesDiary
            fits(gene,changeind) = producerParentFitnessesDiary(producerParent,gen);
        else
            tempProducerParentGenes = origGenes;
            tempProducerParentGenes(1,gene) = tempProducerParentGenes(1,gene)+change;
            producerParentGenes(producerParent,:) = tempProducerParentGenes;
            perceiverParentCorrectness = zeros(numIndividuals,numIndividuals);
            perceiverParentInputsDiary = cell(gen,numSignals);
            perceiverParentOutputsDiary = cell(gen,numSignals);
            producerMuscleOutputsDiary = cell(gen,numSignals);
            producerParentOutputsDiary = cell(gen,numSignals);
            [~,~,~,~,~,perceiverParentCorrectness]...
                = getProducerParentSounds(producerParent,producerParentGenes,numProducerNetInputs,...
                numProducerNetHidden,numProducerNetOutputs,numSignals,producerInputs,PraatDir,gen,...
                useVocalTract,timestep,duration,producerParentOutputsDiary,producerMuscleOutputsDiary,...
                melfcc_wintime,melfcc_nbands,melfcc_hoptime,numIndividuals,perGenes,...
                numPerceiverNetInputs,numPerceiverNetHidden,numPerceiverNetOutputs,perceiverTargets,...
                perceiverParentOutputsDiary,perceiverParentInputsDiary,perceiverParentCorrectness,0);
            fits(gene,changeind) =  mean(perceiverParentCorrectness(:,producerParent))+.01;
        end
    end
    display(['Gene ',num2str(gene),' of ',num2str(numGenes)]);
end

landscapefig = figure('visible','off');
colormap(flipud(gray));
imagedata = fits-producerParentFitnessesDiary(producerParent,gen);
%imagedata = padarray(imagedata,[1,1],0,'both');
%imagedata = interp2(imagedata,7,'nearest');
image(imagedata,'CDataMapping','scaled');
caxis([-.5,.5]);
set (gca,'XTick',[]); set(gca,'YTick',[]);
print(gcf,'-dtiff',[runDir,'landscape.tif']);
close(landscapefig);

csvwrite([runDir,'landscape.csv'],fits);
