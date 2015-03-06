function [areLoud,loudMusc,quietMusc] = loudMuscVals(runDir,generation,thresh)

workspaceFile = strcat(runDir,'NeuralNetVocalControlEvolutionWorkspace.mat');
load(workspaceFile,'numIndividuals','parentSig1LoudnessDiary','parentSig2LoudnessDiary','parentSig3LoudnessDiary','producerParentOutputsDiary');
loudnesses = [parentSig1LoudnessDiary(:,generation); parentSig2LoudnessDiary(:,generation); parentSig3LoudnessDiary(:,generation)];
muscles = [producerParentOutputsDiary{generation,1}(1:6:end,:); producerParentOutputsDiary{generation,1}(3:6:end,:); producerParentOutputsDiary{generation,1}(5:6:end,:)];

for sigNum=1:3

	areLoud = NaN(numIndividuals*3,1);
	loudMusc = [];
	quietMusc = [];
	
	for ind=numIndividuals*(sigNum-1)+1:numIndividuals*sigNum
	    if mean(loudnesses(ind,:)) > thresh
	        areLoud(ind) = 1;
			loudMusc = [loudMusc;muscles(ind,:)];
	    else
	        areLoud(ind) = 0;
			quietMusc = [quietMusc;muscles(ind,:)];
	    end
	end

	close all;
	figure(1); hold on; xlim([0.5,6.5]); ylim([-1,1]);
	for musNum = 1:6
	    plot(musNum*ones(size(quietMusc,1),1),quietMusc(:,musNum),'x','Color','cyan');
	    plot(musNum*ones(size(loudMusc,1),1),loudMusc(:,musNum),'o','Color','black');
	end
	ca = gca;
	ca.XTickLabel = {'Interarytenoid','Cricothyroid','Vocalis','Thyroarytenoid','Posterior Cricoarytenoid','Lateral Cricoarytenoid'};
	rotateXLabels(ca,45);
	ca.FontSize = 14;
	ylabel('Muscle activation','FontSize',16);
	saveas(1,strcat(runDir,'musVals_gen',num2str(generation),'_sigNum',num2str(sigNum),'.png'));
	legend('quiet sounds','loud sounds','location','eastoutside');
	saveas(1,strcat(runDir,'musVals_gen',num2str(generation),'_sigNum',num2str(sigNum),'_legend','.png'));
	close(1);

end