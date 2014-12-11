function [producerParentOutputsDiary,aspectra,perceiverParentOutputsDiary,perceiverParentInputsDiary,perceiverParentCorrectness] = getProducerParentSounds(producerParent,producerParentGenes,numProducerNetInputs,numProducerNetHidden,numProducerNetOutputs,numSignals,producerInputs,PraatDir,generation,useVocalTract,timestep,duration,producerParentOutputsDiary,melfcc_wintime,melfcc_nbands,melfcc_hoptime,numIndividuals,perceiverParentGenes,numPerceiverNetInputs,numPerceiverNetHidden,numPerceiverNetOutputs,perceiverTargets,perceiverParentOutputsDiary,perceiverParentInputsDiary,perceiverParentCorrectness,saveSounds,noiseAmnt,maxFreq)

% Anne S. Warlaumont

% get producer's connection weight values from its genome:
producerNetWeightsIH = producerParentGenes(producerParent,1:numProducerNetInputs*numProducerNetHidden);
producerNetWeightsHO = producerParentGenes(producerParent,((numProducerNetInputs*numProducerNetHidden)+1):((numProducerNetInputs*numProducerNetHidden)+(numProducerNetHidden*numProducerNetOutputs)));

% transform weights from vectors into matrices:
producerNetWeightsIH = reshape(producerNetWeightsIH,numProducerNetInputs,numProducerNetHidden);
producerNetWeightsHO = reshape(producerNetWeightsHO,numProducerNetHidden,numProducerNetOutputs);

% Generate the producer's fixed signals
for signalNum = 1:numSignals
    
    % set the input activations
    producerNetInputs = producerInputs(signalNum,:);
    
    aspectrumFilename = [PraatDir,'generation',num2str(generation),'_producerParent',num2str(producerParent),'_vocalization',num2str(signalNum)];
    
    if useVocalTract
        % choose Praat script and vocalization waveform filenames, then open a new Praat script file
        articulationScriptFilenames{signalNum} = [PraatDir,'generation',num2str(generation),'_producerParent',num2str(producerParent),'_articulationScript',num2str(signalNum),'.praat'];
        vocalizationFilenames{signalNum} = [PraatDir,'generation',num2str(generation),'_producerParent',num2str(producerParent),'_vocalization',num2str(signalNum),'.wav'];
        fid = fopen(articulationScriptFilenames{signalNum},'w');
        fprintf(fid,'Create Speaker... speaker Female 2\n');
        %         fprintf(fid,'Create Speaker... speaker Child 2\n');
        %         fprintf(fid,'Create Speaker... speaker Male 2\n');
        fprintf(fid,'Create Artword... vocalization 0.5\n');
        fprintf(fid,'select Artword vocalization\n');
        
        fprintf(fid,['Set target... ',num2str(0,'%.3f'),' ',num2str(1,'%.2f'),' Lungs\n']);
        fprintf(fid,['Set target... ',num2str(duration,'%.3f'),' ',num2str(0,'%.2f'),' Lungs\n']);
    end
    
    for time=0:timestep:duration
        
        producerNetHidden = tanh(producerNetInputs*producerNetWeightsIH);
        producerNetOutputs = satlins(producerNetHidden*producerNetWeightsHO);
        producerParentOutputsDiary{generation,1}=[producerParentOutputsDiary{generation,1};producerNetOutputs];
        
        if useVocalTract
            fprintf(fid,['Set target... ',num2str(time,'%.3f'),' ',num2str(producerNetOutputs(1,1),'%.2f'),' Interarytenoid\n']);
            fprintf(fid,['Set target... ',num2str(time,'%.3f'),' ',num2str(producerNetOutputs(1,2),'%.2f'),' Cricothyroid\n']);
            fprintf(fid,['Set target... ',num2str(time,'%.3f'),' ',num2str(producerNetOutputs(1,3),'%.2f'),' Vocalis\n']);
            fprintf(fid,['Set target... ',num2str(time,'%.3f'),' ',num2str(producerNetOutputs(1,4),'%.2f'),' Thyroarytenoid\n']);
            fprintf(fid,['Set target... ',num2str(time,'%.3f'),' ',num2str(producerNetOutputs(1,5),'%.2f'),' PosteriorCricoarytenoid\n']);
            fprintf(fid,['Set target... ',num2str(time,'%.3f'),' ',num2str(producerNetOutputs(1,6),'%.2f'),' LateralCricoarytenoid\n']);
        end
        
    end
    
    if useVocalTract
        fprintf(fid,'select Speaker speaker\n');
        fprintf(fid,'plus Artword vocalization\n');
        fprintf(fid,'To Sound... 22050 25 0 0 0 0 0 0 0 0 0\n');
        fprintf(fid,'select Sound vocalization_speaker\n');
        %         fprintf(fid,'Play\n');
        fprintf(fid,['\tWrite to WAV file... ',vocalizationFilenames{signalNum},'\n']);
        fclose(fid);
        
        % Run the Praat script:
        [~,~]=system(['/Applications/Praat.app/Contents/MacOS/Praat ',articulationScriptFilenames{signalNum}]);
        
        % Get information about this vocalization
        [sounds{signalNum},samplingRate] = audioread(vocalizationFilenames{signalNum});
        [~,aspectra{signalNum},~]=melfcc(sounds{signalNum},samplingRate,'minfreq',20,'maxfreq',maxFreq,'wintime',melfcc_wintime,'nbands',melfcc_nbands,'hoptime',melfcc_hoptime); 
        aspectra{signalNum} = (aspectra{signalNum})/.5e12; % used in the actual runs
        plot_aspectrum = sqrt(sqrt(aspectra{signalNum})); % used for plotting for paper
        caxisBounds = [0,1];
    else
        aspectra{signalNum} = reshape(producerNetOutputs,3,2);
        aspectra{signalNum} = (aspectra{signalNum}+1)/2;
        plot_aspectrum = aspectra{signalNum};
        caxisBounds = [0,1];
    end
    
    aspectrum_fig = figure('visible','off'); colormap(flipud(gray)); image(flipud(plot_aspectrum),'CDataMapping','scaled'); caxis(caxisBounds); set (gca,'XTick',[]); set(gca,'YTick',[]); print(aspectrum_fig,'-dtiff',aspectrumFilename); close(aspectrum_fig);
    
    % see how well each of the perceivers does at figuring out what the producer is trying to signal
    for perceiverParent = 1:numIndividuals
        
        % get perceiver's connection weight values from its genome:
        perceiverNetWeightsIH = perceiverParentGenes(perceiverParent,1:numPerceiverNetInputs*numPerceiverNetHidden);
        perceiverNetWeightsHO = perceiverParentGenes(perceiverParent,((numPerceiverNetInputs*numPerceiverNetHidden)+1):((numPerceiverNetInputs*numPerceiverNetHidden)+(numPerceiverNetHidden*numPerceiverNetOutputs)));
        
        % transform weights from vectors into matrices:
        perceiverNetWeightsIH = reshape(perceiverNetWeightsIH,numPerceiverNetInputs,numPerceiverNetHidden);
        perceiverNetWeightsHO = reshape(perceiverNetWeightsHO,numPerceiverNetHidden,numPerceiverNetOutputs);
        
        % randomly decide the timestep when the vocalization should start
        % startTimeBin = round(rand*size(aspectra{signalNum},2));
        startTimeBin = 0;
        
        % Get this perceiver parent's outputs:
        perceiverNetInputs = reshape(aspectra{signalNum},size(aspectra{signalNum},1)*size(aspectra{signalNum},2),1)+noiseAmnt*rand(size(aspectra{signalNum},1)*size(aspectra{signalNum},2),1);
        perceiverTarget = perceiverTargets(signalNum,:);
        perceiverNetHidden = tanh(perceiverNetInputs'*perceiverNetWeightsIH);
        perceiverNetOutput = satlin(perceiverNetHidden*perceiverNetWeightsHO);
        perceiverParentOutputsDiary{generation,1} = [perceiverParentOutputsDiary{generation,1};perceiverNetOutput];
        perceiverParentInputsDiary{generation,signalNum} = [perceiverParentInputsDiary{generation,signalNum};perceiverNetInputs'];
        % perceiverParentCorrectness(perceiverParent,producerParent) = perceiverParentCorrectness(perceiverParent,producerParent) + 1/((sqrt(sum((perceiverNetOutput-perceiverTarget).^2)))+1);
        % perceiverParentCorrectness(perceiverParent,producerParent) = perceiverParentCorrectness(perceiverParent,producerParent) +  log(1/((sqrt(sum((perceiverNetOutput-perceiverTarget).^2)))+1)+1);
        [~,perceiverNetWinner] = max(perceiverNetOutput); % What if there was a tie?
        if perceiverNetWinner == signalNum
            perceiverParentCorrectness(perceiverParent,producerParent) = perceiverParentCorrectness(perceiverParent,producerParent) + 1;
        end
        
    end
    
    if ~saveSounds
        delete([PraatDir,'generation',num2str(generation),'_producerParent*']);
    end
    
    % save([PraatDir,'generation',num2str(generation),'perceiverParentInputsDiary.mat'],'perceiverParentInputsDiary');

end