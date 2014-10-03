% Create figures of some of the sound signals
%
% Anne S. Warlaumont

wavfiles = {'/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160354/generation1_producerParent89_vocalization1.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160354/generation1_producerParent89_vocalization2.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160354/generation1_producerParent89_vocalization3.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160354/generation500_producerParent61_vocalization1.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160354/generation500_producerParent61_vocalization2.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160354/generation500_producerParent61_vocalization3.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160359/generation1_producerParent75_vocalization1.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160359/generation1_producerParent75_vocalization2.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160359/generation1_producerParent75_vocalization3.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160359/generation500_producerParent67_vocalization1.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160359/generation500_producerParent67_vocalization2.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160359/generation500_producerParent67_vocalization3.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160404/generation1_producerParent23_vocalization1.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160404/generation1_producerParent23_vocalization2.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160404/generation1_producerParent23_vocalization3.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160404/generation500_producerParent62_vocalization1.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160404/generation500_producerParent62_vocalization2.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160404/generation500_producerParent62_vocalization3.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160409/generation1_producerParent22_vocalization1.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160409/generation1_producerParent22_vocalization2.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160409/generation1_producerParent22_vocalization3.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160409/generation500_producerParent91_vocalization1.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160409/generation500_producerParent91_vocalization2.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160409/generation500_producerParent91_vocalization3.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160414/generation1_producerParent85_vocalization1.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160414/generation1_producerParent85_vocalization2.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160414/generation1_producerParent85_vocalization3.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160414/generation500_producerParent85_vocalization1.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160414/generation500_producerParent85_vocalization2.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160414/generation500_producerParent85_vocalization3.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160425/generation1_producerParent88_vocalization1.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160425/generation1_producerParent88_vocalization2.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160425/generation1_producerParent88_vocalization3.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160425/generation500_producerParent5_vocalization1.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160425/generation500_producerParent5_vocalization2.wav',...
    '/Volumes/Storage/NeuralNetVocalControlEvolutionRuns/RealisticReady/run20131126T160425/generation500_producerParent5_vocalization3.wav',...
    };

for wavnum = 1:length(wavfiles)
    
    [vocwav,sr] = audioread(wavfiles{wavnum});
	waveformfile = regexprep(wavfiles{wavnum},'.wav','_waveform.eps');
    specfile = regexprep(wavfiles{wavnum},'.wav','_spec.eps');
    melspecfile = regexprep(wavfiles{wavnum},'.wav','_melspec.tif');
    
	waveformfig = figure('visible','off');
	plot(vocwav,'k');
	xlim([1,length(vocwav)]);
	set(gca,'xtick',[1,length(vocwav)]);
	set(gca,'xticklabel',[0,length(vocwav)/sr])
	ylim([-1,1]);
	xlabel('Time (s)');
	ylabel('Amplitude');
	print(waveformfig,'-deps',waveformfile);
	close(waveformfig);
	
    specfig = figure('visible','off');
    specgram(vocwav,512,sr);
    colormap(flipud(gray));
    caxis([-50,100]);
    print(specfig,'-deps',specfile);
    close(specfig);
    
    melspecfig = figure('visible','off');
    [~,vocspec,~] = melfcc(vocwav,sr,'minfreq',20,'maxfreq',2000);
    image(flipud(vocspec/.5e12),'CDataMapping','scaled');
    colormap(flipud(gray));
    caxis([0,.00001]);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    print(melspecfig,'-dtiff',melspecfile);
    close(melspecfig);
    
end
