function dBAtt_2_SNR = plot_EHINrlfs(RLFpics,CALpic)
% File: plot_EHINrlfs.m
% Modified from: quick_EHINrlfs.m (to clean up for ISH)
% Modified (MHeinz 18Apr2005) to plot vs noise Attenuation
% M.Heinz: 11Nov2004 (From GE: quick_vowel)
% For NOHR
%
% USAGE:quick_vowel(RLFpics,CALpic)  [e.g., quick_vowel([2 3],1)]
% Plots RLFs for a set of vowel features
%
% RLFpics: vector of vowel RLFs for different features [F1_pic, T1_pic]
% CALpic: calib file (e.g., 1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO (4/27/06) _ ADD excludelines from PIC %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TextFontSize=18;
AnnotFonotSize=16;
LineWidth=3;

params = [];
params.spikeChan = 1;
params.figNum = 101;
params.TriFiltWidth=5;
params.colors={'b','r','g','k','m','c','y','b','r','g','k','m','c','y'};

Npicts=length(RLFpics);
x=cell(1,Npicts);
driv=cell(1,Npicts);
spont=cell(1,Npicts);
legtext=cell(1,Npicts);

figure(params.figNum); clf
if length(RLFpics)==1
	nameStr = sprintf ('picture p%04d rate-level', RLFpics);
else
	nameStr = sprintf ('pictures p%04d - p%04d rate-level', RLFpics(1), RLFpics(end));
end
set(gcf, 'Name', nameStr,'pos',[838   540   560   420]);

CAL=loadpic(CALpic);
CAL.CalibData(:,2)=trifilt(CAL.CalibData(:,2)',5)';  % Do some smoothing

for PICind=1:Npicts
	x{PICind} = loadPic(RLFpics(PICind));
	if isfield(x{PICind}.Stimuli,'feature')
		legtext{PICind}=x{PICind}.Stimuli.feature;
	else
		legtext{PICind}='TN';
	end
	[driv{PICind} spont{PICind}] = compute_driv_spont_rates(x{PICind}, params);

	if isfield(x{PICind}.Stimuli,'feature')
		if strcmp(x{PICind}.Stimuli.feature(1),'F')
			linetype=strcat(params.colors{PICind},'-');
		else
			linetype=strcat(params.colors{PICind},'--');
		end
	else
		linetype=strcat(params.colors{PICind},':');
	end

	level_dBSPL(PICind)=x{PICind}.Stimuli.Condition.Level_dBSPL;
	noiseAtten_dB=x{PICind}.Line.attens.Rlist(:,2);

	%% Calibrate Noise Level (as Noise energy out to T3)
	if isfield(x{1}.Stimuli,'BASELINE')

		maxTONE_dBSPL=CalibInterp(x{1}.Stimuli.BASELINE.TargetFreq_Hz/1000,CAL.CalibData(:,1:2));
		dBreTONE_noise = x{1}.Stimuli.BASELINE.dBreTONE_noise;
		BWcalibNoise_Hz = x{1}.Stimuli.BASELINE.FeatFreqs_Hz(7);  % Take noise power out to T3 frequency

		OALnoise_dBSPL = maxTONE_dBSPL + dBreTONE_noise - noiseAtten_dB;

		NoiseSpL_dBSPL = OALnoise_dBSPL - 10*log10(x{1}.Stimuli.BASELINE.Fs_Hz_noise/2);
		CalibNoiseLevels_dBSPL = NoiseSpL_dBSPL + 10*log10(BWcalibNoise_Hz);

		SNR_dB = level_dBSPL(PICind) - CalibNoiseLevels_dBSPL;
	else
		beep
		error(sprintf('Problem with NOISE Calibration '))
	end

	dBAtt_2_SNR = -(noiseAtten_dB(1) - SNR_dB(1)); % Conversion from dBatten to SNR
	plot(SNR_dB(1:length(driv{PICind})), trifilt(driv{PICind},params.TriFiltWidth),linetype,'LineWidth',LineWidth);
	hold on
end

xlabel('Signal to Noise Ratio (dB)','FontSize',TextFontSize)
ylabel('Driven rate (sp/sec)','FontSize',TextFontSize)
set(gca,'FontSize',AnnotFonotSize)
% legend(legtext,4)
% title(nameStr)
set(gca,'XDir','reverse')
text(.05,.25,sprintf('Vowel Level = %.f dB SPL',x{1}.Stimuli.Condition.Level_dBSPL),'HorizontalAlignment','left','units','norm','FontSize',TextFontSize)
text(.94,.94,'B','FontSize',24,'units','norm')

%%#######################################################################
function [driv, spont] = compute_driv_spont_rates(x, params)
% compute driven and spontaneous rates for each line
spikeTimes = x.spikes{params.spikeChan};
driven_dur_sec = x.Hardware.Trigger.StmOn / 1000;
spont_dur_sec = x.Hardware.Trigger.StmOff / 1000;
for line = 1:x.Stimuli.fully_presented_lines
	spikeIndices = find( (spikeTimes(:,1) == line ) );
	driv(line) = length (find(spikeTimes(spikeIndices,2) <= driven_dur_sec));
	spont(line) = length(spikeIndices) - driv(line);
end
driv = driv / driven_dur_sec;
spont = spont / spont_dur_sec;
