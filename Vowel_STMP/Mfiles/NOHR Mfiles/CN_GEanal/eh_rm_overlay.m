function eh_rm_overlay(trackNum, unitNum, rm_pic, data_folder)
% Written by GE 29Mar2004.


global uinfo FIG Data;

if ~exist('data_folder')
   uinfo.dataFolder = [cd '\'];
else
   uinfo.dataFolder = data_folder;
end
uinfo.TNum = trackNum;
uinfo.UNum = unitNum;
uinfo.unitName = sprintf('%d.%02d', trackNum, unitNum);

error = get_unit_info;
if (error)
   break;
end
layout_figure;

overlay_eh_spectra;

%%##################################################################################################
function error = get_unit_info()
global uinfo Data;

error = 0;

% Temporarily switch to experiment data directory and load "info file":
origDir = cd; eval (['cd ''' uinfo.dataFolder '''']); load DataInfoFile; eval (['cd ''' origDir '''']);

% Setup appropriate parameters based on loaded unit info data:

     % General info:
     
if (uinfo.TNum > size(Data.Info, 1) |  uinfo.UNum > size(Data.Info, 2))
	error_report = 'Invalid track or unit number.'
	error = -1;
	return
end
infoField = Data.Info{uinfo.TNum, uinfo.UNum};
if (isempty(infoField))
   error_report = 'Invalid track or unit number.'
   error = -1;
   return
end
uinfo.text = sprintf('UNIT %s (experiment date %s)  \n%s %s unit: BF = %.1fkHz, threshold = -%ddB ', ...
               uinfo.unitName, ...
               strrep(Data.General.date,'_', '\_'), ...
               infoField.Location, ...
               infoField.UnitType, ...
               infoField.BF_kHz, ...
               infoField.Threshold_dBatten ...
             );
uinfo.BF_kHz = infoField.BF_kHz;
          
     % /eh/ vowel smp info:
uinfo.SMP.EH = Data.Vowels_SMP{uinfo.TNum,uinfo.UNum}.EH;


%%##################################################################################################
function layout_figure()
global uinfo FIG;

FIG.fontSize = 8;
FIG.handles.main = figure(4); clf; hold on;  % ge debug 19Mar2004, figure number hard-coded, for now
set(gcf, 'Name', 'modulated /eh/ spectra (eh_rm_overlay.m)');

axisLeft = 0.1; axisBottom = 0.1; axisWidth = 0.8; axisHeight = 0.8;
FIG.handles.eh_overlay = subplot('Position',[axisLeft axisBottom axisWidth axisHeight]);
hold on;
set(gca, 'FontSize', FIG.fontSize);
set(gca, 'YAxisLocation', 'left');
set(gca, 'XTickLabelMode', 'auto');
YLabel('dB atten for response map');
XLabel('freq (kHz)');
% set(gca, 'XScale', 'log');
freqTicks = [.1 .2 .3 .5 1 2 3 5 10 20 30 50 100];
set(gca, 'XTick', 1000*freqTicks, 'XTickLabel', freqTicks, 'XMinorTick', 'off');
% set(gca, 'YTickLabel', '');

FIG.handles.figBox = axes('Position',[0 0 1 1],'Visible','off');
text(0.02, 0.95, uinfo.text, 'Units', 'normalized', 'FontSize', FIG.fontSize-1);


%%##################################################################################################
function error = overlay_eh_spectra()
global uinfo Data FIG;


   
% Parameters specific for vowel /eh/:
feature_names = {'T0' 'F1' 'T1' 'F2' 'T2' 'F3' 'T3'};
feature_base_freq_Hz = [300 501 1202 1703 2203 2504 3005];
feature_levels_dB = [-11.2 0 -30.6 -15.5 -33.6 -28.7 -42.7];
feature_line_colors = ['k' 'r' 'r' 'g' 'g' 'b' 'b'];


subplot(FIG.handles.eh_overlay);

% Plot the shifted vowel spectra on copies of the response map for the unit:
origDir = cd; 	eval (['cd ''' uinfo.dataFolder '\Signals\GE_MH']); 	[eh, fs] = wavread('vow17_MH10k.wav');	eval (['cd ''' origDir '''']);
eh_fft = abs(fft(eh));
eh_fft = eh_fft(1:floor(length(eh_fft)/2));
eh_fft = eh_fft / max(eh_fft);
eh_freq = fs*[1:length(eh_fft)]/length(eh);
% figure(4); plot(eh_freq, 20*log10(eh_fft/max(eh_fft)));
% eh_fft(find(eh_fft<10^(-1.75))) = NaN;
figure(4); plot(eh_freq, 20*log10(eh_fft));
[eh_fft_peaks eh_fft_peaks_f] = get_peaks(eh_fft, eh_freq);
% figure(4); plot(eh_freq, eh_fft);
figure(4); plot(eh_fft_peaks_f, 20*log10(eh_fft_peaks), 'r');
for feature_index = 1:7
   
end

%%##################################################################################################
function [out1, out2] = get_peaks(in1, in2)

diffFFT = diff(in1);
out1 = in1(find(diffFFT > 0.005)+1);
out2 = in2(find(diffFFT > 0.005));
