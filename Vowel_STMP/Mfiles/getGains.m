function gains = getGains(thisDir,gains)

freq_Hz = 1700;

if nargin<1
    thisDir = 'C:\Research\MATLAB\Vowel_STMP\ExpData';
end
if nargin<2
    gains.linear = [];
    gains.nonlinear = [];
end

% get gains based on the file names in this dir
files = dir(fullfile(thisDir,'*_linearAid_*'));
if ~isempty(files)
    files_lin = arrayfun(@(x) x.name,files,'UniformOutput',false);
    ind = strfind(files_lin,'dBEH');
    for ii=1:numel(files_lin)
        gains.linear = [gains.linear getLinearGain(freq_Hz)];
    end
end

files = dir(fullfile(thisDir,'*_nonlinearAid_*'));
if ~isempty(files)
    files_nonlin = arrayfun(@(x) x.name,files,'UniformOutput',false);
    ind = strfind(files_nonlin,'dBEH');
    for ii=1:numel(files_nonlin)
        atten = str2double(files_nonlin{ii}(ind{ii}+(-3:-1)));
        gains.nonlinear = [gains.nonlinear getNonLinearGain(atten,freq_Hz)];
    end
end

% now, do subdirs
listing = dir(thisDir);
subDirs = arrayfun(@(x) isdir(fullfile(thisDir,x.name)),listing);
subDirs = listing(subDirs);
for ii=3:numel(subDirs)
    gains = getGains(fullfile(thisDir,subDirs(ii).name),gains);
end
end



function gain = getLinearGain(freq_Hz)
audiogram = [16   18   20    9    9]; % based on avg animal data (500OBN exposure)
freqs_Hz = [500,1000,2000,4000,6000];
freqs_Hz_std=[250,500,1000,2000,3000,4000,6000];

dBLoss = abs(audiogram);
% interpolate/extrapolate audiogram at standard frequencies
dBLoss = interp1(freqs_Hz,dBLoss,freqs_Hz_std,'nearest','extrap');

k_NAL = [-17 -8 1 -1 -2 -2 -2]; % dB
H_3FA = sum(dBLoss(2:4)) / 3; % sum up loss at 500Hz,1kHz,2kHz
X = 0.15*H_3FA;
R = 0.31; % NAL-R formula  (not quite half-gain rule)
NAL_IG = X + R.*dBLoss + k_NAL; % insertion gain

gain = interp1(freqs_Hz_std,NAL_IG,freq_Hz); % gain @ freq_Hz
end

function gain = getNonLinearGain(atten,freq_Hz)
MaxSPL = 105;
audiogram = [16   18   20    9    9]; % based on avg animal data (500OBN exposure)
freqs_Hz = [500,1000,2000,4000,6000];
InputRMS = 0.4259;
InputSPL = (MaxSPL-atten) + 20*log10(InputRMS/(1/sqrt(2)));

[Gain, Thresh, Ratio, TargetREAR_60dBSPL] = readDSLfile('ChinchillaDSLtargets.csv');
indx=1;

% Initialize gains based on compression parameters
PrescribedGain = Gain(indx,:) + ...
    max(0,(InputSPL-Thresh(indx,:))).*Ratio(indx,:) - ...
    max(0,(InputSPL-Thresh(indx,:)));
LTASS_RMS = -12.49; % avg rms(dBFS) of LTASS
LTASS_dBFS_sb = [-12.55, -33.21]; % avg subband levels (dBFS), fc=2.1kHz
LTASS_dBFS_to = [-24.08, -42.37]; % 1/3-octave levels at 500Hz & 4kHz
LTASS_dBSPL_sb = 60 + (LTASS_dBFS_sb-LTASS_RMS); % get subband SPL (total RMS=60dB SPL)
LTASS_dBSPL_to = 60 + (LTASS_dBFS_to-LTASS_RMS);
PrescribedGain2 = Gain(indx,:) + ...
    max(0,(LTASS_dBSPL_sb-Thresh(indx,:))).*Ratio(indx,:) - ...
    max(0,(LTASS_dBSPL_sb-Thresh(indx,:))); % get prescribed gain for LTASS
%PrescribedGain2 = max(0,PrescribedGain2); % minimum 0dB gain
PredictedREAR = LTASS_dBSPL_to + PrescribedGain2;
adjustment = TargetREAR_60dBSPL(indx,:) - PredictedREAR;

PrescribedGain = PrescribedGain + adjustment; % to reach REAR
gain = PrescribedGain(2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [TKgain_out,TK_out,CR_out,Target_out] = readDSLfile(filename)  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TKgain_out,TK_out,CR_out,Target_out] = readDSLfile(filename)
a = csvread(filename,0,1);
length = size(a,2);

% Program 1 (in quiet)
freqs(1,:) = a(1,1:length-1);
% Line 2: thresholds entered into the DSL m(I/O) v5.0a GUI
Thresh(1,:) = a(2,1:length-1);
% (note that line 9 was empty, so didn't get read)
% Line 10: thresholds in dB SPL (Line 2 will be the same if real ear SPL thresholds are entered)
ThreshSPL(1,:) = a(9,1:length-1);
% Line 15: Output at each compression threshold
CTout(1,:) = a(14,1:length-1);
% Line 17:  Channelized input at each compression threshold
TK(1,:) = a(16,1:length-1);
% Line 12: BOLT (Broadband Output Limiting Targets)
BOLT(1,:) = a(11,1:length-1);
% Line 19: Compression ratios
CR(1,:) = a(18,1:length-1);
% [The difference between lines 17 and 15 is used to initialize gain]
TKgain(1,:) = CTout(1,:) - TK(1,:);
% Lines 25: Real-ear aided response for low speech (52 dB SPL)
TargetLo(1,:) = a(23,1:length-1);
% Lines 25: Real-ear aided response for average speech (60 dB SPL)
TargetAvg(1,:) = a(24,1:length-1);
% Lines 25: Real-ear aided response for high speech (74 dB SPL)
TargetHi(1,:) = a(25,1:length-1);

% Program 2 (in noise)
freqs(2,:) = a(1,1:length-1);
% Line 2: thresholds entered into the DSL m(I/O) v5.0a GUI
Thresh(2,:) = a(2,1:length-1);
% (note that line 9 was empty, so didn't get read)
% Line 10: thresholds in dB SPL (Line 2 will be the same if real ear SPL thresholds are entered)
ThreshSPL(2,:) = a(9+23,1:length-1);
% Line 15: Output at each compression threshold
CTout(2,:) = a(14+23,1:length-1);
% Line 17:  Channelized input at each compression threshold
TK(2,:) = a(16+23,1:length-1);
% Line 12: BOLT (Broadband Output Limiting Targets)
BOLT(2,:) = a(11+23,1:length-1);
% Line 19: Compression ratios
CR(2,:) = a(18+23,1:length-1);
% [The difference between lines 17 and 15 is used to initialize gain]
TKgain(2,:) = CTout(2,:) - TK(2,:);
% Lines 25: Real-ear aided response for low speech (52 dB SPL)
TargetLo(2,:) = a(23+23,1:length-1);
% Lines 25: Real-ear aided response for average speech (60 dB SPL)
TargetAvg(2,:) = a(24+23,1:length-1);
% Lines 25: Real-ear aided response for high speech (74 dB SPL)
TargetHi(2,:) = a(25+23,1:length-1);


indices = [5,14]; %[500,4000]Hz
% Compression Parameters
% [in_quiet_band1,in_quiet_band2;
%  in_noise_band1,in_noise_band2;]
TKgain_out = [TKgain(1,indices(1)),TKgain(1,indices(end));...
    TKgain(2,indices(1)),TKgain(2,indices(end))];
TK_out = [TK(1,indices(1)),TK(1,indices(end));...
    TK(2,indices(1)),TK(2,indices(end))];
CR_out = [CR(1,indices(1)),CR(1,indices(end));...
    CR(2,indices(1)),CR(2,indices(end))];
% REAR Targets
Target_out = [TargetAvg(1,indices(1)),TargetAvg(1,indices(end));...
    TargetAvg(2,indices(1)),TargetAvg(2,indices(end))];

end %function
