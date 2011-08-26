function runANmodel(ExpDate,picNum)
% File: runANmodel.m
% 10Feb2005 (M Heinz)
% Takes a saved picture from NEL (see NELsoftware/SC091704) to get desired params,
% creates appropriate stimulus for AN model, runs AN model to get X reps, then substitutes
% AN model data into picture and saves new and original (e.g., op0005_EHreBF.m)
% Assumes EH_reBF, so that there is only 1 condition per file.

global NOHR_dir NOHR_ExpList
global FeaturesText FormsAtHarmonicsText InvertPolarityText
global x

ANmodel_dir= 'C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\ANModel\ARLO';

path(ANmodel_dir,path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify parameters and experiment, unit are valid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Verify in passed parameters if needed
if ~exist('ExpDate','var')
   %%% HARD CODE FOR NOW
   ExpDate='021005'
   error('No ExpDate')
end
if ~exist('picNum','var')
   error('No picNum')
end

%%%% Find the full Experiment Name 
ExpDateText=strcat('20',ExpDate(end-1:end),'_',ExpDate(1:2),'_',ExpDate(3:4));
for i=1:length(NOHR_ExpList)
   if ~isempty(strfind(NOHR_ExpList{i},ExpDateText))
      ExpName=NOHR_ExpList{i};
      break;
   end
end
if ~exist('ExpName','var')
   disp(sprintf('***ERROR***:  Experiment: %s not found\n   Experiment List:',ExpDate))
   disp(strvcat(NOHR_ExpList))
   beep
   break
end

%%%% Parse out the Track and Unit Number 
[temp]=getTrackUnit(getFileName(picNum));
TrackNum=temp(1);
UnitNum=temp(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_dir=fullfile(NOHR_dir,'ExpData');
anal_dir=fullfile(NOHR_dir,'Data Analysis');
stim_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\NELsoftware\SC091704\Signals\MH\EHvowels';

eval(['cd ''' fullfile(data_dir,ExpName) ''''])
disp(sprintf('Running AN model for picture %d (Experiment: ''%s''; Unit: %d.%02d)',picNum,ExpName,TrackNum,UnitNum))

FileName=strcat('DataList_',ExpDateText,'.mat');
disp(['   *** Loading file: "' FileName '" to get unit BF'])
eval(['load ' FileName])

% Load picture
x=loadPic(picNum);

% Load stimulus, and create original BASELINE wavefile - full duration
[pathstr,name,ext,versn]=fileparts(x.Stimuli.BASELINE.FileName);
stimfilename=strcat(name,ext);
[vowel,BASELINE_Fs_Hz] = wavread(fullfile(stim_dir,stimfilename));

% Invert Stimulus Polarity if desired (MH 07July2004)
if strcmp(x.Stimuli.Condition.InvertPolarity,'yes')
   vowel=-vowel;
end

% Repeat or truncate waveform to fit requested stimulus duration:
stim = refit_waveform(vowel,x.Stimuli.updateRate_Hz);

% Window waveform using linear rise/fall:
stim = window_waveform(stim,x.Stimuli.updateRate_Hz);

% Rescale data to Pascals for ANmodel 
stim = rescale_waveform(stim);

% % % % Get UpdateRate for this index
% % % ur_Hz = static_bi.updateRate_Hz(ind);


% RUN MODEL
tdres_sec=1/x.Stimuli.updateRate_Hz;  
cf_Hz=DataList.Units{TrackNum,UnitNum}.Info.BF_kHz*1000;
SR_sps=50;
model=1; %1 is normal
species=9; % 9 is cat
ifspike=1; % 1 means yes, get spikes


disp(sprintf('--- BF = %.4f kHz;  FeatureTargetFreq = %.4f kHz;  Urate = %.2f Hz',cf_Hz/1000,x.Stimuli.featureTargetFreq_Hz/1000,1/tdres_sec))

sout = an_arlo([tdres_sec,cf_Hz,SR_sps,model,species,ifspike],stim'); % cat
nrep = 30;
[sptime,nspikes] = sgmodel([tdres_sec, nrep],sout);
mysptime = sptime(1:nspikes);
NELspikes=NaN*ones(nspikes,2);
NELspikes(:,2)=mysptime;
REPendINDs(1)=0;
REPendINDs(2:nrep)=find(diff(mysptime)<0);
REPendINDs(nrep+1)=nspikes;

% Format spikes into NEL spikes format, and plug into x data struct
for REPind=1:nrep
   NELspikes(REPendINDs(REPind)+1:REPendINDs(REPind+1),1)=REPind;
end





%%%%%%%%% Re-write NEL data
% Save as p0005....., and save original as op0005....
pfileName=getfileName(picNum);
ofileName=pfileName;
ofileName(1)='o';
if ~isempty(dir(ofileName))
   uiwait(warndlg('ofileName already exists','non-modal'))
end
eval(['!copy ' pfileName ' ' ofileName])

% Fix x.spikes,fully_presented_lines,fully_presented_Stim
xORIG=x;
x=rmfield(x,'spikes');  % Take out old spikes data
x.Stimuli.fully_presented_stimuli=nrep;
x.Stimuli.fully_presented_lines=nrep;
global spikes
spikes.times{1}=NELspikes;  % This seems to be what NEL wants
spikes.last=nspikes;

% Check at beginning for presence of op, and if so, don't run this because its been run before
fname=fullfile(data_dir,ExpName,pfileName);
rc = write_nel_data(fname,x,1);

if rc<0
   beep
   disp('Error writing data file')
end

   
rmpath(ANmodel_dir)

return;


%##########################################################################################
function out_wave = refit_waveform(in_wave,ur_Hz)
   global x
   stimOn_sec = x.Hardware.Trigger.StmOn / 1000;
   required_nSamples = floor(ur_Hz*stimOn_sec);
   current_nSamples = size(in_wave, 1);
   if (current_nSamples > required_nSamples) % Waveform needs to be truncated.
      out_wave = in_wave(1:required_nSamples);
      
      disp(['In function ''runANmodel'': Input waveform has been truncated ' ...
            'to fit requested duration.']);
   elseif (current_nSamples < required_nSamples) % Waveform needs to be "repeated".
      nRepeats = ceil(required_nSamples / current_nSamples);
      out_wave = repmat(in_wave, nRepeats, 1);
      out_wave = out_wave(1:required_nSamples); % truncate any extra of the final repeat.
      disp(['In function ''runANmodel'': Input waveform has been repeated ' ...
                 'to fill requested duration.']);
   else
      out_wave = in_wave;
   end


%##########################################################################################
function out_wave = window_waveform(in_wave,ur_Hz)
   global x 
   stimOn_sec = x.Hardware.Trigger.StmOn / 1000;
   rf_time_sec = default_rise_time(stimOn_sec * 1000)/1000;
   ramp_nSamples = ceil(ur_Hz * rf_time_sec);
   end_sample = size(in_wave, 1);
   windowing_vector = zeros(size(in_wave));
   windowing_vector(1:ramp_nSamples) = [0:(ramp_nSamples-1)]/ramp_nSamples;
   windowing_vector((ramp_nSamples+1):(end_sample-ramp_nSamples)) = 1;
   windowing_vector(((end_sample-ramp_nSamples)+1):end) = [(ramp_nSamples-1):-1:0]/ramp_nSamples;
   out_wave = windowing_vector .* in_wave;
   return;
   
%##########################################################################################
function out_wave = rescale_waveform(in_wave)
   global x
   % This assumes WAV file input, e.g. peak =1, and uses dBreTONE for this vowel to scale to Pascals
   VowelWAV_dBSPL = 20*log10(sqrt(2)/2/20e-6) + x.Stimuli.BASELINE.dBreTONE;
   scale_factor = (10^(x.Stimuli.Condition.Level_dBSPL/20))/(10^(VowelWAV_dBSPL/20));  
   out_wave = scale_factor * in_wave;
   return;