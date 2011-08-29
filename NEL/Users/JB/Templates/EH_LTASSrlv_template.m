function [tmplt,DAL,stimulus_vals,units,errstr] = EH_LTASSrlv_template(fieldname,stimulus_vals,units)
% MH 10Nov2004: for NOHR project
%
% Modified 11-Apr-2005 to add background noise
%   - this is rate vs. Noise level with a fixed EH level in dB SPL
%
%    Rate-level function for features at BF for a BASELINE EH with F2 at BF and F0=75 Hz
%
% From EH_reBF_template and EH_template (CNexps)
%   - kept same params as for EH_reBF, but hardwired OCTshift to 0
%
% Template for EH-vowel RLFs, using NI board to allow resampling
%
% MH 07Nov2003, modified from 'nel_NI_wavefile_template'
%
%  Written by GE, adapted from 'nel_wavefile_template' written by AF (11/26/01).
%  Modification dates: 03oct2003, 04oct2003.

% persistent   prev_playdur  prev_min_period  % GE debug: not really necessary???
% We use the persistent variables to detect a change that requires some fields update.
% For example, if thplay duration is changed we would like to update the gating information.
% We restict the automatic updates to allow the user the overide them.

used_devices.Llist         = 'L3';
used_devices.Rlist         = 'R3';
tmplt = template_definition(fieldname);
if (exist('stimulus_vals','var') == 1)

   %%% Make sure the needed BASELINE-EH wavfile exists.  If not, create it before proceeding
   %%%  BASELINE-EH has F2 at BaseFreq and F0=75Hz;
   %%%  wavfile is of 1 cycle of the vowel
   global signals_dir

   % MH 17Nov2004
   % If no BF is defined, e.g., on startup, need to catch this before all the calcs get done
   if isempty(stimulus_vals.Inloop.BaseFrequency)
      stimulus_vals.Inloop.BaseFrequency=.01;  % Will get through, and let you change later
   end
   
   EHsignals_dir=strcat(signals_dir,'JB\EHvowels');
   BASELINE_TargetFreq_Hz=stimulus_vals.Inloop.BaseFrequency*1000;
   BASELINE_F0_Hz=100;
%    % MH: 12Nov2004 - keeps basic shape the same for all BFs, meeting our constraints
%    if BASELINE_TargetFreq_Hz>=500
%       BASELINE_F0_Hz=75;  % BFs: >=500
%    elseif (BASELINE_TargetFreq_Hz>=300)&(BASELINE_TargetFreq_Hz<500)
%       BASELINE_F0_Hz=45;  % BFs:300-500
%    elseif BASELINE_TargetFreq_Hz<300
%       BASELINE_F0_Hz=15;  % BFs:100-300
%    else
%       BASELINE_F0_Hz=75;  % in case BASELINE_TargetFreq_Hz is empty or something, e.g, on startup
%    end
   BASELINE_Feature='F2';
   Fix2Harms=strcmp(stimulus_vals.Inloop.FormsAtHarmonics,'yes'); % Set formants at nearest harmonic: 0:no, 1:yes
   PolarityFact=(strcmp(stimulus_vals.Inloop.InvertPolarity,'yes')-.5)*-2;  % For inverting waveform if necessary
   % Get filename and FormFreqs: mode 3 returns empty stim,Fs,dBreTONE
   [Xstim,XFs,filename,XdBreTONE,BASELINE_FormFreqs_Hz]= ...
      synth_BASELINE_eh(BASELINE_TargetFreq_Hz,BASELINE_F0_Hz,BASELINE_Feature,Fix2Harms,3);
   %%% MH: 11Nov2004 Check here for problems with stimulus design (e.g., formants at 0 freq, or equal to other formants 
   if length(unique(BASELINE_FormFreqs_Hz))~=length(BASELINE_FormFreqs_Hz)|sum(BASELINE_FormFreqs_Hz==0)
      BADstim=1;
   else
      BADstim=0;
   end
   
   stimulus_vals.Inloop.Compiled_FileName=fullfile(EHsignals_dir,filename);
   stimulus_vals.Inloop.Noise_FileName=fullfile(EHsignals_dir,'baseNOISE_LTASS.wav');
   if ~BADstim   
      % If file does not exist, synthesize the new BASELINE vowel and save as a new wavfile
      if isempty(dir(fullfile(EHsignals_dir,filename)))
         % stim returned is 1 cycle of the vowel (Mode 2: synthesizes vowel, but does not show it)
         [stim,Fs,filename,dBreTONE,BASELINE_FormFreqs_Hz]= ...
            synth_BASELINE_eh(BASELINE_TargetFreq_Hz,BASELINE_F0_Hz,BASELINE_Feature,Fix2Harms,2);
         wavwrite(stim,Fs,fullfile(EHsignals_dir,filename))
      end
      
      %% From here, we just read the wavfile and proceed like ORIGINAL EH
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (exist(stimulus_vals.Inloop.Compiled_FileName,'file') ~= 0)
         Llist = {stimulus_vals.Inloop.Compiled_FileName};
         Rlist = {stimulus_vals.Inloop.Noise_FileName};
         [vowel BASELINE_Fs] = wavread(stimulus_vals.Inloop.Compiled_FileName);
         [noise BASENOISE_Fs] = wavread(stimulus_vals.Inloop.Noise_FileName);
         if BASELINE_Fs~=BASENOISE_Fs
            nelerror('In EHrlv_IN_template: BASE sampling rates dont match between EH and NOISE!!');
         end
      else
         Llist = {};
         Rlist = {};
      end
      vowel=vowel*PolarityFact;  % Invert if necessary      
      noise=noise*PolarityFact;  % Invert if necessary
      [stimulus_vals units] = NI_check_gating_params(stimulus_vals, units);
   end
   
   % 12Apr2005: M. Heinz: not sure if this is necessary, but we'll keep the mixing recompute
   %
   %% In this template we have to change the Gating and Mix parameters according to 
   %% the Llist, Rlist and the playback duration of the wav files.
   %% We do this by first, updating the relevant template definitions, and second,
   %% by calling 'structdlg' in its non-interactive invisible mode, to recalculated 'stimulus_vals' fields.
   if (isempty(Llist))
      tmplt.IO_def.Mix = rmfield(tmplt.IO_def.Mix,'Llist');
   end
   if (isempty(Rlist))
      tmplt.IO_def.Mix = rmfield(tmplt.IO_def.Mix,'Rlist');
   end
   [stimulus_vals.Mix units.Mix] = structdlg(tmplt.IO_def.Mix,'',stimulus_vals.Mix,'off');
   
   %% Set maximum sampling rate based on how many wavfiles to be played
   if ~isempty(Llist)&~isempty(Rlist)  % Two files used on NI board, max sampling rate is HALF!
      MAX_NI_SAMPLING_RATE=NI6052UsableRate_Hz(Inf)/2;
   else  % Only 1 file, can use full rate
      MAX_NI_SAMPLING_RATE=NI6052UsableRate_Hz(Inf);
   end   
   
   featureNames = {'T0','F1','T1','F2','T2','F3','T3'};
   featureIND=find(strcmp(featureNames,stimulus_vals.Inloop.Feature));

   if ~BADstim
      % Call here to get params, later call to view vowel is needed
      [BASELINE_FeatFreqs_Hz,BASELINE_FeatLevs_dB,dBreTONE]= ...
         getVowelParams(vowel,BASELINE_Fs,stimulus_vals.Gating.Duration/1000,BASELINE_FormFreqs_Hz,0);
      dBreTONE_noise=20*log10(sqrt(mean(noise.^2))/.707);
   else
      BASELINE_FeatFreqs_Hz=NaN*ones(size(featureNames));
      BASELINE_FeatLevs_dB=NaN*ones(size(featureNames));
      dBreTONE=NaN;
      dBreTONE_noise=NaN;
      BASELINE_Fs=33000;
      Llist={};
      Rlist={};
   end

   % Calculate Target Frequency from Condition Params
   Direction=1;   % Hardwired
   Offset_Direction='above';
   FreqOffset=0;  % Hardwired
   stimulus_vals.Inloop.Computed_FeatureTarget_Hz=stimulus_vals.Inloop.BaseFrequency*1000*2^(Direction*FreqOffset);
   
   if ~BADstim
      stimulus_vals.Inloop.Computed_UpdateRate_Hz= ...
         BASELINE_Fs*(stimulus_vals.Inloop.Computed_FeatureTarget_Hz/BASELINE_FeatFreqs_Hz(featureIND)); %Shift to feature
   else
      stimulus_vals.Inloop.Computed_UpdateRate_Hz=BASELINE_Fs;
   end
   
   if (stimulus_vals.Inloop.Computed_UpdateRate_Hz> MAX_NI_SAMPLING_RATE)
      stimulus_vals.Inloop.Used_UpdateRate_Hz=MAX_NI_SAMPLING_RATE;
      ding
      nelerror(sprintf('In EH_template: Requested sampling rate (%.0f Hz) greater than MAX rate allowed by NI board (%.0f Hz)!', ...
         stimulus_vals.Inloop.Computed_UpdateRate_Hz,MAX_NI_SAMPLING_RATE));
   end

   stimulus_vals.Inloop.Used_UpdateRate_Hz = NI6052UsableRate_Hz(stimulus_vals.Inloop.Computed_UpdateRate_Hz); % GE/MH 04Nov2003:
   % Template forces use of a rate that is valid for the NI6052e board in the
   %  mode in which it is called (see 'd2a.c').
   stimulus_vals.Inloop.Used_FeatureTarget_Hz=stimulus_vals.Inloop.Computed_FeatureTarget_Hz* ...
      stimulus_vals.Inloop.Used_UpdateRate_Hz/stimulus_vals.Inloop.Computed_UpdateRate_Hz;
   if ~BADstim
      if strcmp(stimulus_vals.Inloop.ViewVowel,'baseline')
         [BASELINE_FeatFreqs_Hz,BASELINE_FeatLevs_dB,dBreTONE]= ...
            getVowelParams(vowel,BASELINE_Fs,stimulus_vals.Gating.Duration/1000,BASELINE_FormFreqs_Hz,1);
      elseif strcmp(stimulus_vals.Inloop.ViewVowel,'shifted')
         [SHIFTED_FeatFreqs_Hz,SHIFTED_FeatLevs_dB,SHIFTED_dBreTONE]= ...
            getVowelParams(vowel,stimulus_vals.Inloop.Used_UpdateRate_Hz, ...
            stimulus_vals.Gating.Duration/1000, ...
            BASELINE_FormFreqs_Hz*stimulus_vals.Inloop.Used_UpdateRate_Hz/BASELINE_Fs,2);
      end
   end
   stimulus_vals.Inloop.ViewVowel='no';  % Reset each time
   
   Inloop.Name                         = 'DALinloop_NI_SCC_wavfiles';
   Inloop.params.list                  = Llist;
   Inloop.params.Rlist                 = Rlist; 

   Inloop.params.Condition.BaseFrequency_kHz           = stimulus_vals.Inloop.BaseFrequency;
   Inloop.params.Condition.FreqOffset_octs             = FreqOffset;
   Inloop.params.Condition.Offset_Direction            = Offset_Direction;
   Inloop.params.Condition.Feature                     = stimulus_vals.Inloop.Feature; 
   Inloop.params.Condition.FormsAtHarmonics            = stimulus_vals.Inloop.FormsAtHarmonics;
   Inloop.params.Condition.InvertPolarity              = stimulus_vals.Inloop.InvertPolarity;
   Inloop.params.Condition.Level_dBSPL                 = stimulus_vals.Inloop.Signal_Level;
   Inloop.params.CalibPicNum                           = stimulus_vals.Inloop.CalibPicNum;
   Inloop.params.Condition.High_Attenuation_noise      = stimulus_vals.Inloop.High_Attenuation_noise;
   Inloop.params.Condition.Low_Attenuation_noise       = stimulus_vals.Inloop.Low_Attenuation_noise;
   Inloop.params.Condition.dBstep_Atten_noise          = stimulus_vals.Inloop.dBstep_Atten_noise;
   
   Inloop.params.BASELINE.F0_Hz              = BASELINE_F0_Hz;
   Inloop.params.BASELINE.Feature            = BASELINE_Feature; 
   Inloop.params.BASELINE.TargetFreq_Hz      = BASELINE_TargetFreq_Hz;
   Inloop.params.BASELINE.FormsAtHarmonics   = stimulus_vals.Inloop.FormsAtHarmonics;
   Inloop.params.BASELINE.FileName           = stimulus_vals.Inloop.Compiled_FileName;
   Inloop.params.BASELINE.FormFreqs_Hz       = BASELINE_FormFreqs_Hz;
   Inloop.params.BASELINE.dBreTONE           = dBreTONE;
   Inloop.params.BASELINE.Fs_Hz              = BASELINE_Fs;
   Inloop.params.BASELINE.FeatFreqs_Hz       = BASELINE_FeatFreqs_Hz;
   Inloop.params.BASELINE.FeatLevs_dB        = BASELINE_FeatLevs_dB;
   Inloop.params.BASELINE.FileName_noise     = stimulus_vals.Inloop.Noise_FileName;
   Inloop.params.BASELINE.Fs_Hz_noise        = BASENOISE_Fs;
   Inloop.params.BASELINE.dBreTONE_noise     = dBreTONE_noise;
   Inloop.params.featureNames                = featureNames;

   Inloop.params.Computed.FeatureTarget_Hz   = stimulus_vals.Inloop.Computed_FeatureTarget_Hz;
   Inloop.params.Computed.UpdateRate_Hz   = stimulus_vals.Inloop.Computed_UpdateRate_Hz;
   
   %%% Account for Calibration to set Level in dB SPL
   if ~isempty(stimulus_vals.Inloop.CalibPicNum)
      if stimulus_vals.Inloop.CalibPicNum==0
         max_dBSPL=-999;
      else
         cdd
         if ~isempty(dir(sprintf('p%04d_calib.m',Inloop.params.CalibPicNum)))
            x=loadpic(stimulus_vals.Inloop.CalibPicNum);
            CalibData=x.CalibData(:,1:2);
            CalibData(:,2)=trifilt(CalibData(:,2)',5)';
            max_dBSPL=CalibInterp(stimulus_vals.Inloop.BaseFrequency,CalibData);
         else
            max_dBSPL=[];
            Inloop.params.CalibPicNum=NaN;
         end
         rdd
      end
   else
      max_dBSPL=[];
   end
   Inloop.params.attens                 = max_dBSPL-stimulus_vals.Inloop.Signal_Level+dBreTONE;
   stimulus_vals.Inloop.Computed_Attenuation_dB          = Inloop.params.attens;
   Inloop.params.Rattens               = stimulus_vals.Inloop.High_Attenuation_noise:-stimulus_vals.Inloop.dBstep_Atten_noise:stimulus_vals.Inloop.Low_Attenuation_noise;  % For noise
   % Need to make this the same length for DALinloop
   if length(Inloop.params.attens)==1
      Inloop.params.attens                = Inloop.params.attens*ones(size(Inloop.params.Rattens));
   end
   
   Inloop.params.repetitions           = stimulus_vals.Inloop.Repetitions;
                                       
   Inloop.params.origstimUpdateRate_Hz = BASELINE_Fs;
   Inloop.params.feature               = featureNames{featureIND};
   Inloop.params.featureFreq_Hz        = BASELINE_FeatFreqs_Hz(featureIND);
   Inloop.params.featureTargetFreq_Hz  = stimulus_vals.Inloop.Used_FeatureTarget_Hz;
   Inloop.params.updateRate_Hz         = stimulus_vals.Inloop.Used_UpdateRate_Hz;
     
   DAL.funcName = 'data_acquisition_loop_NI'; % added by GE 30oct2003.
   DAL.Inloop = Inloop;
   DAL.Gating = stimulus_vals.Gating;
   DAL.Mix         = mix_params2devs(stimulus_vals.Mix,used_devices); % GE debug: see 'used_devices.File' line at beginning of function
   DAL.short_description   = 'EHINrlv';
   
%   DAL.endLinePlotParams                  = nel_plot_pst_params(DAL.Gating.Period/1000, DAL.Gating.Duration/1000);  % GE 04Nov2003.

   DAL.description = build_description(DAL,stimulus_vals);
   errstr = check_DAL_params(DAL,fieldname);
end

%----------------------------------------------------------------------------------------
function str = build_description(DAL,stimulus_vals)
p = DAL.Inloop.params;
[fpath,file] = fileparts(stimulus_vals.Inloop.Compiled_FileName);
str{1} = sprintf('IN_RLF %dreps:', p.repetitions);
str{1} = sprintf('%s ''%s''', str{1},file);
str{1} = sprintf('%s(%s),', str{1}, stimulus_vals.Mix.Llist(1));
if ~isempty(p.attens)
   str{1} = sprintf('%s %1.1f dB Attn.', str{1}, p.attens(1));
else
   str{1} = sprintf('%s %1.1f dB Attn.', str{1}, -999);
end
if (length(p.Rattens) > 1)
   str{1} = sprintf('%s noise@ %1.f-%1.f dBatt (%1.1f step)', str{1}, p.Rattens(1), p.Rattens(end), diff(p.Rattens(1:2)));
else
   str{1} = sprintf('%s noise@ %1.1f dB Attn.', str{1}, p.Rattens);
end
str{1} = sprintf('%s(%s),', str{1}, stimulus_vals.Mix.Rlist(1));
str{1} = sprintf('%s Urate: %.0f Hz,', str{1}, stimulus_vals.Inloop.Used_UpdateRate_Hz);
str{1} = sprintf('%s (%s @ %.f Hz)', str{1}, stimulus_vals.Inloop.Feature,stimulus_vals.Inloop.Used_FeatureTarget_Hz);
if strcmp(stimulus_vals.Inloop.InvertPolarity,'yes')
   str{1} = sprintf('%s (Pol:-)', str{1});
else
   str{1} = sprintf('%s (Pol:+)', str{1});
end

%----------------------------------------------------------------------------------------
function errstr = check_DAL_params(DAL,fieldname)
% Some extra error checks
errstr = '';
if (isequal(fieldname,'Inloop'))
   if isnan(DAL.Inloop.params.BASELINE.dBreTONE)  % Indicates BADstim!!!
      errstr=sprintf('Problem with stimulus design!!  BASELINE formants were: %.1f, %.1f, %.1f, %.1f, %.1f Hz ',DAL.Inloop.params.BASELINE.FormFreqs_Hz);
   else
      if (isempty(DAL.Inloop.params.attens))
         errstr = 'Attenuations are not set correctly! (high vs. low mismatch?)';
      elseif (DAL.Inloop.params.attens<0)
         errstr = sprintf('Negative Attenuation of %.1f dB cannot be set, Level must be lowered',DAL.Inloop.params.attens);
      elseif (DAL.Inloop.params.attens>120)
         errstr = sprintf('Attenuation of %.1f dB cannot be set, Level must be raised',DAL.Inloop.params.attens);
      elseif (isempty(DAL.Inloop.params.Rattens))
         errstr = 'Noise Attenuations are not set correctly! (high vs. low mismatch?)';
      elseif (DAL.Inloop.params.Rattens<0)
         errstr = sprintf('Negative Noise Attenuation of %.1f dB cannot be set, Level must be lowered',DAL.Inloop.params.Rattens);
      elseif (DAL.Inloop.params.Rattens>120)
         errstr = sprintf('Noise Attenuation of %.1f dB cannot be set, Level must be raised',DAL.Inloop.params.Rattens);
      elseif isnan(DAL.Inloop.params.Computed.UpdateRate_Hz)   
         errstr = 'Computed UpdateRate is NaN, this feature is undefined!!';
      end
   end
end

%----------------------------------------------------------------------------------------
function tmplt = template_definition(fieldname)
global signals_dir
persistent prev_unit_bf
%%%%%%%%%%%%%%%%%%%%
%% Inloop Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Inloop.BaseFrequency     =  {'current_unit_bf'   'kHz'      [0.04  50] 0 0};
IO_def.Inloop.CalibPicNum  =  {[]   ''       [0 6000]};
IO_def.Inloop.Feature           =  {'T0|{F1}|T1|F2|T2|F3|T3'};
IO_def.Inloop.Signal_Level             =  {65 'dB SPL'       [-50    150]   0  0}; 
IO_def.Inloop.FormsAtHarmonics  =  {'no|{yes}'};
IO_def.Inloop.InvertPolarity  =  {'{no}|yes'};

IO_def.Inloop.Low_Attenuation_noise   = { 1             'dB'    [0    120]      };  % Here, these are for the noise
IO_def.Inloop.High_Attenuation_noise  = { 80           'dB'    [0    120]      };
IO_def.Inloop.dBstep_Atten_noise      = { 2           'dB'    [1    120]      };

IO_def.Inloop.Repetitions  =  {1   ''       [1 600]}; 
IO_def.Inloop.ViewVowel  =  {'{no}|baseline|shifted'};

if (~isequal(current_unit_bf, prev_unit_bf) & isequal(fieldname,'Inloop'))
   IO_def.Inloop.BaseFrequency{5}            = 1; % ignore dflt. Always recalculate.
   prev_unit_bf = current_unit_bf;
end


%%%%%%%%%%%%%%%%%%%%
%% Gating Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Gating.Duration             = {400       'ms'    [20 2000]};
IO_def.Gating.Period               = {'default_period(this.Duration)'    'ms'   [50 5000]};
IO_def.Gating.Rise_fall_time       = {'default_rise_time(this.Duration)' 'ms'   [0  1000]};

%%%%%%%%%%%%%%%%%%%%
%% Mix Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Mix.Llist        =  {'Left|Both|{Right}'};
IO_def.Mix.Rlist        =  {'Left|Both|{Right}'};

tmplt.tag               = 'EH_LTASSrlv_tmplt';
tmplt.IO_def = IO_def;
