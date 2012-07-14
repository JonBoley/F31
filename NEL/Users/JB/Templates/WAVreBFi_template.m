function [tmplt,DAL,stimulus_vals,units,errstr] = WAV_reBFi_template(fieldname,stimulus_vals,units)
% MH 15Dec2008 - for ISH2009
% copied from BBNrBFi (SK), made more general for any WAV file list, and multiple levels (Still needs +/- as separate stim files)
%
% MH 10Nov2004: for NOHR project
% Modified version to allow interleaving of conditions (
% From EH_reBF_template
%
% Initial creation LIMITATIONS:
%     1) limited to one Polarity
%     2) assumes 1 filename
%     3) probably limited somewhere to 1 level, if more than 1 ur??  THIS SHOULD BE FIXED EVENTUALLY
%
%
%    places features near BF for a BASELINE EH with F2 at BF and F0=75 Hz
%
% From EH_template (CNexps)
% Template for EH-vowel RLFs, using NI board to allow resampling
%
% MH 07Nov2003, modified from 'nel_NI_wavefile_template'
%
%  Written by GE, adapted from 'nel_wavefile_template' written by AF (11/26/01).
%  Modification dates: 03oct2003, 04oct2003.

% persistent   prev_playdur  prev_min_period  % GE debug: not really necessary???
% We use the persistent variables to detect a change that requires some fields update.
% For example, if the play duration is changed we would like to update the gating information.
% We restict the automatic updates to allow the user the overide them.
persistent prev_maxlen

% used_devices.File         = 'L3';    % ge debug: what if R3 is used also?
used_devices.Llist         = 'L3';   % added by GE 26Jul2002
tmplt = template_definition(fieldname);
if (exist('stimulus_vals','var') == 1)
   if isempty(stimulus_vals.Inloop.BaseUpdateRate)
      stimulus_vals.Inloop.BaseUpdateRate=33000;  % Will get through, and let you change later
   end
   
   %PolarityFact=(strcmp(stimulus_vals.Inloop.InvertPolarity,'yes')-.5)*-2;  % For inverting waveform if necessary
   %     OctShifts=[-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.30,0.40,0.50];  % 11 CFs: NOHR
%    OctShifts=[-0.25,-0.15,-0.05,0,0.05,0.15,0.30];  % 7 CFs: R01 (122208)
%    OctShifts=[-0.5 -0.35,-0.15,-0.05,0,0.05,0.15,0.40 0.50];  % 9 CFs: R01 (011909-add 2 and shift few)
   OctShifts=[-.75 -.50 -.25 -.15 -.05 0 .05 .15 .25 .35 .50]; %JDB (added July 12 2012)
   stimulus_vals.Inloop.OctShifts=OctShifts;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% FIX Duration to allow for proper STMP shifting of WAV files
   BASELINE_Duration = 2000; % (based on specific WAV files used)
   EXTENDED_Duration = BASELINE_Duration * 2^(-min(stimulus_vals.Inloop.OctShifts));
   OFFtime = stimulus_vals.Gating.Period - stimulus_vals.Gating.Duration;
   stimulus_vals.Gating.Duration = EXTENDED_Duration;
   stimulus_vals.Gating.Period = EXTENDED_Duration + OFFtime;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % GENERATE STIMULI
   
   % Hearing Aid Prescriptions
   % strategy = 0 (no gain)
   %        or {1 or 'linear'}
   %        or {2 or 'nonlinear_quiet'}
   %        or {3 or 'nonlinear_noise'}
   strategy=0;
   if strcmp(stimulus_vals.Inloop.HearingAid_Linear,'yes'), strategy = 1;  end
   if strcmp(stimulus_vals.Inloop.HearingAid_Nonlinear,'yes'), strategy = 2;  end
   % only use 3 in noise (never specify it here)
   
   % Write WAV files for F1,F2, and 3 noise levels (this may take some time...)
   newMixedAtten = writeVowelSTMPwavs(stimulus_vals.Inloop.CalibPicNum,...
       stimulus_vals.Inloop.BaseFrequency, stimulus_vals.Inloop.Signal_Level,...
       stimulus_vals.Inloop.Noise_Atten_mid,strategy);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   if (exist(stimulus_vals.Inloop.List_File,'file') ~= 0)
      [Llist,Rlist] = read_rotate_list_file(stimulus_vals.Inloop.List_File);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [data fs] = wavread(Llist{1});
      [stimulus_vals units] = NI_check_gating_params(stimulus_vals, units);%optional??
      [stimulus_vals.Mix units.Mix] = structdlg(tmplt.IO_def.Mix,'',stimulus_vals.Mix,'off');
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
   else
      Llist = [];
      Rlist = [];
      prev_maxlen = 0;
   end

   %%%%%%%%%%%%%%%%%%%%
   % Setup levels to run (allow different attens for diferent files)
   % LATER : ADD list of levels
   %    AttenOFFsets=[10 30 50];
   %    stimulus_vals.Inloop.NOISE_thr;
   FileAttenList=zeros(size(Llist));
   for i=1:length(Llist)
      if ~isempty(findstr(Llist{i},'EHLTASS'))
         FileAttenList(i)=stimulus_vals.Inloop.NOISE_thr-stimulus_vals.Inloop.Thr_Offset;
      else
         FileAttenList(i)=NaN;
         error('File Mismatch in assigning level (not "EHLTASS")')
      end
   end
   
   
   % TODO 121508:
   % - Atten List is all set, now add to attens, and make sure DALinloop works and saves correctly
   % **- then, set based on thr SP/BBN,
   % LATER:- then add multiple levels
  
   
   Computed_FeatureTarget_Hz_List=zeros(1,length(Llist)*length(OctShifts));
   Computed_UpdateRate_Hz_List=zeros(1,length(Llist)*length(OctShifts));
   Used_UpdateRate_Hz_List=zeros(1,length(Llist)*length(OctShifts));
   Used_FeatureTarget_Hz_List=zeros(1,length(Llist)*length(OctShifts));
   Features_List=cell(1,length(Llist)*length(OctShifts));
   OctShifts_List=zeros(1,length(Llist)*length(OctShifts));
   Attens_List=zeros(1,length(Llist)*length(OctShifts));

   stimulus_vals.Inloop.BestFrequency=current_unit_bf;  % Don't need to specify this in IOdef

   %%%%%%%% Generate Lists here
   condIND=0;
   for FeatIND=1:length(Llist)
      for shiftIND=1:length(OctShifts)
         condIND=condIND+1;
         Features_List{condIND}=Llist{FeatIND};
         Attens_List(condIND)=FileAttenList(FeatIND);
         
         OctShifts_List(condIND)=OctShifts(shiftIND);
         
         Computed_FeatureTarget_Hz_List(condIND)=stimulus_vals.Inloop.BestFrequency*1000*2^(OctShifts(shiftIND));
         Computed_UpdateRate_Hz_List(condIND)=stimulus_vals.Inloop.BaseUpdateRate*2^(OctShifts(shiftIND));        
         %since it's a frozen noise and not F1 or F2 it does not have specific feature frequency. Rather noise is presented at BF 
         %originally hence it is being taken and shifted to new frequency by changing the sampling rate       
         
         if (Computed_UpdateRate_Hz_List(condIND)> NI6052UsableRate_Hz(Inf))
            Used_UpdateRate_Hz_List(condIND)=NI6052UsableRate_Hz(Inf);
            nelerror('In EH_template: Requested sampling rate greater than MAX rate allowed by NI board!!');
         end
         Used_UpdateRate_Hz_List(condIND) = NI6052UsableRate_Hz(Computed_UpdateRate_Hz_List(condIND)); % GE/MH 04Nov2003:
         
         % Template forces use of a rate that is valid for the NI6052e board in the
         %  mode in which it is called (see 'd2a.c').
         Used_FeatureTarget_Hz_List(condIND)=Computed_FeatureTarget_Hz_List(condIND)* ...
            Used_UpdateRate_Hz_List(condIND)/Computed_UpdateRate_Hz_List(condIND);
      end
   end
   Llist=Features_List;
   
   stimulus_vals.Inloop.Computed_FeatureTarget_Hz_List=Computed_FeatureTarget_Hz_List;
   stimulus_vals.Inloop.Computed_UpdateRate_Hz_List=Computed_UpdateRate_Hz_List;
   stimulus_vals.Inloop.Used_UpdateRate_Hz_List=Used_UpdateRate_Hz_List;
   stimulus_vals.Inloop.Used_FeatureTarget_Hz_List=Used_FeatureTarget_Hz_List;
   stimulus_vals.Inloop.File_List=Features_List;
   
   Inloop.Name                         = 'DALinloop_NI_SCCi_wavfiles';
   Inloop.params.list                  = Llist;
   Inloop.params.Condition.BaseUpdateRate           = stimulus_vals.Inloop.BaseUpdateRate;
   Inloop.params.Condition.BestFrequency           = stimulus_vals.Inloop.BestFrequency;
   Inloop.params.Condition.OctShifts               = stimulus_vals.Inloop.OctShifts;
   Inloop.params.Condition.InvertPolarity          = 'no';
   
   Inloop.params.Condition.BaseFrequency        = stimulus_vals.Inloop.BaseFrequency;
   Inloop.params.Condition.CalibPicNum          = stimulus_vals.Inloop.CalibPicNum;
   Inloop.params.Condition.Signal_Level         = stimulus_vals.Inloop.Signal_Level;
   Inloop.params.Condition.Noise_Atten_mid      = stimulus_vals.Inloop.Noise_Atten_mid;
   Inloop.params.Condition.use_F1               = stimulus_vals.Inloop.use_F1;
   Inloop.params.Condition.use_F2               = stimulus_vals.Inloop.use_F2;
   Inloop.params.Condition.FormsAtHarmonics     = stimulus_vals.Inloop.FormsAtHarmonics;
   Inloop.params.Condition.HearingAid_None      = stimulus_vals.Inloop.HearingAid_None;
   Inloop.params.Condition.HearingAid_Linear    = stimulus_vals.Inloop.HearingAid_Linear;
   Inloop.params.Condition.HearingAid_Nonlinear = stimulus_vals.Inloop.HearingAid_Nonlinear;
   
   Inloop.params.Computed.FeatureTarget_Hz_List   = stimulus_vals.Inloop.Computed_FeatureTarget_Hz_List;
   Inloop.params.Computed.UpdateRate_Hz_List   = stimulus_vals.Inloop.Computed_UpdateRate_Hz_List;
   Inloop.params.Used.FeatureTarget_Hz_List   = stimulus_vals.Inloop.Used_FeatureTarget_Hz_List;
   Inloop.params.Used.UpdateRate_Hz_List       = stimulus_vals.Inloop.Used_UpdateRate_Hz_List;
   
   Inloop.params.attens                = Attens_List;
   Inloop.params.Rlist                 = []; % 
   Inloop.params.Rattens               = []; %      
   Inloop.params.repetitions           = stimulus_vals.Inloop.Repetitions;
   
   Inloop.params.Used.Features_List           = Features_List;
   Inloop.params.Used.OctShifts_List           = OctShifts_List;
   Inloop.params.Used.Attens_List           = Attens_List;
   Inloop.params.updateRate_Hz         = stimulus_vals.Inloop.Used_UpdateRate_Hz_List;
   
   DAL.funcName = 'data_acquisition_loop_NI'; % added by GE 30oct2003.
   DAL.Inloop = Inloop;
   DAL.Gating = stimulus_vals.Gating;
   DAL.Mix         = mix_params2devs(stimulus_vals.Mix,used_devices); % GE debug: see 'used_devices.File' line at beginning of function
   DAL.short_description   = 'WAVreBFi';
      
   DAL.description = build_description(DAL,stimulus_vals);
   errstr = check_DAL_params(DAL,fieldname);
end

%----------------------------------------------------------------------------------------
function str = build_description(DAL,stimulus_vals)
p = DAL.Inloop.params;
[listpath,listfile] = fileparts(stimulus_vals.Inloop.List_File);
str{1} = sprintf('%d reps:', p.repetitions);
str{1} = sprintf('%s List ''%s'' (%d files x %d CFs) ', str{1},listfile, length(unique(p.list)), length(p.Condition.OctShifts));
str{1} = sprintf('%s @THR+%1.1fdB,', str{1}, p.Condition.Thr_Offset);

%----------------------------------------------------------------------------------------
function errstr = check_DAL_params(DAL,fieldname)
% Some extra error checks
errstr = '';
% if (isequal(fieldname,'Inloop'))
%    if isnan(DAL.Inloop.params.BASELINE.dBreTONE)  % Indicates BADstim!!!
%       errstr=sprintf('Problem with stimulus design!!  BASELINE formants were: %.1f, %.1f, %.1f, %.1f, %.1f Hz ',DAL.Inloop.params.BASELINE.FormFreqs_Hz);
%    else
%       if isempty(DAL.Inloop.params.CalibPicNum)
%          errstr = 'Need to set Calibration PicNum!';
%       elseif isnan(DAL.Inloop.params.CalibPicNum)
%          errstr = 'Not a valid Calibration PicNum! (use 0 to escape)';      
%       elseif DAL.Inloop.params.CalibPicNum==0
%          errstr = '';      % Allows you to escape the structdlg to go find a valid calibration file
%       else
%          if (isempty(DAL.Inloop.params.attens))
%             errstr = 'Attenuations are not set correctly! (high vs. low mismatch?)';
%          elseif (DAL.Inloop.params.attens<0)
%             errstr = sprintf('Negative Attenuation of %.1f dB cannot be set, Level must be lowered',DAL.Inloop.params.attens);
%          elseif (DAL.Inloop.params.attens>120)
%             errstr = sprintf('Attenuation of %.1f dB cannot be set, Level must be raised',DAL.Inloop.params.attens);
%          elseif sum(isnan(DAL.Inloop.params.Computed.UpdateRate_Hz_List))
%             errstr = 'Computed UpdateRate is NaN, one chosen feature is undefined!!';
%          end
%       end
%    end
% end

%----------------------------------------------------------------------------------------
function tmplt = template_definition(fieldname)
global signals_dir
persistent prev_unit_bf
%%%%%%%%%%%%%%%%%%%%
%% Inloop Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Inloop.BaseFrequency    =  {'current_unit_bf'   'kHz'      [0.04  50] 0 0};
IO_def.Inloop.CalibPicNum  =  {[]   ''       [0 6000]};
IO_def.Inloop.Signal_Level     =  {65 'dB SPL'       [-50    150]   0  0}; 
IO_def.Inloop.Noise_Atten_mid  =  {30 'dB atten'       [0    120]   0  0};
IO_def.Inloop.use_F1           =  {'no|{yes}'};
IO_def.Inloop.use_F2           =  {'no|{yes}'};
IO_def.Inloop.FormsAtHarmonics =  {'{no}|yes'};
IO_def.Inloop.HearingAid_None       =  {'no|{yes}'};
IO_def.Inloop.HearingAid_Linear     =  {'{no}|yes'};
IO_def.Inloop.HearingAid_Nonlinear  =  {'{no}|yes'};

IO_def.Inloop.List_File             = { {['uigetfile(''' signals_dir 'Lists\JB\VowelLTASS\VowelLTASS.m'')']} };
IO_def.Inloop.Repetitions           = { 25                        ''      [1    Inf]      };
IO_def.Inloop.BaseUpdateRate        = { 33000                  'Hz'      [1    NI6052UsableRate_Hz(Inf)]      };


%%%%%%%%%%%%%%%%%%%%
%% Gating Section 
%%%%%%%%%%%%%%%%%%%%
% NOTE: these are for the BASELINE stim (WAV files should be desired duration, with equal zeros at end of WAV file to allow STMP shifting!)
IO_def.Gating.Duration             = {2000       'ms'    [20 4000]};
IO_def.Gating.Period               = {3000    'ms'   [50 5000]};
IO_def.Gating.Rise_fall_time       = {'default_rise_time(this.Duration)' 'ms'   [0  1000]}; 
%%%%%%%%%%%%%%%%%%%%
%% Mix Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Mix.Llist        =  {'Left|Both|{Right}'};
tmplt.tag               = 'WAVreBFi_tmplt';
tmplt.IO_def = IO_def;
