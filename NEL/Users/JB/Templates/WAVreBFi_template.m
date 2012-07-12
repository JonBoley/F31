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
      if ~isempty(findstr(Llist{i},'LTASS'))
         FileAttenList(i)=stimulus_vals.Inloop.NOISE_thr-stimulus_vals.Inloop.Thr_Offset;
      elseif ~isempty(findstr(Llist{i},'EH'))
         FileAttenList(i)=stimulus_vals.Inloop.SIGNAL_thr-stimulus_vals.Inloop.Thr_Offset;
      else
         FileAttenList(i)=NaN;
         error('File Mismatch in assigning level (not "LTASS" or "EH"')
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
   Inloop.params.Condition.NOISE_thr               = stimulus_vals.Inloop.NOISE_thr;
   Inloop.params.Condition.SIGNAL_thr               = stimulus_vals.Inloop.SIGNAL_thr;
   Inloop.params.Condition.Thr_Offset               = stimulus_vals.Inloop.Thr_Offset;
   
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
IO_def.Inloop.List_File             = { {['uigetfile(''' signals_dir 'Lists\MH\ISH2009\ISH2009.m'')']} };
IO_def.Inloop.NOISE_thr           = {'max(0,current_unit_thresh-40)'  'dB'    [0    120]      };
IO_def.Inloop.SIGNAL_thr           = {'max(0,current_unit_thresh-50)'  'dB'    [0    120]      };
IO_def.Inloop.Thr_Offset           = {10  'dB'    [0    120]      };
IO_def.Inloop.Repetitions           = { 25                        ''      [1    Inf]      };
IO_def.Inloop.BaseUpdateRate        = { 33000                  'Hz'      [1    NI6052UsableRate_Hz(Inf)]      };


%%%%%%%%%%%%%%%%%%%%
%% Gating Section 
%%%%%%%%%%%%%%%%%%%%
% NOTE: these are for the BASELINE stim (WAV files should be desired duration, with equal zeros at end of WAV file to allow STMP shifting!)
IO_def.Gating.Duration             = {1700       'ms'    [20 4000]};
IO_def.Gating.Period               = {2200    'ms'   [50 5000]};
IO_def.Gating.Rise_fall_time       = {'default_rise_time(this.Duration)' 'ms'   [0  1000]}; 
%%%%%%%%%%%%%%%%%%%%
%% Mix Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Mix.Llist        =  {'Left|Both|{Right}'};
tmplt.tag               = 'WAVreBFi_tmplt';
tmplt.IO_def = IO_def;
