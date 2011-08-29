function [tmplt,DAL,stimulus_vals,units,errstr] = BBN_reBFi_template(fieldname,stimulus_vals,units)
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
    OctShifts=[-0.25,-0.15,-0.05,0,0.05,0.15,0.30];
    stimulus_vals.Inloop.OctShifts=OctShifts;
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
   
   %featureIND_List=zeros(1,length(Features)*length(OctShifts));
%    Computed_FeatureTarget_Hz_List=zeros(1,length(Features)*length(OctShifts));
   Computed_FeatureTarget_Hz_List=zeros(1,length(Llist)*length(OctShifts));
%    Computed_UpdateRate_Hz_List=zeros(1,length(Features)*length(OctShifts));
   Computed_UpdateRate_Hz_List=zeros(1,length(Llist)*length(OctShifts));
%    Used_UpdateRate_Hz_List=zeros(1,length(Features)*length(OctShifts));
   Used_UpdateRate_Hz_List=zeros(1,length(Llist)*length(OctShifts));
%    Used_FeatureTarget_Hz_List=zeros(1,length(Features)*length(OctShifts));
Used_FeatureTarget_Hz_List=zeros(1,length(Llist)*length(OctShifts));
%    Features_List=cell(1,length(Features)*length(OctShifts));
Features_List=cell(1,length(Llist)*length(OctShifts));
%    OctShifts_List=zeros(1,length(Features)*length(OctShifts));
  OctShifts_List=zeros(1,length(Llist)*length(OctShifts));
   
   %%%%%%%% Generate Lists here
   condIND=0;
   for FeatIND=1:length(Llist)
      for shiftIND=1:length(OctShifts)
         condIND=condIND+1;
         Features_List{condIND}=Llist{FeatIND};
         OctShifts_List(condIND)=OctShifts(shiftIND);
         
         %          stimulus_vals.Inloop.Feature='F1';
         %featureIND_List(condIND)=find(strcmp(featureNames,Features(FeatIND)));
%          Computed_FeatureTarget_Hz_List(condIND)=stimulus_vals.Inloop.BaseFrequency*1000*2^(OctShifts(shiftIND));
         Computed_FeatureTarget_Hz_List(condIND)=stimulus_vals.Inloop.BestFrequency *1000*2^(OctShifts(shiftIND));
         Computed_UpdateRate_Hz_List(condIND)= ...
           stimulus_vals.Inloop.BaseUpdateRate*(Computed_FeatureTarget_Hz_List(condIND)/(stimulus_vals.Inloop.BestFrequency*1000)); 
       %since it's a frozen noise and not F1 or F2 it does not have specific feature frequency. Rather noise is presented at BF 
       %originally hence it is being taken and shifted to new frequency by changing the sampling rate

         
%          if ~BADstim
%             Computed_UpdateRate_Hz_List(condIND)= ...
%                BASELINE_Fs*(Computed_FeatureTarget_Hz_List(condIND)/BASELINE_FeatFreqs_Hz(featureIND_List(condIND))); %Shift to feature
%          else
%             Computed_UpdateRate_Hz_List(condIND)=BASELINE_Fs;
%          end
         
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
%    Inloop.params.Condition.Level_dBSPL                 = stimulus_vals.Inloop.Level;
%    Inloop.params.CalibPicNum                 = stimulus_vals.Inloop.CalibPicNum;
   
%    Inloop.params.BASELINE.F0_Hz              = BASELINE_F0_Hz;
%    Inloop.params.BASELINE.Feature            = BASELINE_Feature; 
%    Inloop.params.BASELINE.TargetFreq_Hz      = BASELINE_TargetFreq_Hz;
%    Inloop.params.BASELINE.FormsAtHarmonics   = stimulus_vals.Inloop.FormsAtHarmonics;
%    Inloop.params.BASELINE.FileName           = stimulus_vals.Inloop.Compiled_FileName;
%    Inloop.params.BASELINE.FormFreqs_Hz       = BASELINE_FormFreqs_Hz;
%    Inloop.params.BASELINE.dBreTONE           = dBreTONE;
%    Inloop.params.BASELINE.Fs_Hz              = stimulus_vals.Inloop.BaseFs;
%    Inloop.params.BASELINE.FeatFreqs_Hz       = BASELINE_FeatFreqs_Hz;
%    Inloop.params.BASELINE.FeatLevs_dB        = BASELINE_FeatLevs_dB;
%    Inloop.params.featureNames                = featureNames;

   Inloop.params.Computed.FeatureTarget_Hz_List   = stimulus_vals.Inloop.Computed_FeatureTarget_Hz_List;
   Inloop.params.Computed.UpdateRate_Hz_List   = stimulus_vals.Inloop.Computed_UpdateRate_Hz_List;
   Inloop.params.Used.FeatureTarget_Hz_List   = stimulus_vals.Inloop.Used_FeatureTarget_Hz_List;
   Inloop.params.Used.UpdateRate_Hz_List       = stimulus_vals.Inloop.Used_UpdateRate_Hz_List;
   
   %%% Account for Calibration to set Level in dB SPL
%    if ~isempty(stimulus_vals.Inloop.CalibPicNum)
%       if stimulus_vals.Inloop.CalibPicNum==0
%          max_dBSPL=-999;
%       else
%          cdd
%          if ~isempty(dir(sprintf('p%04d_calib.m',Inloop.params.CalibPicNum)))
%             x=loadpic(stimulus_vals.Inloop.CalibPicNum);
%             CalibData=x.CalibData(:,1:2);
%             CalibData(:,2)=trifilt(CalibData(:,2)',5)';
%             max_dBSPL=CalibInterp(stimulus_vals.Inloop.BaseFrequency,CalibData);
%          else
%             max_dBSPL=[];
%             Inloop.params.CalibPicNum=NaN;
%          end
%          rdd
%       end
%    else
%       max_dBSPL=[];
%    end
   Inloop.params.attens                = stimulus_vals.Inloop.Attenuation;
   Inloop.params.Rlist                 = []; % 
   Inloop.params.Rattens               = []; %      
   Inloop.params.repetitions           = stimulus_vals.Inloop.Repetitions;
   
   %    Inloop.params.origstimUpdateRate_Hz = BASELINE_Fs;
   Inloop.params.Used.Features_List           = Features_List;
   Inloop.params.Used.OctShifts_List           = OctShifts_List;
   %    Inloop.params.featureFreq_Hz_List        = BASELINE_FeatFreqs_Hz(featureIND);
   %    Inloop.params.featureTargetFreq_Hz_List  = stimulus_vals.Inloop.Used_FeatureTarget_Hz_List;
   Inloop.params.updateRate_Hz         = stimulus_vals.Inloop.Used_UpdateRate_Hz_List;
     
   DAL.funcName = 'data_acquisition_loop_NI'; % added by GE 30oct2003.
   DAL.Inloop = Inloop;
   DAL.Gating = stimulus_vals.Gating;
   DAL.Mix         = mix_params2devs(stimulus_vals.Mix,used_devices); % GE debug: see 'used_devices.File' line at beginning of function
   DAL.short_description   = 'BBNrBFi';
   
%    DAL.endLinePlotParams                  = nel_plot_pst_params(DAL.Gating.Period/1000, DAL.Gating.Duration/1000);  % GE 04Nov2003.

   DAL.description = build_description(DAL,stimulus_vals);
   errstr = check_DAL_params(DAL,fieldname);
end

%----------------------------------------------------------------------------------------
function str = build_description(DAL,stimulus_vals)
p = DAL.Inloop.params;
[listpath,listfile] = fileparts(stimulus_vals.Inloop.List_File);
str{1} = sprintf('%d reps:', p.repetitions);
str{1} = sprintf('%s List ''%s'' (%d files; %d conds) ', str{1},listfile, length(unique(p.list)), length(p.list));
str{1} = sprintf('%s @%1.1f dBatt,', str{1}, p.attens);
   
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
% IO_def.Inloop.BaseF0            =  {75   'Hz'      [0.01  1000] 0 0};
IO_def.Inloop.BestFrequency     =  {'current_unit_bf'   'kHz'      [0.02  20] 0 0};
% IO_def.Inloop.CalibPicNum  =  {[]   ''       [0 6000]};
% IO_def.Inloop.OctShifts_low        =  {-0.25 'octaves'};
% IO_def.Inloop.OctShifts_high        =  {.25 'octaves'};
% IO_def.Inloop.OctShifts_step        =  {.05 'octaves'};
%IO_def.Inloop.Offset_Direction  =  {'{below}|above' 're: BF'};

% IO_def.Inloop.List_File             = { {['uigetfile(''' signals_dir 'Lists\SK\BBN_SCCXCC\*.m'')']} };
IO_def.Inloop.List_File             = { {['uigetfile(''' signals_dir 'Lists\SK\BBN_SCCXCC\BBN_SCCXCC_A2.m'')']} };
IO_def.Inloop.Attenuation           = {'max(0,current_unit_thresh-50)'  'dB'    [0    120]      };
IO_def.Inloop.Repetitions           = { 25                        ''      [1    Inf]      };
IO_def.Inloop.BaseUpdateRate        = { 33000                  'Hz'      [1    NI6052UsableRate_Hz(Inf)]      };

if (~isequal(current_unit_bf, prev_unit_bf) & isequal(fieldname,'Inloop'))
   IO_def.Inloop.BestFrequency{5}            = 1; % ignore dflt. Always recalculate.
   prev_unit_bf = current_unit_bf;
end


%%%%%%%%%%%%%%%%%%%%
%% Gating Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Gating.Duration             = {1700       'ms'    [20 2000]};
IO_def.Gating.Period               = {2200    'ms'   [50 5000]};
IO_def.Gating.Rise_fall_time       = {'default_rise_time(this.Duration)' 'ms'   [0  1000]}; 
%%%%%%%%%%%%%%%%%%%%
%% Mix Section 
%%%%%%%%%%%%%%%%%%%%
% IO_def.Mix.File        =  {'Left|Both|{Right}'};
IO_def.Mix.Llist        =  {'Left|Both|{Right}'};
tmplt.tag               = 'BBNrBFi_tmplt';
tmplt.IO_def = IO_def;
