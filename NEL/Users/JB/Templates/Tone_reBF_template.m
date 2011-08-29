function [tmplt,DAL,stimulus_vals,units,errstr] = Tone_reBF_template(fieldname,stimulus_vals,units)
%
% Modified by M.Heinz 12Nov2003, from nel_rate_level_template.m

% Adapted from "nel_pst_template.m" by GE, 29Mar2002.

used_devices.Tone       = 'RP1.1';
tmplt = template_definition(fieldname);
if (exist('stimulus_vals','var') == 1)
   Inloop.Name                               = 'DALinloop_general_TN';
   Inloop.params.repetitions                 = stimulus_vals.Inloop.Repetitions;
   Inloop.params.Condition.BaseFrequency_kHz           = stimulus_vals.Inloop.BaseFrequency;
   Inloop.params.Condition.FreqOffset_octs             = stimulus_vals.Inloop.FreqOffset;
   Inloop.params.Condition.Offset_Direction            = stimulus_vals.Inloop.Offset_Direction;
   Inloop.params.Condition.Level_dBSPL                 = stimulus_vals.Inloop.Level;
   Inloop.params.CalibPicNum                 = stimulus_vals.Inloop.CalibPicNum;

   Inloop.params.main.source                 = 'Tone';
   clear Direction
   if strcmp(stimulus_vals.Inloop.Offset_Direction,'above')
      Direction=1;
   elseif strcmp(stimulus_vals.Inloop.Offset_Direction,'below')
      Direction=-1;
   else
      error('stimulus_vals.Inloop.Offset_Direction NOT SET CORRECTLY!');
   end
   if ischar(stimulus_vals.Inloop.FreqOffset)
      FreqOffset=str2num(stimulus_vals.Inloop.FreqOffset);
   else
      FreqOffset=stimulus_vals.Inloop.FreqOffset;
   end
   Inloop.params.main.tone.freq                 = stimulus_vals.Inloop.BaseFrequency*1000*2^(Direction*FreqOffset);
   stimulus_vals.Inloop.calculated_frequency_Hz = Inloop.params.main.tone.freq;
   
   Inloop.params.main.tone.bw                = 0;
   Inloop.params.main.noise.low_cutoff       = 0;
   Inloop.params.main.noise.high_cutoff      = 0;

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
            max_dBSPL=CalibInterp(Inloop.params.main.tone.freq/1000,CalibData);
         else
            max_dBSPL=[];
            Inloop.params.CalibPicNum=NaN;
         end
         rdd
      end
   else
      max_dBSPL=[];
   end
   Inloop.params.main.attens                 = max_dBSPL-stimulus_vals.Inloop.Level;
   stimulus_vals.Inloop.calculated_Attenuation_dB   = Inloop.params.main.attens;
%    disp(sprintf('Stim Freq: %.2f kHz; Atten: %.1f dB',Inloop.params.main.tone.freq,Inloop.params.main.attens))
   
   Inloop.params.secondary.source            = 'None';
   Inloop.params.secondary.tone.freq         = 0;
   Inloop.params.secondary.noise.low_cutoff  = 0;
   Inloop.params.secondary.noise.high_cutoff = 0;
   Inloop.params.secondary.noise.gating      = '';
   Inloop.params.secondary.noise.adaptation  = 0;
   Inloop.params.secondary.atten             = [];
   Inloop.params.rise_fall                   = stimulus_vals.Gating.Rise_fall_time;

   DAL.Inloop = Inloop;
   DAL.Gating = stimulus_vals.Gating;
   DAL.short_description   = 'TrBF';
%    DAL.endLinePlotParams                  = nel_plot_pst_params(DAL.Gating.Period/1000, DAL.Gating.Duration/1000);  % GE 04Nov2003.

   % [stimulus_vals.Mix units.Mix] = structdlg(tmplt.IO_def.Mix,'',stimulus_vals.Mix,'off');
   DAL.Mix         = mix_params2devs(stimulus_vals.Mix,used_devices);

   DAL.description = build_description(DAL,stimulus_vals);
   errstr = check_DAL_params(DAL,fieldname);
end

%----------------------------------------------------------------------------------------
function str = build_description(DAL,stimulus_vals)
p = DAL.Inloop.params;
str{1} = sprintf('PST %d Repetitions', p.repetitions);
str{1} = sprintf('%s %1.2f kHz Tone', str{1}, p.main.tone.freq/1000);
str{1} = sprintf('%s @ %1.1f dB Attn.', str{1}, p.main.attens);
str{1} = sprintf('%s (%s)', str{1}, stimulus_vals.Mix.Tone);

%----------------------------------------------------------------------------------------
function errstr = check_DAL_params(DAL,fieldname)
% Some extra error checks
errstr = '';
if (isequal(fieldname,'Inloop'))
   if isempty(DAL.Inloop.params.CalibPicNum)
      errstr = 'Need to set Calibration PicNum!';
   elseif isnan(DAL.Inloop.params.CalibPicNum)
      errstr = 'Not a valid Calibration PicNum! (use 0 to escape)';      
   elseif DAL.Inloop.params.CalibPicNum==0
      errstr = '';      % Allows you to escape the structdlg to go find a valid calibration file
   else
      if (isempty(DAL.Inloop.params.main.attens))
         errstr = 'Attenuations are not set correctly! (high vs. low mismatch?)';
      elseif (DAL.Inloop.params.main.attens<0)
         errstr = sprintf('Negative Attenuation of %.1f dB cannot be set, Level must be lowered',DAL.Inloop.params.main.attens);
      elseif (DAL.Inloop.params.main.attens>120)
         errstr = sprintf('Attenuation of %.1f dB cannot be set, Level must be raised',DAL.Inloop.params.main.attens);
      end
      if (isempty(DAL.Inloop.params.main.tone.freq))
         errstr = 'Tone Frequency is empty!';
      end
   end
end

%----------------------------------------------------------------------------------------
function tmplt = template_definition(fieldname)
persistent prev_unit_bf prev_unit_thresh
%%%%%%%%%%%%%%%%%%%%
%% Inloop Section 
%%%%%%%%%%%%%%%%%%%%

IO_def.Inloop.BaseFrequency     =  {'current_unit_bf'   'kHz'      [0.04  50]   0  0};
IO_def.Inloop.FreqOffset        =  {'{0}|1/4|1/2|3/4|1' 'octaves'};
IO_def.Inloop.Offset_Direction  =  {'{below}|above' 're: BF'};
IO_def.Inloop.CalibPicNum  =  {[]   ''       [0 6000]};
IO_def.Inloop.Level  =  {65 'dB SPL'       [-50    150]   0  0}; 
IO_def.Inloop.Repetitions  =  {100   ''       [1 600]}; 

if (isequal(fieldname,'Inloop'))
   if (~isequal(current_unit_bf, prev_unit_bf))
      IO_def.Inloop.BaseFrequency{5}            = 1; % ignore dflt. Always recalculate.
      prev_unit_bf = current_unit_bf;
   end
end

%%%%%%%%%%%%%%%%%%%%
%% Gating Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Gating.Duration         = {400       'ms'    [20 60000]};
IO_def.Gating.Period           = {'default_period(this.Duration)'    'ms'   [50 120000]};
IO_def.Gating.Rise_fall_time   = {'default_rise_time(this.Duration)' 'ms'   [0  1000]}; 

%%%%%%%%%%%%%%%%%%%%
%% Mix Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Mix.Tone      =  {'Left|Both|{Right}'};

tmplt.tag         = 'TrBF_tmplt';
tmplt.IO_def = IO_def;
