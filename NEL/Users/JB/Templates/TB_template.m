function [tmplt,DAL,stimulus_vals,units,errstr] = TB_template(fieldname,stimulus_vals,units)
%
% Template for Recruitment Stimulus TB: BF tone, rate-level, 200-ms duration
%
% MH 12/17/01

used_devices.Tone       = 'RP1.1';
tmplt = template_definition(fieldname);
if (exist('stimulus_vals','var') == 1)
   Inloop.Name                               = 'DALinloop_general_TN';
   Inloop.params.main.source                 = 'Tone';
   Inloop.params.main.tone.freq              = stimulus_vals.Inloop.Frequency*1000;
   Inloop.params.main.tone.bw                = 0;
   Inloop.params.main.noise.low_cutoff       = 0;
   Inloop.params.main.noise.high_cutoff      = 0;
   Inloop.params.main.attens                 = stimulus_vals.Inloop.High_Attenuation :-1:stimulus_vals.Inloop.Low_Attenuation;
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
   DAL.short_description   = 'TB';

   % [stimulus_vals.Mix units.Mix] = structdlg(tmplt.IO_def.Mix,'',stimulus_vals.Mix,'off');
   DAL.Mix         = mix_params2devs(stimulus_vals.Mix,used_devices);

   DAL.description = build_description(DAL,stimulus_vals);
   errstr = check_DAL_params(DAL,fieldname);
   
   %%%%%%%
   % If parameters are NOT correct for this template, Take away this template name
   if (length(DAL.Inloop.params.main.attens) == 1)
      DAL.short_description   = '';
   end
end

%----------------------------------------------------------------------------------------
function str = build_description(DAL,stimulus_vals)
p = DAL.Inloop.params;
str{1} = sprintf('%1.2f kHz Tone', p.main.tone.freq/1000);
if (length(p.main.attens) > 1)
   str{1} = sprintf('%s @ %1.1f - %1.1f dB Attn.', str{1}, p.main.attens(1), p.main.attens(end));
else
   if (~isempty(p.main.attens))
      str{1} = sprintf('%s @ %1.1f dB Attn.', str{1}, p.main.attens(1));
   end
end
str{1} = sprintf('%s (%s)', str{1}, stimulus_vals.Mix.Tone);

%----------------------------------------------------------------------------------------
function errstr = check_DAL_params(DAL,fieldname)
% Some extra error checks
errstr = '';
if (isequal(fieldname,'Inloop'))
   if (isempty(DAL.Inloop.params.main.attens))
      errstr = 'Attenuations are not set correctly! (high vs. low mismatch?)';
   end
   if (isempty(DAL.Inloop.params.main.tone.freq))
      errstr = 'Tone Frequency is empty!)';
   end
end

%----------------------------------------------------------------------------------------
function tmplt = template_definition(fieldname)
persistent prev_unit_bf prev_unit_thresh
%% DEFS: {Value Units Allowed-range Locked ??}

%%%%%%%%%%%%%%%%%%%%
%% Inloop Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Inloop.Frequency         =  {'round(current_unit_bf*1e4)/1e4'   'kHz'      [0.04  45]   0  0}; 
IO_def.Inloop.Low_Attenuation   =  { 1                  'dB'       [0    120]};
IO_def.Inloop.High_Attenuation  =  {'min(120,current_unit_thresh+20)'  'dB'   [0    120] 0 0};

if (~isequal(current_unit_bf, prev_unit_bf) & (isempty(fieldname) | isequal(fieldname,'Inloop')))
   IO_def.Inloop.Frequency{5}            = 1; % ignore dflt. Always recalculate.
   prev_unit_bf = current_unit_bf;
end
if (~isequal(current_unit_thresh, prev_unit_thresh) & (isempty(fieldname) | isequal(fieldname,'Inloop')))
   IO_def.Inloop.High_Attenuation{5}            = 1;
   prev_unit_thresh = current_unit_thresh;
end


%%%%%%%%%%%%%%%%%%%%
%% Gating Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Gating.Duration         = {200       'ms'    [20 2000] 1};
IO_def.Gating.Period           = {'default_period(this.Duration)'    'ms'   [50 5000] 1};
IO_def.Gating.Rise_fall_time   = {'default_rise_time(this.Duration)' 'ms'   [0  1000] 1}; 

%%%%%%%%%%%%%%%%%%%%
%% Mix Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Mix.Tone      =  {'Left|Both|{Right}' '' [] 1};

tmplt.tag         = 'MH_TB_tmplt';
tmplt.IO_def = IO_def;
