function [tmplt,DAL,stimulus_vals,units,errstr] = RM2T_template(fieldname,stimulus_vals,units)
%
% Modified by M.Heinz 03Dec2003, from nel_TT_resp_map_template.m

% AF 11/26/01

used_devices.Tone       = 'RP1.1';
used_devices.Fixed_Tone   = 'RP2.1';
tmplt = template_definition(fieldname);
if (exist('stimulus_vals','var') == 1)
   Inloop.Name                               = 'DALinloop_general_TN';
   Inloop.params.main.source                 = 'Tone';
   Inloop.params.main.tone.freq              = stimulus_vals.Inloop.Frequency*1000;
   Inloop.params.main.tone.bw                = stimulus_vals.Inloop.Bandwidth;;
   Inloop.params.main.noise.low_cutoff       = 0;
   Inloop.params.main.noise.high_cutoff      = 0;
   Inloop.params.main.attens                 = stimulus_vals.Inloop.Attenuation;
   Inloop.params.secondary.source            = 'Tone';
   Inloop.params.secondary.tone.freq         = stimulus_vals.Inloop.Fixed_Tone_Frequency*1000;
   Inloop.params.secondary.noise.low_cutoff  = 0;
   Inloop.params.secondary.noise.high_cutoff = 0;
   Inloop.params.secondary.noise.gating      = '';
   Inloop.params.secondary.noise.adaptation  = 0;
   Inloop.params.secondary.atten             = stimulus_vals.Inloop.Fixed_Tone_Attenuation;
   Inloop.params.rise_fall                   = stimulus_vals.Gating.Rise_fall_time;

   DAL.Inloop = Inloop;
   DAL.Gating = stimulus_vals.Gating;
   DAL.short_description   = 'RM_2T';

   DAL.Mix         = mix_params2devs(stimulus_vals.Mix,used_devices);

   DAL.description = build_description(DAL,stimulus_vals);
   errstr = check_DAL_params(DAL,fieldname);
   if (DAL.Inloop.params.main.tone.bw == 0)
      DAL.short_description   = '';
   end
end

%----------------------------------------------------------------------------------------
function str = build_description(DAL,stimulus_vals)
p = DAL.Inloop.params;
str{1} = sprintf('%1.1f oct. sweep around %1.2f kHz', p.main.tone.bw, p.main.tone.freq/1000);
str{1} = sprintf('%s @ %1.1f dB Attn.', str{1}, p.main.attens);
str{1} = sprintf('%s (%s)', str{1}, stimulus_vals.Mix.Tone);
str{2} = sprintf('Background %1.2f kHz Tone @ %1.1f dB Attn.', p.secondary.tone.freq/1000, p.secondary.atten);
str{2} = sprintf('%s (%s)', str{2}, stimulus_vals.Mix.Fixed_Tone);

%----------------------------------------------------------------------------------------
function errstr = check_DAL_params(DAL,fieldname)
% Some extra error checks
errstr = '';
if (isequal(fieldname,'Inloop'))
   if (isempty(DAL.Inloop.params.main.attens))
      errstr = 'You must set an Attenuation levels';
   end
   if (isempty(DAL.Inloop.params.main.tone.freq))
      errstr = 'Tone Frequency is empty!)';
   end
   if (isempty(DAL.Inloop.params.secondary.tone.freq))
      errstr = 'Fixed Tone Frequency is empty!)';
   end
   if (isempty(DAL.Inloop.params.secondary.atten))
      errstr = 'Fixed Tone Attenuation is not set!';
   end
end

%----------------------------------------------------------------------------------------
function tmplt = template_definition(fieldname)
persistent prev_unit_bf prev_unit_thresh
%%%%%%%%%%%%%%%%%%%%
%% Inloop Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Inloop.Frequency              = {'current_unit_bf'        'kHz'   [0.04  50]   0  0}; % ignore dflt. Always recalculate.
IO_def.Inloop.Bandwidth              = { 4                       'oct'   [0    Inf]       };
IO_def.Inloop.Attenuation            = {'current_unit_thresh-20' 'dB'    [0    120]   0  0}; 
IO_def.Inloop.Fixed_Tone_Frequency   = {'current_unit_bf'        'kHz'   [0.04  50]   0  0}; % ignore dflt. Always recalculate.
IO_def.Inloop.Fixed_Tone_Attenuation = {'current_unit_thresh-3'  'dB'    [0    120]   0  0};

if (isequal(fieldname,'Inloop'))
   if (~isequal(current_unit_bf, prev_unit_bf))
      IO_def.Inloop.Frequency{5}            = 1; % ignore dflt. Always recalculate.
      IO_def.Inloop.Fixed_Tone_Frequency{5} = 1;
      prev_unit_bf = current_unit_bf;
   end
   if (~isequal(current_unit_thresh, prev_unit_thresh))
      IO_def.Inloop.Attenuation{5}            = 1;
      IO_def.Inloop.Fixed_Tone_Attenuation{5} = 1;% ignore dflt. Always recalculate.
      prev_unit_thresh = current_unit_thresh;
   end
end


%%%%%%%%%%%%%%%%%%%%
%% Gating Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Gating.Duration         = {50       'ms'    [20 2000]};
IO_def.Gating.Period           = {'default_period(this.Duration)'    'ms'   [50 5000]};
% IO_def.Gating.Rise_fall_time   = {'default_rise_time(this.Duration)' 'ms'   [0  1000]}; 
IO_def.Gating.Rise_fall_time   = {2 'ms'   [0  1000]}; 

%%%%%%%%%%%%%%%%%%%%
%% Mix Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Mix.Tone        =  {'Left|Both|{Right}'};
IO_def.Mix.Fixed_Tone  =  {'Left|Both|{Right}'};

tmplt.tag         = 'NELresponse_map_2T_tmplt';
tmplt.IO_def = IO_def;
