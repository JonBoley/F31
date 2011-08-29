function [tmplt,DAL,stimulus_vals,units,errstr] = TT_rate_level_template(fieldname,stimulus_vals,units)
%

% AF 11/26/01

used_devices.Tone         = 'RP1.1';
used_devices.Fixed_Tone   = 'RP2.1';
tmplt = template_definition(fieldname);
if (exist('stimulus_vals','var') == 1)
   Inloop.Name                               = 'DALinloop_general_TN';
   Inloop.params.main.source                 = 'Tone';
   Inloop.params.main.tone.freq              = stimulus_vals.Inloop.Frequency*1000;
   Inloop.params.main.tone.bw                = 0;
   Inloop.params.main.noise.low_cutoff       = 0;
   Inloop.params.main.noise.high_cutoff      = 0;
   Inloop.params.main.attens                 = stimulus_vals.Inloop.High_Attenuation :-1:stimulus_vals.Inloop.Low_Attenuation;
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
   DAL.short_description   = 'RLV_2T';

   DAL.Mix         = mix_params2devs(stimulus_vals.Mix,used_devices);

   DAL.description = build_description(DAL,stimulus_vals);
   errstr = check_DAL_params(DAL,fieldname);
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
   str{1} = sprintf('%s @ %1.1f dB Attn.', str{1}, p.main.attens(1));
end
str{1} = sprintf('%s (%s)', str{1}, stimulus_vals.Mix.Tone);
str{2} = sprintf('Background %1.2f kHz Tone @ %1.1f dB Attn.', p.secondary.tone.freq/1000, p.secondary.atten);
str{2} = sprintf('%s (%s)', str{2}, stimulus_vals.Mix.Fixed_Tone);
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
IO_def.Inloop.Frequency              = {'current_unit_bf'       'kHz'   [0.04  50]  0  0}; 
IO_def.Inloop.Low_Attenuation        = { 1                      'dB'    [0    120]      };
IO_def.Inloop.High_Attenuation       = { 100                    'dB'    [0    120]      };
IO_def.Inloop.Fixed_Tone_Frequency   = {'current_unit_bf'       'kHz'   [0.04  50]  0  0};
IO_def.Inloop.Fixed_Tone_Attenuation = {'current_unit_thresh-3' 'dB'    [0    120]  0  0};

if (isequal(fieldname,'Inloop'))
   if (~isequal(current_unit_bf, prev_unit_bf))
      IO_def.Inloop.Frequency{5}            = 1; % ignore dflt. Always recalculate.
      IO_def.Inloop.Fixed_Tone_Frequency{5} = 1;
      prev_unit_bf = current_unit_bf;
   end
   if (~isequal(current_unit_thresh, prev_unit_thresh))
      IO_def.Inloop.Fixed_Tone_Attenuation{5} = 1;% ignore dflt. Always recalculate.
      prev_unit_thresh = current_unit_thresh;
   end
end

%%%%%%%%%%%%%%%%%%%%
%% Gating Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Gating.Duration         = {200       'ms'    [20 2000]};
IO_def.Gating.Period           = {'default_period(this.Duration)'    'ms'   [50 5000]};
IO_def.Gating.Rise_fall_time   = {'default_rise_time(this.Duration)' 'ms'   [0  1000]}; 

%%%%%%%%%%%%%%%%%%%%
%% Mix Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Mix.Tone        =  {'{Left}|Both|Right'};
IO_def.Mix.Fixed_Tone  =  {'{Left}|Both|Right'};

tmplt.tag         = 'LRAFrate_level_2T_tmplt';
tmplt.IO_def = IO_def;
