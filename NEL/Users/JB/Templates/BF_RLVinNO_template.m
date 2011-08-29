function [tmplt,DAL,stimulus_vals,units,errstr] = BF_RLVinNO_template(fieldname,stimulus_vals,units)
%
% AF 11/26/01
%
% Modified by M.Heinz 04Dec2003, from nel_TT_rate_level_template.m

used_devices.Tone         = 'RP1.1';
used_devices.Fixed_Noise   = 'RP2.1';
tmplt = template_definition(fieldname);
if (exist('stimulus_vals','var') == 1)
   Inloop.Name                               = 'DALinloop_general_TN';
   Inloop.params.main.source                 = 'Tone';
   Inloop.params.main.tone.freq              = stimulus_vals.Inloop.Frequency*1000;
   Inloop.params.main.tone.bw                = 0;
   Inloop.params.main.noise.low_cutoff       = 0;
   Inloop.params.main.noise.high_cutoff      = 0;
   Inloop.params.main.attens                 = stimulus_vals.Inloop.High_Attenuation :-1:stimulus_vals.Inloop.Low_Attenuation;
   Inloop.params.secondary.source            = 'Noise';
   Inloop.params.secondary.tone.freq         = 0;
   Inloop.params.secondary.noise.low_cutoff  = stimulus_vals.Inloop.Fixed_Noise_Low_Cutoff*1000;
   Inloop.params.secondary.noise.high_cutoff = stimulus_vals.Inloop.Fixed_Noise_High_Cutoff*1000;
   Inloop.params.secondary.noise.gating      = 'Positive';
   Inloop.params.secondary.noise.adaptation  = 0;
   Inloop.params.secondary.atten             = stimulus_vals.Inloop.Fixed_Noise_Attenuation;
   Inloop.params.rise_fall                   = stimulus_vals.Gating.Rise_fall_time;

   DAL.Inloop = Inloop;
   DAL.Gating = stimulus_vals.Gating;
   DAL.short_description   = 'BFNO';

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
str{2} = sprintf('Fixed-Noise Masker:  %1.2f - %1.2f @ %1.1f dB Attn.', p.secondary.noise.low_cutoff/1000,p.secondary.noise.high_cutoff/1000, p.secondary.atten);
str{2} = sprintf('%s (%s)', str{2}, stimulus_vals.Mix.Fixed_Noise);
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
    if (DAL.Inloop.params.secondary.noise.low_cutoff >= DAL.Inloop.params.secondary.noise.high_cutoff)
        errstr = 'Fixed-Noise-Masker Frequency Band not set correctly!)';
    end
    if (isempty(DAL.Inloop.params.secondary.atten))
        errstr = 'Fixed-Noise-Masker Attenuation is not set!';
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
IO_def.Inloop.Fixed_Noise_Low_Cutoff   = { 0      'kHz'   [0  50]  1};
IO_def.Inloop.Fixed_Noise_High_Cutoff   = { 50       'kHz'    [0  50]  1};
IO_def.Inloop.Fixed_Noise_Attenuation = { 'current_unit_thresh-30' 'dB'    [0    120] 0  0};

if (isequal(fieldname,'Inloop'))
   if (~isequal(current_unit_bf, prev_unit_bf))
      IO_def.Inloop.Frequency{5}            = 1; % ignore dflt. Always recalculate.
      prev_unit_bf = current_unit_bf;
   end
   if (~isequal(current_unit_thresh, prev_unit_thresh))
      IO_def.Inloop.Fixed_Noise_Attenuation{5} = 1;% ignore dflt. Always recalculate.
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
IO_def.Mix.Fixed_Noise  =  {'Left|Both|{Right}'};

tmplt.tag         = 'BF_RLVinNO_tmplt';
tmplt.IO_def = IO_def;
