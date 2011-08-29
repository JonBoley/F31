function [tmplt,DAL,stimulus_vals,units,errstr] = B2_template(fieldname,stimulus_vals,units)
%
% Template for Recruitment Stimulus B2: 2-kHz tone in BP (1.8-3.0 kHz) noise, rate-level, 200-ms duration
%
% MH 12/17/01

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
   DAL.short_description   = 'B2';

   DAL.Mix         = mix_params2devs(stimulus_vals.Mix,used_devices);

   DAL.description = build_description(DAL,stimulus_vals);
   errstr = check_DAL_params(DAL,fieldname);

   %%%%%%%
   % If parameters are NOT correct for this template, Take away this template name
   if((stimulus_vals.Inloop.Frequency ~= 2)| ...
           (length(DAL.Inloop.params.main.attens) == 1)| ...
           (stimulus_vals.Inloop.Fixed_Noise_Low_Cutoff ~= 1.8)| ...
           (stimulus_vals.Inloop.Fixed_Noise_High_Cutoff ~= 3))
      DAL.short_description   = '';
   end
   latest_user_attn(stimulus_vals.Inloop.High_Attenuation);
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
%persistent prev_unit_bf prev_unit_thresh
%%%%%%%%%%%%%%%%%%%%
%% Inloop Section 
%%%%%%%%%%%%%%%%%%%%
hi_attn = latest_user_attn;
if (isempty(hi_attn))
    hi_attn = 100;
end
IO_def.Inloop.Frequency              = { 2       'kHz'   [0.04  45]  1}; 
IO_def.Inloop.Low_Attenuation        = { 1                      'dB'    [0    120]      };
IO_def.Inloop.High_Attenuation       = { hi_attn                'dB'    [0    120]      };
IO_def.Inloop.Fixed_Noise_Low_Cutoff   = { 1.8      'kHz'   [0  50]  1};
IO_def.Inloop.Fixed_Noise_High_Cutoff   = { 3       'kHz'    [0  50]  1};
IO_def.Inloop.Fixed_Noise_Attenuation = { [] 'dB'    [0    120]};

%%%%%%%%%%%%%%%%%%%%
%% Gating Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Gating.Duration         = {200       'ms'    [20 2000] 1};
IO_def.Gating.Period           = {'default_period(this.Duration)'    'ms'   [50 5000] 1};
IO_def.Gating.Rise_fall_time   = {'default_rise_time(this.Duration)' 'ms'   [0  1000] 1}; 

%%%%%%%%%%%%%%%%%%%%
%% Mix Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Mix.Tone        =  {'Left|Both|{Right}' '' [] 1};
IO_def.Mix.Fixed_Noise  =  {'Left|Both|{Right}' '' [] 1};

tmplt.tag         = 'MH_B2_tmplt';
tmplt.IO_def = IO_def;
