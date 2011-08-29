function [tmplt,DAL,stimulus_vals,units,errstr] = SAMtone_IN_template(fieldname,stimulus_vals,units)


% used_devices.Tone       = 'RP1.1';
used_devices.Llist         = 'L3';
used_devices.Rlist         = 'R3';
tmplt = template_definition(fieldname);
if (exist('stimulus_vals','var') == 1)
   
   global signals_dir
   
   SAMsignals_dir=strcat(signals_dir,'MH\SAMtones');
   %% later also pass the modulation depth as a parameter
   [tone,fs,filename]=amtone(stimulus_vals.Inloop.Carrfreq*1000,stimulus_vals.Inloop.Modfrequency*1000,stimulus_vals.Inloop.Moddepth);
   
   stimulus_vals.Inloop.Compiled_FileName=fullfile(SAMsignals_dir,filename);
   stimulus_vals.Inloop.Noise_FileName=fullfile(SAMsignals_dir,'baseNOISE.wav');
   
   wavwrite(tone,fs,fullfile(SAMsignals_dir,filename));
   Llist = {stimulus_vals.Inloop.Compiled_FileName};
   Rlist = {stimulus_vals.Inloop.Noise_FileName};
   
   [stimulus_vals units] = NI_check_gating_params(stimulus_vals, units);
   [stimulus_vals.Mix units.Mix] = structdlg(tmplt.IO_def.Mix,'',stimulus_vals.Mix,'off');
   
   if ~isempty(Llist)&~isempty(Rlist)  % Two files used on NI board, max sampling rate is HALF!
      MAX_NI_SAMPLING_RATE=NI6052UsableRate_Hz(Inf)/2;
   else  % Only 1 file, can use full rate
      MAX_NI_SAMPLING_RATE=NI6052UsableRate_Hz(Inf);
   end   
   
   if (isempty(Llist))
      tmplt.IO_def.Mix = rmfield(tmplt.IO_def.Mix,'Llist');
   end
   if (isempty(Rlist))
      tmplt.IO_def.Mix = rmfield(tmplt.IO_def.Mix,'Rlist');
   end
   [stimulus_vals.Mix units.Mix] = structdlg(tmplt.IO_def.Mix,'',stimulus_vals.Mix,'off');
   
%    %% Set maximum sampling rate based on how many wavfiles to be played
%    if ~isempty(Llist)&~isempty(Rlist)  % Two files used on NI board, max sampling rate is HALF!
%       MAX_NI_SAMPLING_RATE=NI6052UsableRate_Hz(Inf)/2;
%    else  % Only 1 file, can use full rate
%       MAX_NI_SAMPLING_RATE=NI6052UsableRate_Hz(Inf);
%    end   
   
   stimulus_vals.Inloop.Used_UpdateRate_Hz = NI6052UsableRate_Hz(fs);
%    stimulus_vals.Inloop.Used_UpdateRate_Hz=fs;
   
   %Inloop.Name                               = 'DALinloop_wavfiles2';
   Inloop.Name                                = 'DALinloop_NI_SCC_wavfiles2';
   
   Inloop.params.list                        = Llist;
   Inloop.params.Rlist                       = Rlist;

   Inloop.params.Rattens                      = stimulus_vals.Inloop.High_Attenuation_noise :-stimulus_vals.Inloop.dBstep_Atten_noise:stimulus_vals.Inloop.Low_Attenuation_noise;
   Inloop.params.attens                      = stimulus_vals.Inloop.Tone_Attenuation*ones(size(Inloop.params.Rattens));  % Needed to run all levels
   %Inloop.params.resample_ratio              = [];
   %Inloop.params.playback_slowdown           = [];
   Inloop.params.condition.Modfrequency      = stimulus_vals.Inloop.Modfrequency;
   Inloop.params.condition.Moddepth          = stimulus_vals.Inloop.Moddepth;
   Inloop.params.condition.Carrfreq          = stimulus_vals.Inloop.Carrfreq;
   Inloop.params.repetitions                 = stimulus_vals.Inloop.Repetitions;
   Inloop.params.updateRate_Hz               = stimulus_vals.Inloop.Used_UpdateRate_Hz;

   DAL.funcName = 'data_acquisition_loop_NI'; 
   DAL.Inloop = Inloop;
   DAL.Gating = stimulus_vals.Gating;
   DAL.short_description   = 'SAM_IN';

  
   DAL.Mix         = mix_params2devs(stimulus_vals.Mix,used_devices);

   DAL.description = build_description(DAL,stimulus_vals);
   errstr = check_DAL_params(DAL,fieldname);
   
   %%%%%%%
   % If parameters are NOT correct for this template, Take away this template name
%    if((stimulus_vals.Inloop.Frequency ~= 2)|(length(DAL.Inloop.params.main.attens) == 1))
%       DAL.short_description   = '';
%    end
   %latest_user_attn(stimulus_vals.Inloop.High_Attenuation);
end

%----------------------------------------------------------------------------------------
function str = build_description(DAL,stimulus_vals)
% p = DAL.Inloop.params;
p=stimulus_vals.Inloop;
str{1} = sprintf('%1.2f kHz Tone', p.Carrfreq);
% if (length(p.main.attens) > 1)
%    str{1} = sprintf('%s @ %1.1f - %1.1f dB Attn.', str{1}, p.main.attens(1), p.main.attens(end));
% else
%    str{1} = sprintf('%s @ %1.1f dB Attn.', str{1}, p.main.attens(1));
% end
% str{1} = sprintf('%s (%s)', str{1}, stimulus_vals.Mix.Tone);

%----------------------------------------------------------------------------------------
function errstr = check_DAL_params(DAL,fieldname)
% Some extra error checks
errstr = '';
% if (isequal(fieldname,'Inloop'))
%    if (isempty(DAL.Inloop.params.main.attens))
%       errstr = 'Attenuations are not set correctly! (high vs. low mismatch?)';
%    end
%    if (isempty(DAL.Inloop.params.main.tone.freq))
%       errstr = 'Tone Frequency is empty!)';
%    end
% end

%----------------------------------------------------------------------------------------
function tmplt = template_definition(fieldname)
%% DEFS: {Value Units Allowed-range Locked ??}

%%%%%%%%%%%%%%%%%%%%
%% Inloop Section 
%%%%%%%%%%%%%%%%%%%%

IO_def.Inloop.Modfrequency      =  { 0.05             'kHz'      [0.01  2] 0}; 
IO_def.Inloop.Moddepth          =  { 1             ''         [0.01  1]   0}; 
IO_def.Inloop.Carrfreq          =  {'current_unit_bf' 'kHz'      [0.04  50]  0 0}; 
% IO_def.Inloop.Level             =  { 65               'dB SPL'   [-10  120]  0}; 
IO_def.Inloop.Tone_Attenuation  =  {'min(120,current_unit_thresh-20)'  'dB'   [0    120] 0 0};  % 20 dB above threshold

IO_def.Inloop.Low_Attenuation_noise   = { 30             'dB'    [0    120]      };  % Here, these are for the noise
IO_def.Inloop.High_Attenuation_noise  = { 80           'dB'    [0    120]      };
IO_def.Inloop.dBstep_Atten_noise      = {10           'dB'    [1    120]      };
% IO_def.Inloop.Noise_Attenuation =  {10               'dB'    [0    120]      };

IO_def.Inloop.Repetitions       =  {25                 ''       [1 600]};
%%%%%%%%%%%%%%%%%%%%
%% Gating Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Gating.Duration         = {400       'ms'    [20 2000] 1};
IO_def.Gating.Period           = {'default_period(this.Duration)'    'ms'   [50 5000] 1};
IO_def.Gating.Rise_fall_time   = {'default_rise_time(this.Duration)' 'ms'   [0  1000] 1}; 

%%%%%%%%%%%%%%%%%%%%
%% Mix Section 
%%%%%%%%%%%%%%%%%%%%
%IO_def.Mix.Tone      =  {'Left|Both|{Right}' '' [] 1};
IO_def.Mix.Llist        =  {'Left|Both|{Right}'};
IO_def.Mix.Rlist        =  {'Left|Both|{Right}'};

tmplt.tag         = 'SK_SAM_IN_tmplt';
tmplt.IO_def = IO_def;
