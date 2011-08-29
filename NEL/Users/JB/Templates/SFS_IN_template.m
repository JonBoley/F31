function [tmplt,DAL,stimulus_vals,units,errstr] = SFS_IN_template(fieldname,stimulus_vals,units)


% used_devices.Tone       = 'RP1.1';
used_devices.Llist         = 'L3';
used_devices.Rlist         = 'R3';
tmplt = template_definition(fieldname);
if (exist('stimulus_vals','var') == 1)
   
   global signals_dir
   
   SFSsignals_dir=strcat(signals_dir,'SK\SFvowel');
   fs=81920;
   [vowel,dBreTONE,filename]= sfs_synth(stimulus_vals.Inloop.Carrfreq*1000,stimulus_vals.Inloop.f0*1000,stimulus_vals.Gating.Duration/1000,fs);
   
   stimulus_vals.Inloop.Compiled_FileName=fullfile(SFSsignals_dir,filename);
   stimulus_vals.Inloop.Noise_FileName=fullfile(SFSsignals_dir,'baseNOISE.wav');
   
   wavwrite(vowel,fs,fullfile(SFSsignals_dir,filename));
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
   
   Inloop.params.Rattens                     = stimulus_vals.Inloop.High_Attenuation_noise :-stimulus_vals.Inloop.dBstep_Atten_noise:stimulus_vals.Inloop.Low_Attenuation_noise;
   Inloop.params.attens                      = stimulus_vals.Inloop.Vowel_Attenuation*ones(size(Inloop.params.Rattens));  % Needed to run all levels
   
   
   Inloop.params.condition.f0                = stimulus_vals.Inloop.f0;
   Inloop.params.condition.Carrfreq          = stimulus_vals.Inloop.Carrfreq;
   
   Inloop.params.repetitions                 = stimulus_vals.Inloop.Repetitions;
   Inloop.params.updateRate_Hz               = stimulus_vals.Inloop.Used_UpdateRate_Hz;

   DAL.funcName = 'data_acquisition_loop_NI'; 
   DAL.Inloop = Inloop;
   DAL.Gating = stimulus_vals.Gating;
   DAL.short_description   = 'SFS_IN';

  
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
% str{1} = sprintf('%1.2f kHz Tone', p.Carrfreq);
str{1} = sprintf('%1.2f kHz fc,%1.2f kHz f0,%.f dBatten', p.Carrfreq,p.f0,p.Vowel_Attenuation);
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
global signals_dir
persistent prev_unit_bf prev_unit_thresh
%%%%%%%%%%%%%%%%%%%%
%% Inloop Section 
%%%%%%%%%%%%%%%%%%%%

IO_def.Inloop.f0                =  { 0.1             'kHz'      [0.1  0.3] 0}; 
IO_def.Inloop.Carrfreq          =  {'current_unit_bf' 'kHz'      [0.04  50]  0 0}; 
% IO_def.Inloop.Duration          =  {200  'ms'      [50  500]  0 0}; 

IO_def.Inloop.Vowel_Attenuation       = {'min(120,current_unit_thresh-30)'  'dB'   [0    120] 0 0};  % 20 dB above threshold
IO_def.Inloop.Low_Attenuation_noise   = { 30             'dB'    [0    120]      };  % Here, these are for the noise
IO_def.Inloop.High_Attenuation_noise  = { 80           'dB'    [0    120]      };
IO_def.Inloop.dBstep_Atten_noise      = {5           'dB'    [1    120]      };

IO_def.Inloop.Repetitions       =  {25                 ''       [1 600]};

if (~isequal(current_unit_bf, prev_unit_bf) & isequal(fieldname,'Inloop'))
   IO_def.Inloop.Carrfreq{5}            = 1; % ignore dflt. Always recalculate.
   prev_unit_bf = current_unit_bf;
end
if (~isequal(current_unit_thresh, prev_unit_thresh) & isequal(fieldname,'Inloop'))
   IO_def.Inloop.Vowel_Attenuation{5}            = 1; % ignore dflt. Always recalculate.
   prev_unit_thresh = current_unit_thresh;
end

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

tmplt.tag         = 'SK_SFS_IN_tmplt';
tmplt.IO_def = IO_def;
