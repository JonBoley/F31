function [tmplt,DAL,stimulus_vals,units,errstr] = SAM_TMDF_template(fieldname,stimulus_vals,units)



% used_devices.Tone       = 'RP1.1';
used_devices.list         = 'L3';
tmplt = template_definition(fieldname);
if (exist('stimulus_vals','var') == 1)
   
   global signals_dir
   SAMsignals_dir=strcat(signals_dir,'SK\SAMtones');
   
   mdMin=stimulus_vals.Inloop.mdMin; mdMax=stimulus_vals.Inloop.mdMax; 
   no_of_mdPoints=stimulus_vals.Inloop.no_of_mdPoints;
   
   if strcmp(stimulus_vals.Inloop.mdScale_type,'logarithmic')
       mdScale=logspace(log10(mdMin),log10(mdMax),no_of_mdPoints);
   elseif strcmp(stimulus_vals.Inloop.mdScale_type,'linear') 
       mdScale=linspace(mdMin,mdMax,no_of_mdPoints);
   end    
stimulus_vals.Inloop.mdScale=mdScale;
   for i=1:length(mdScale)
       [tone,fs,filename]=amtone(stimulus_vals.Inloop.Carrfreq*1000,stimulus_vals.Inloop.Modfreq*1000,mdScale(i));
       stimulus_vals.Inloop.Compiled_FileName=fullfile(SAMsignals_dir,filename);
       list{i}=stimulus_vals.Inloop.Compiled_FileName;
       wavwrite(tone,fs,fullfile(SAMsignals_dir,filename));
   end    
  
   
   [stimulus_vals units] = NI_check_gating_params(stimulus_vals, units);
   [stimulus_vals.Mix units.Mix] = structdlg(tmplt.IO_def.Mix,'',stimulus_vals.Mix,'off');
   

   
   if (isempty(list))
      tmplt.IO_def.Mix = rmfield(tmplt.IO_def.Mix,'list');
   end
%    if (isempty(Rlist))
%       tmplt.IO_def.Mix = rmfield(tmplt.IO_def.Mix,'Rlist');
%    end
   [stimulus_vals.Mix units.Mix] = structdlg(tmplt.IO_def.Mix,'',stimulus_vals.Mix,'off');
   
   %% Set maximum sampling rate based on how many wavfiles to be played
%    if ~isempty(Llist)&~isempty(Rlist)  % Two files used on NI board, max sampling rate is HALF!
%       MAX_NI_SAMPLING_RATE=NI6052UsableRate_Hz(Inf)/2;
%    else  % Only 1 file, can use full rate
%       MAX_NI_SAMPLING_RATE=NI6052UsableRate_Hz(Inf);
%    end   
   
   stimulus_vals.Inloop.Used_UpdateRate_Hz = NI6052UsableRate_Hz(fs);
%    stimulus_vals.Inloop.Used_UpdateRate_Hz=fs;
   
   %Inloop.Name                               = 'DALinloop_wavfiles2';
   Inloop.Name                                = 'DALinloop_NI_SCC_wavfiles2';
   
   Inloop.params.list                        = list;
   Inloop.params.Rlist                       = []; 
   Inloop.params.Rattens                     = []; 
   Inloop.params.attens                      = stimulus_vals.Inloop.Tone_Attenuation;
   
   %Inloop.params.Condition.Level_dBSPL       = stimulus_vals.Inloop.Level; %changed on 06/25/2007
   Inloop.params.condition.mdScale           = stimulus_vals.Inloop.mdScale;
   Inloop.params.condition.Modfreq          = stimulus_vals.Inloop.Modfreq;
   Inloop.params.condition.Carrfreq          = stimulus_vals.Inloop.Carrfreq;
   Inloop.params.repetitions                 = stimulus_vals.Inloop.Repetitions;
   Inloop.params.updateRate_Hz               = stimulus_vals.Inloop.Used_UpdateRate_Hz;

   DAL.funcName = 'data_acquisition_loop_NI'; 
   DAL.Inloop = Inloop;
   DAL.Gating = stimulus_vals.Gating;
   DAL.short_description   = 'SAM_TMDF';

  
   DAL.Mix         = mix_params2devs(stimulus_vals.Mix,used_devices);

   DAL.description = build_description(DAL,stimulus_vals);
   errstr = check_DAL_params(DAL,fieldname);
   
  
end

%----------------------------------------------------------------------------------------
function str = build_description(DAL,stimulus_vals)
% p = DAL.Inloop.params;
p=stimulus_vals.Inloop;
str{1} = sprintf('%1.2f kHz fc, %1.2f kHz fm, %.2f-%.2f md, %.f dBatten', p.Carrfreq,p.Modfreq,p.mdMin,p.mdMax,p.Tone_Attenuation);
%str{1} = sprintf('%1.2f kHz fm,%1.2f kHz fc,%1.2f m', p.Carrfreq,p.fmScale.p.Moddepth);
% str{1} = sprintf('%1.2f kHz fm', p.Modfrequency);
% str{1} = sprintf('%1.2f kHz m', p.Moddepth);
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
global signals_dir
persistent prev_unit_bf prev_unit_thresh
%%%%%%%%%%%%%%%%%%%%
%% Inloop Section 
%%%%%%%%%%%%%%%%%%%%

IO_def.Inloop.Tone_Attenuation  =  {'min(120,current_unit_thresh-10)'  'dB'   [0    120] 0 0};  % 20 dB above threshold
IO_def.Inloop.Modfreq          =  { 0.05             'kHz'      [0.01  1.5]   0}; 
IO_def.Inloop.Carrfreq          =  {'current_unit_bf' 'kHz'     [0.04  50]  0 0}; 
IO_def.Inloop.Repetitions       =  {25                 ''       [1 600]};
IO_def.Inloop.mdScale_type      =  {'{logarithmic}|linear'};
IO_def.Inloop.mdMin             =  {0.1         ''    [0 1]};
IO_def.Inloop.mdMax             =  {1        ''       [0 1]};
IO_def.Inloop.no_of_mdPoints    =  {4          ''     [1 100]};

%IO_def.Inloop.Noise_Attenuation =  {10               'dB'    [0    120]      };
if (~isequal(current_unit_bf, prev_unit_bf) & isequal(fieldname,'Inloop'))
   IO_def.Inloop.Carrfreq{5}            = 1; % ignore dflt. Always recalculate.
   prev_unit_bf = current_unit_bf;
end
if (~isequal(current_unit_thresh, prev_unit_thresh) & isequal(fieldname,'Inloop'))
   IO_def.Inloop.Tone_Attenuation{5}            = 1; % ignore dflt. Always recalculate.
   prev_unit_thresh = current_unit_thresh;
end
%%%%%%%%%%%%%%%%%%%%
%% Gating Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Gating.Duration         = {600       'ms'    [20 2000] 1};
% IO_def.Gating.Period           = {'default_period(this.Duration)'    'ms'   [50 5000] 1};
IO_def.Gating.Period           = {1000    'ms'   [50 5000] 1};
IO_def.Gating.Rise_fall_time   = {'default_rise_time(this.Duration)' 'ms'   [0  1000] 1}; 

%%%%%%%%%%%%%%%%%%%%
%% Mix Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Mix.list        =  {'Left|Both|{Right}'};

tmplt.tag         = 'SK_SAM_TMDF_tmplt';
tmplt.IO_def = IO_def;
