function [tmplt,DAL,stimulus_vals,units,errstr] = SP_template(fieldname,stimulus_vals,units)
%
% Template for Recruitment Stimulus SP: besh97k.wav, rate-level, 300-ms (TOTAL) duration
%
% MH 12/17/01

used_devices.File         = 'RP1.1';
tmplt = template_definition(fieldname);
if (exist('stimulus_vals','var') == 1)  
   Inloop.Name                         = 'DALinloop_wavfiles';
   Inloop.params.list                  = {stimulus_vals.Inloop.File};
   Inloop.params.attens                = stimulus_vals.Inloop.High_Attenuation :-1:stimulus_vals.Inloop.Low_Attenuation;
   Inloop.params.Rlist                 = [];
   Inloop.params.Rattens               = [];
   Inloop.params.repetitions           = 1;
   Inloop.params.resample_ratio        = NaN;
   Inloop.params.playback_slowdown     = 1;
   
   DAL.Inloop = Inloop;
   DAL.Gating = stimulus_vals.Gating;
   DAL.Mix         = mix_params2devs(stimulus_vals.Mix,used_devices);
   DAL.short_description   = 'SP';
   DAL.description = build_description(DAL,stimulus_vals);
   errstr = check_DAL_params(DAL,fieldname);
   
   %%%%%%%
   % If parameters are NOT correct for this template, Take away this template name
   [Xdir,Xfile,Xext]=fileparts(stimulus_vals.Inloop.File);
   if((~isequal(fullfile('',[Xfile Xext]),'besh97k.wav'))| ...
           (length(DAL.Inloop.params.attens) == 1))
      DAL.short_description   = '';
   end
end

%----------------------------------------------------------------------------------------
function str = build_description(DAL,stimulus_vals)
p = DAL.Inloop.params;
[fpath,file] = fileparts(stimulus_vals.Inloop.File);
str{1} = sprintf('File ''%s'' ', file);
if (length(p.attens) > 1)
   str{1} = sprintf('%s @ %1.1f - %1.1f dB Attn.', str{1}, p.attens(1), p.attens(end));
else
   if (~isempty(p.attens))
      str{1} = sprintf('%s @ %1.1f dB Attn.', str{1}, p.attens(1));
   end
end
str{1} = sprintf('%s (%s)', str{1}, stimulus_vals.Mix.File);

%----------------------------------------------------------------------------------------
function errstr = check_DAL_params(DAL,fieldname)
% Some extra error checks
errstr = '';
if (isequal(fieldname,'Inloop'))
   if (isempty(DAL.Inloop.params.attens))
      errstr = 'Attenuations are not set correctly! (high vs. low mismatch?)';
   end
end

%----------------------------------------------------------------------------------------
function tmplt = template_definition(fieldname)
global signals_dir
persistent prev_unit_thresh
%%%%%%%%%%%%%%%%%%%%
%% Inloop Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Inloop.File              = { fullfile(signals_dir,'besh97k.wav') '' [] 1 };
IO_def.Inloop.Low_Attenuation   = { 1             'dB'    [0    120]      };
IO_def.Inloop.High_Attenuation  = { 'current_unit_thresh'     'dB'    [0    120]  0  0 };

if (~isequal(current_unit_thresh, prev_unit_thresh) & (isempty(fieldname) | isequal(fieldname,'Inloop')))
   IO_def.Inloop.High_Attenuation{5}            = 1;
   prev_unit_thresh = current_unit_thresh;
end

%%%%%%%%%%%%%%%%%%%%
%% Gating Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Gating.Duration             = {300       'ms'    [20 2000] 1};
IO_def.Gating.Period               = {'default_period(this.Duration)'    'ms'   [50 5000] 1};


%%%%%%%%%%%%%%%%%%%%
%% Mix Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Mix.File        =  {'Left|Both|{Right}' '' [] 1};

tmplt.tag               = 'MH_SP_tmplt';
tmplt.IO_def = IO_def;
