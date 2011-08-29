function [tmplt,DAL,stimulus_vals,units,errstr] = FW_template(fieldname,stimulus_vals,units)
%
% Template for T. Ji Exps, "Five Women ..." sentence: fivewo.wav 
%
% MH: 7/31/02

persistent   prev_playdur  prev_min_period
% We use the persistent variables to detect a change that requires some fields update.
% For example, of the play duration is changed we would like to update the gating information.
% We restict the automatic updates to allow the user the overide them.

used_devices.File         = 'RP1.1';
tmplt = template_definition(fieldname);
if (exist('stimulus_vals','var') == 1)    
   if (exist(stimulus_vals.Inloop.File,'file') ~= 0)
      list = {stimulus_vals.Inloop.File};
      [data fs] = wavread(stimulus_vals.Inloop.File);
      playback_rate = 97656.25/(round(97656.25/fs)); % 97656.25 is the rco's sampling rate
      playdur = round(length(data)/playback_rate*1000);	%compute file duration based on sampling rate
      min_period = max(1000,ceil(1.7*playdur/100)*100);
      
      %% In this template we have to change the Gating parameters according to 
      %% the playback duration of the wav files.
      %% We do this by first, updating the relevant template definitions, and second,
      %% by calling 'structdlg' in its non-interactive invisible mode, to recalculated 'stimulus_vals' fields.
      %
      if (isequal(fieldname,'Inloop') | ~isequal(playdur,prev_playdur) | ~isequal(min_period,prev_min_period))
         tmplt.IO_def.Gating.Duration{1}  = playdur;
         tmplt.IO_def.Gating.Period{1}    = playdur+1000;
         %         tmplt.IO_def.Gating.Period{1}    = ['max(' num2str(min_period) ',default_period(this.Duration))'];
         [stimulus_vals.Gating units.Gating] = structdlg(tmplt.IO_def.Gating,'',[],'off');
         prev_playdur = playdur;
         prev_min_period = min_period;
      end
   else
      list = {};
   end
     
   Inloop.Name                         = 'DALinloop_wavfiles';
   Inloop.params.list                  = {stimulus_vals.Inloop.File};
   Inloop.params.attens                = stimulus_vals.Inloop.Attenuation;
   Inloop.params.Rlist                 = [];
   Inloop.params.Rattens               = [];
   Inloop.params.repetitions           = stimulus_vals.Inloop.Repetitions;
   Inloop.params.resample_ratio        = NaN;
   Inloop.params.playback_slowdown     = 1;
   
   DAL.Inloop = Inloop;
   DAL.Gating = stimulus_vals.Gating;
   DAL.Mix         = mix_params2devs(stimulus_vals.Mix,used_devices);
   DAL.short_description   = 'FW';
   DAL.description = build_description(DAL,stimulus_vals);
   errstr = check_DAL_params(DAL,fieldname);
   
   %%%%%%%
   % If parameters are NOT correct for this template, Take away this template name
   [Xdir,Xfile,Xext]=fileparts(stimulus_vals.Inloop.File);
   if((~isequal(fullfile('',[Xfile Xext]),'fivewo.wav'))| ...
         (length(DAL.Inloop.params.attens) ~= 1))
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
      errstr = 'Attenuation is not set!';
   end
   if (length(DAL.Inloop.params.attens) > 1)
      errstr = 'Only one attenuation please';
   end
end

%----------------------------------------------------------------------------------------
function tmplt = template_definition(fieldname)
global signals_dir
%%%%%%%%%%%%%%%%%%%%
%% Inloop Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Inloop.File              = { fullfile(signals_dir,'fivewo.wav') '' [] 1 };
IO_def.Inloop.Attenuation   = { []             'dB'    [0    120]      };
IO_def.Inloop.Repetitions            = { 100                        ''      [1    Inf]      };


%%%%%%%%%%%%%%%%%%%%
%% Gating Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Gating.Duration             = {1601       'ms'    [20 2000] 1};
%IO_def.Gating.Period               = {'default_period(this.Duration)'    'ms'   [50 5000] 1};
IO_def.Gating.Period               = {2600    'ms'   [50 5000] 1};


%%%%%%%%%%%%%%%%%%%%
%% Mix Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Mix.File        =  {'Left|Both|{Right}' '' [] 1};

tmplt.tag               = 'FW_tmplt';
tmplt.IO_def = IO_def;
