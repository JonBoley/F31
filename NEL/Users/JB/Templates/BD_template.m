function [tmplt,DAL,stimulus_vals,units,errstr] = BD_template(fieldname,stimulus_vals,units)
%

% AF 11/26/01
persistent   prev_playdur  prev_min_period  prev_maxlen
% We use the persistent variables to detect a change that requires some fields update.
% For example, of the play duration is changed we would like to update the gating information.
% We restict the automatic updates to allow the user the overide them.

used_devices.Llist         = 'RP1.1';
used_devices.Rlist         = 'RP2.1';
tmplt = template_definition;
if (exist('stimulus_vals','var') == 1)
   if (exist(stimulus_vals.Inloop.List_File,'file') ~= 0)
      [Llist,Rlist] = read_rotate_list_file(stimulus_vals.Inloop.List_File);
      if (~isempty(Llist))
         [data fs] = wavread(Llist{1});
      elseif (~isempty(Rlist))
         [data fs] = wavread(Rlist{1});
      end
      playback_rate = 97656.25/(round(97656.25/fs)); % 97656.25 is the rco's sampling rate
      if (~isnan(stimulus_vals.Inloop.Resample_Ratio) & ~isempty(stimulus_vals.Inloop.Resample_Ratio))
         playdur = round(stimulus_vals.Inloop.Resample_Ratio*length(data)/playback_rate*1000);	%compute file duration based on sampling rate
      else
         playdur = round(length(data)/playback_rate*1000);	%compute file duration based on sampling rate
      end
      min_period = max(1000,ceil(1.7*playdur/100)*100);
      if (~isempty(Rlist))
         min_period = min_period + ceil(0.8*playdur/100)*100;
      end
      %% In this template we have to change the Gating and Mix parameters according to 
      %% the Llist, Rlist and the playback duration of the wav files.
      %% We do this by first, updating the relevant template definitions, and second,
      %% by calling 'structdlg' in its non-interactive invisible mode, to recalculated 'stimulus_vals' fields.
      if (isempty(Llist))
         tmplt.IO_def.Mix = rmfield(tmplt.IO_def.Mix,'Llist');
      end
      if (isempty(Rlist))
         tmplt.IO_def.Mix = rmfield(tmplt.IO_def.Mix,'Rlist');
      end
      [stimulus_vals.Mix units.Mix] = structdlg(tmplt.IO_def.Mix,'',stimulus_vals.Mix,'off');
      %
      if (isequal(fieldname,'Inloop') | ~isequal(playdur,prev_playdur) | ~isequal(min_period,prev_min_period))
         tmplt.IO_def.Gating.Duration{1}  = playdur;
         tmplt.IO_def.Gating.Period{1}    = playdur*2;
         %         tmplt.IO_def.Gating.Period{1}    = ['max(' num2str(min_period) ',default_period(this.Duration))'];
         [stimulus_vals.Gating units.Gating] = structdlg(tmplt.IO_def.Gating,'',[],'off');
         prev_playdur = playdur;
         prev_min_period = min_period;
      end
      %
      if (~isequal(prev_maxlen, max(length(Llist),length(Rlist))))
         tmplt.IO_def.Inloop.Repetitions{1}  = round(100/max(length(Llist),length(Rlist)));
         new_stimulus_vals.Inloop  = structdlg(tmplt.IO_def.Inloop,'',[],'off'); % change only the repetitions!!
         stimulus_vals.Inloop.Repetitions = new_stimulus_vals.Inloop.Repetitions;
         prev_maxlen = max(length(Llist),length(Rlist));
      end
   else
      Llist = [];
      Rlist = [];
      prev_maxlen = 0;
   end
   
   Inloop.Name                         = 'DALinloop_wavfiles';
   Inloop.params.list                  = Llist;
   Inloop.params.Rlist                 = Rlist;
   Inloop.params.attens                = stimulus_vals.Inloop.Attenuation;
   Inloop.params.Rattens               = stimulus_vals.Inloop.Attenuation;
   Inloop.params.repetitions           = stimulus_vals.Inloop.Repetitions;
   Inloop.params.resample_ratio        = stimulus_vals.Inloop.Resample_Ratio;
   Inloop.params.playback_slowdown     = 1;
   

   DAL.Inloop = Inloop;
   DAL.Gating = stimulus_vals.Gating;
   DAL.Mix         = mix_params2devs(stimulus_vals.Mix,used_devices);
   DAL.short_description   = 'ROT';
   DAL.description = build_description(DAL,stimulus_vals);
   errstr = check_DAL_params(DAL,fieldname);
end

%----------------------------------------------------------------------------------------
function str = build_description(DAL,stimulus_vals)
p = DAL.Inloop.params;
[listpath,listfile] = fileparts(stimulus_vals.Inloop.List_File);
str{1} = sprintf('List ''%s'' (%d,%d files) ', listfile, length(p.list), length(p.Rlist));
if (~isempty(p.attens));
   str{1} = sprintf('%s @ %1.1f dB Attn.', str{1}, p.attens(1));
end
if (isfield(stimulus_vals.Mix,'Llist'))
   str{1} = sprintf('%s (L->%s)', str{1}, stimulus_vals.Mix.Llist);
end
if (isfield(stimulus_vals.Mix,'Rlist'))
   str{1} = sprintf('%s (R->%s)', str{1}, stimulus_vals.Mix.Rlist);
end

%----------------------------------------------------------------------------------------
function errstr = check_DAL_params(DAL,fieldname)
% Some extra error checks
errstr = '';
if (isequal(fieldname,'Inloop'))
   if (isempty(DAL.Inloop.params.attens))
      errstr = 'Attenuation is not set';
   end
   if (length(DAL.Inloop.params.attens) > 1)
      errstr = 'Only one attenuation please';
   end
end

%----------------------------------------------------------------------------------------
function tmplt = template_definition
global signals_dir
%%%%%%%%%%%%%%%%%%%%
%% Inloop Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Inloop.List_File             = { {['uigetfile(''' signals_dir 'Lists\SB\bd.m'')']} };
IO_def.Inloop.Attenuation           = {  []  'dB'    [0    120]      };
IO_def.Inloop.Repetitions            = { 20                        ''      [1    Inf]      };
IO_def.Inloop.Resample_Ratio        = { 1                        ''      [0.05    4]      };


%%%%%%%%%%%%%%%%%%%%
%% Gating Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Gating.Duration             = {410       'ms'    [20 2000]};
IO_def.Gating.Period               = {820    'ms'   [50 5000]};
%IO_def.Gating.Period               = {'default_period(this.Duration)'    'ms'   [50 5000]};
% IO_def.Gating.Rise_fall_time       = {'default_rise_time(this.Duration)' 'ms'   [0  1000]}; 

%%%%%%%%%%%%%%%%%%%%
%% Mix Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Mix.Llist        =  {'Left|Both|{Right}'};
IO_def.Mix.Rlist        =  {'Left|Both|{Right}'};

tmplt.tag               = 'BD_tmplt';
tmplt.IO_def = IO_def;
