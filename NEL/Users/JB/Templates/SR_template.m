function [tmplt,DAL,stimulus_vals,units,errstr] = PST_template(fieldname,stimulus_vals,units)
%
% Modified by M.Heinz 12Nov2003, from nel_rate_level_template.m

% Adapted from "nel_pst_template.m" by GE, 29Mar2002.

used_devices.Tone       = 'RP1.1';
tmplt = template_definition(fieldname);
if (exist('stimulus_vals','var') == 1)
%   Inloop.Name                               = 'DALinloop_general_TN_delay';  % modified by GE (07Jan2003)  %Removed MH(12Nov2003)
   Inloop.Name                               = 'DALinloop_general_TN';
   Inloop.params.repetitions                 = stimulus_vals.Inloop.Repetitions;
%    Inloop.params.main.source                 = 'Tone';
   Inloop.params.main.source                 = 'None';
   Inloop.params.main.tone.freq              = 100;
   Inloop.params.main.tone.bw                = 0;
   Inloop.params.main.noise.low_cutoff       = 0;
   Inloop.params.main.noise.high_cutoff      = 0;
   Inloop.params.main.attens                 = 120;    % shut off to get SR
   Inloop.params.secondary.source            = 'None';
   Inloop.params.secondary.tone.freq         = 0;
   Inloop.params.secondary.noise.low_cutoff  = 0;
   Inloop.params.secondary.noise.high_cutoff = 0;
   Inloop.params.secondary.noise.gating      = '';
   Inloop.params.secondary.noise.adaptation  = 0;
   Inloop.params.secondary.atten             = 120;
   Inloop.params.rise_fall                   = 20;
%    stimulus_vals.Gating.Period = 1000;
   stimulus_vals.Gating.Duration = 100;

   DAL.Inloop = Inloop;
   DAL.Gating = stimulus_vals.Gating;
   DAL.short_description   = 'SR';
   DAL.endLinePlotParams                  = nel_plot_pst_params(DAL.Gating.Period/1000, DAL.Gating.Duration/1000);  % GE 04Nov2003.

   % [stimulus_vals.Mix units.Mix] = structdlg(tmplt.IO_def.Mix,'',stimulus_vals.Mix,'off');
   DAL.Mix         = mix_params2devs(stimulus_vals.Mix,used_devices);

   DAL.description = build_description(DAL,stimulus_vals);
   errstr = check_DAL_params(DAL,fieldname);
end

%----------------------------------------------------------------------------------------
function str = build_description(DAL,stimulus_vals)
p = DAL.Inloop.params;
str{1} = sprintf('PST %d Repetitions', p.repetitions);
str{1} = sprintf('%s %1.2f kHz Tone', str{1}, p.main.tone.freq/1000);
str{1} = sprintf('%s @ %1.1f dB Attn.', str{1}, p.main.attens);
str{1} = sprintf('%s (%s)', str{1}, stimulus_vals.Mix.Tone);

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
end

%----------------------------------------------------------------------------------------
function tmplt = template_definition(fieldname)
persistent prev_unit_bf prev_unit_thresh
%%%%%%%%%%%%%%%%%%%%
%% Inloop Section 
%%%%%%%%%%%%%%%%%%%%

% IO_def.Inloop.Frequency    =  {0.01   'kHz'      [0.04  50]   0  0};
% IO_def.Inloop.Attenuation  =  {120 'dB'       [0    120]   0  0};

IO_def.Inloop.Repetitions  =  {20   ''       [1 600]}; 

% if (isequal(fieldname,'Inloop'))
%    if (~isequal(current_unit_bf, prev_unit_bf))
%       IO_def.Inloop.Frequency{5}            = 1; % ignore dflt. Always recalculate.
%       prev_unit_bf = current_unit_bf;
%    end
%    if (~isequal(current_unit_thresh, prev_unit_thresh))
%       IO_def.Inloop.Attenuation{5}     = 1;
%       prev_unit_thresh = current_unit_thresh;
%    end
% end

%%%%%%%%%%%%%%%%%%%%
%% Gating Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Gating.Period           = {1000    'ms'   [50 120000]};
% IO_def.Gating.Rise_fall_time   = {20 'ms'   [0  1000]}; 
IO_def.Gating.Duration         = {100       'ms'    [20 60000]};

%%%%%%%%%%%%%%%%%%%%
%% Mix Section 
%%%%%%%%%%%%%%%%%%%%
IO_def.Mix.Tone      =  {'{Left}|Both|Right'};

tmplt.tag         = 'SR_tmplt';
tmplt.IO_def = IO_def;
