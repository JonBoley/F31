function varargout = DALinloop_NI_SCC_wavfiles2(varargin)
% Adapted by MH 07July2004 from 'DALinloop_NI_wavefiles' (GE 7/26/02)
%   adds the ability to invert the waveform (based on static_bi.Condition.InvertPolarity)
%
% Adapted by GE 26Jul2002 from 'DALinloop_wavfiles' (AF 9/22/01).
%  Further modifications: GE 03oct2003, 04oct2003, 06oct2003.
%
%   The input and output arguments are structure with the following fields
%   varargin{1} (common)   : index, left, right
%   varargin{2} (specific) : see list below
%
%   varargout{1} (stim_info)  : attens_devices 
%   varargout{2} (block_info) : nstim, nlines, stm_lst,  list,attens,Rlist,Rattens
%   varargout{3} (plot_info) : var_name, var_unit, var_vals, var_frmt, XYprops
%
%  specific->
%               list: []
%             attens: []
%     resample_ratio: []
%     playback_slowdown: []
%              Rlist: []
%            Rattens: []
%        repetitions: []


global static_bi static_pi; % Static variable
global root_dir % Used to find correct control 'rco' file.
persistent Chan0data Chan1data  % To avoid re-reading file, when waveform is 
                                %  same as for previous line.

rc = 1;  % set return code to "good" at start of function.

% Initialization (if no input arguments) or setting of input arguments into
%  local variables 'common' and 'specific'.
if (nargin == 0)    % "clean-up" only (at end of data_acquisition_loop)
   static_bi = [];
   static_pi = [];
   clear static_bi static_pi;
   [statCodes retVals] = d2a(2);  % Terminate/Clean-up NI 6052 board (see file "d2a.c")
   rc_d2a = check_error_codes(statCodes, retVals, 37);
   return;
elseif (nargin >=1)
   common = varargin{1};
end
if (nargin == 2)
   specific = varargin{2};
end

% if common.index==1                % ge debug
%    sprintf('\n\n***************')
% end
% disp(sprintf('Index=%d',common.index)) 


% Choose action based on call-mode ('common.index'):
if (common.index == 0)  % initialization (Pre-Loop Call, do nothing)
   return;  % "clean-up" is done -- return.
end
if (common.index == 1)  % first picture
%    [statCodes retVals] = d2a(0);     % Initialize NI 6052 board. (see "d2a.c").
%    rc_d2a = check_error_codes(statCodes, retVals, 51);
   
    %%%% Initialize and set specific values to the static block info structure
   static_pi = default_inloop_plot_info;
   static_bi = specific;
   if (~isfield(specific,'repetitions'))
      specific.repetitions = 1;
   end

   % Determine number of distinct output waveforms that will be delivered during
   %  the picture (with different attenuations of the same wavefile counting as
   %  distinct outputs):
   
   % 12Apr2005 MHeinz
   nstim = max(length(static_bi.list),length(static_bi.attens));
   
   % Check to make sure that # of attenuations matches number of files, if using list.
   if (min(length(static_bi.list),length(static_bi.attens)) ~= 1 & length(static_bi.list) ~= length(static_bi.attens))
      error(sprintf('Inconsistent number of attenuations (%d) and file names (%d)', ...
         length(static_bi.attens), length(static_bi.list)));
   end
   
   if (nstim == 1 & static_bi.repetitions == 1)
      nlines = 100;
      static_pi.var_name = 'Presentation number';
   else
      nlines = nstim * static_bi.repetitions;
   end
   if (length(static_bi.attens) > 1)
      if (static_bi.repetitions == 1)
         static_pi.var_name = 'Attenuation';
         static_pi.var_unit = 'dB';
         static_pi.var_frmt = '%d';

         if length(unique(static_bi.attens))>1
            static_pi.var_vals = static_bi.attens;
            static_pi.XYprops.Lim  = [min(static_bi.attens)-1 max(static_bi.attens)+1];
         else  % This is a RLF in the 2nd channel
            static_pi.var_vals = static_bi.Rattens;
            static_pi.XYprops.Lim  = [min(static_bi.Rattens)-1 max(static_bi.Rattens)+1];
         end

         static_pi.XYprops.Tick = [0:10:120];
         static_pi.XYprops.Dir  = 'reverse';
      else
         var_labels = cell(1,nlines);
         counter = 1;
         for ii = 1:static_bi.repetitions
            for jj = 1:length(static_bi.attens)
               var_labels{counter} = sprintf('%d (Rep. #%d) ', static_bi.attens(jj), ii);
               counter = counter+1;
            end
         end
         static_pi.var_name = 'Attenuation';
         static_pi.var_unit = 'dB';
         static_pi.var_frmt = '%s';
         shift              = repmat([0:static_bi.repetitions-1]/(static_bi.repetitions*2),length(static_bi.attens),1);
         static_pi.var_vals     = repmat(static_bi.attens,1,static_bi.repetitions) + shift(:)';
         static_pi.var_labels   = var_labels;
         static_pi.XYprops.Lim  = [min(static_bi.attens)-1 max(static_bi.attens)+1];
         static_pi.XYprops.Tick = [0:10:120];
         static_pi.XYprops.Dir  = 'reverse';
      end
   else
      if (nstim > 1)
         static_pi.var_name = 'File';
         static_pi.var_unit = '';
         static_pi.var_frmt = '%s';
         var_labels = cell(1,nlines);
         counter = 1;
         for ii = 1:specific.repetitions
            for jj = 1:nstim
               [dummypath fname] = fileparts(static_bi.list{jj});
               var_labels{counter} = fname; 
               if (~isempty(static_bi.Rlist))
                  [dummypath fname] = fileparts(static_bi.Rlist{jj});
                  var_labels{counter} = [var_labels{counter} ' + ' fname]; 
               end
               if (ii > 1)
                  var_labels{counter} = [var_labels{counter} ' (Repetition #' int2str(ii) ')']; 
               end
               counter = counter+1;
            end   
         end
         static_pi.var_labels   = var_labels;
         static_pi.var_vals     = 1:nlines;
         static_pi.XYprops.Lim  = [0 nlines+1];
         % static_pi.XYprops.Tick = [0:10:120];
         static_pi.XYprops.Dir  = 'normal';
      else
         static_pi.var_name = 'Repetition #';
      end
   end
   
   % Load 'control.rco' to RP (for triggering and switchbox control):
   rconame = [root_dir 'stimulate\object\control.rco'];
   if (isempty(static_bi.Rlist))
      rc = RPload_rco(rconame) & (rc==1);   %Creates click: MH 10/20/03
       % use L3 as source for Chan0data:
      dev_description = nel_devices_vector('L3','list');
   else
      rc = RPload_rco({rconame, rconame}) & (rc==1);
       % Use L3 and R3, respectively, as sources for Chan0data and Chan1data:
      dev_description = nel_devices_vector({'L3','R3'},{'list','Rlist'});
   end
   if (rc == 0)
      return;
   end

   short_description = 'wav_NIboard';
   if (isfield(common,'short_description'))
      short_description = common.short_description;
   end
   if (isfield(common,'description'))
      description = common.description;
   end
   %%%%% These fields should be set in ANY inloop function %%%%%%
   static_bi.nstim           = nstim;
   static_bi.nlines          = nlines;
   static_bi.description     = description;
   static_bi.short_description  = short_description;
   static_bi.dev_description = dev_description;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else  %MH: debug
   %% If not the first time through, reset NI-board
   [statCodes retVals] = d2a(2);  % Terminate/Clean-up NI 6052 board (see file "d2a.c")
   rc_d2a = check_error_codes(statCodes, retVals, 37);
end

% If first line or using list file, need to load waveform from file:
if ((common.index == 1) | (length(static_bi.list) > 1))

   Chan0data = []; % initialize and empty the data vector for channel 0.
   Chan1data = []; % initialize and empty the data vector for channel 1. (14Apr2005 MH: need this because of persistent Chan1data, stays around if you've run Rlist previously

   ind = mod(common.index-1,length(static_bi.list))+1;

   % Load the waveform from the file and place in 'Chan0data':
   [data,sr,rc] = nel_wavread(static_bi.list{min(ind,end)});
   if (rc == 0)
      nelerror(['Can''t read wavfile ''' static_bi.list{min(ind,end)} '''']);
   else
      Chan0data = data;
   end
   
   % Invert Stimulus Polarity if desired (MH 07July2004)
%    if strcmp(static_bi.Condition.InvertPolarity,'yes')
%       Chan0data=-Chan0data;
%    end
   
   % Repeat or truncate waveform to fit requested stimulus duration:
   Chan0data = refit_waveform(Chan0data);   
   
   % Window waveform using linear rise/fall:
   Chan0data = window_waveform(Chan0data);
   
   % Rescale data to provide peak-to-peak voltage output that matches TDT implementation:
   Chan0data = rescale_waveform(Chan0data);
   
   
   if (length(static_bi.list) > 1)
      stim_info.file{1}                    = static_bi.list{min(ind,end)};
      stim_info.playback_sampling_rate(1) = static_bi.updateRate_Hz;
      %stim_info.playback_sampling_rate(1)  = sr;
   else                                           
      static_bi.playback_sampling_rate(1)  = static_bi.updateRate_Hz;
      %static_bi.playback_sampling_rate(1)  = sr;
   end
   
end

%A right list file exists AND (either first time through OR more than one element in list)
if (((common.index == 1) | (length(static_bi.Rlist) > 1)) & (length(static_bi.Rlist) >0))
   Chan1data = []; % empty the data vector.

   ind = mod(common.index-1,length(static_bi.Rlist))+1;
   
   % Load the waveform from the file and place in 'Chan1data'.
   [data,sr,rc] = nel_wavread(static_bi.Rlist{min(ind,end)});
   if (rc == 0)
      nelerror(['Can''t read wavfile ''' static_bi.Rlist{min(ind,end)} '''']);
   else
      Chan1data = data;
   end
   
   % Invert Stimulus Polarity if desired (MH 07July2004)
%    if strcmp(static_bi.Condition.InvertPolarity,'yes')
%       Chan1data=-Chan1data;
%    end
   
   % Repeat or truncate waveform to fit requested stimulus duration:
   Chan1data = refit_waveform(Chan1data);   
   
   % Window waveform using linear rise/fall:
   Chan1data = window_waveform(Chan1data);
   
   % Rescale data to provide peak-to-peak voltage output that matches TDT implementation:
   Chan1data = rescale_waveform(Chan1data);

   if (length(static_bi.Rlist) > 1)
      stim_info.file{2}                    = static_bi.Rlist{min(ind,end)};
      stim_info.playback_sampling_rate(2) = static_bi.updateRate_Hz;
   else                                           
      static_bi.playback_sampling_rate(2)  = static_bi.updateRate_Hz;
   end

end


% % GE 04Nov2003: Error checking here using "NI6052UsableRate_Hz" to be sure that valid
% %   board update rate has been specified.  Should be checked in template, but in case it's not...
% if (static_bi.updateRate_Hz ~= NI6052UsableRate_Hz(static_bi.updateRate_Hz))
%    nelerror(['ERROR in <DALinloop_NI_wavfiles>: invalid update rate for NI board!!!!!']);
% end

% Load the samples data to the NI board and enable trigger (has to be done anew for every line):
if (~isempty(Chan1data))    % load both channels.
   samplesMatrix = [Chan0data Chan1data];
elseif (~isempty(Chan0data))    % load channel 0 only.
   samplesMatrix = [Chan0data];
end
% 
% MH/GE debug: return actualUpdateRate from d2a to check
%   [statCodes retVals] = d2a(1, static_bi.updateRate_Hz, samplesMatrix);
[statCodes retVals actualUpdateRate_Hz] = d2a(1, static_bi.updateRate_Hz, samplesMatrix);
% GE/MH 05Nov2003: Error checking here using actualUpdateRate returned from d2a.c to be sure that valid
%   board update rate was specified.  Should be checked in template, but in case it's not...
if actualUpdateRate_Hz~=static_bi.updateRate_Hz
   nelwarn(['WARNING in <DALinloop_NI_wavfiles>: NI board / Matlab update rate mismatch!']);
end
rc_d2a = check_error_codes(statCodes, retVals, 231);



% Set attenuations and devices
ind = mod(common.index-1,length(static_bi.attens))+1;
atten = static_bi.attens(min(ind,end));
left_dev =  atten * nel_devices_vector('L3');  % Using L3 as device for Chan0data.
if (isempty(static_bi.Rlist))
   right_dev = nel_devices_vector([]);
else
   if (~isempty(static_bi.Rattens))
      ind = mod(common.index-1,length(static_bi.Rattens))+1;
      Ratten = static_bi.Rattens(min(ind,end));
   else
      Ratten = atten;
   end
   right_dev = Ratten * nel_devices_vector('R3');  % Using R3 as device for Chan1data.
end
devs     = max([left_dev right_dev],[],2);
stim_info.attens_devices = [devs devs];     


% Set output arguments:
if (nargout >=1)
   varargout{1} = stim_info;
end
if (nargout >=2)
   varargout{2} = static_bi;
end
if (nargout >=3)
   varargout{3} = static_pi;
end


%##########################################################################################
function out_wave = refit_waveform(in_wave)
   global NelData static_bi
   stimOn_sec = NelData.DAL.Gating.Duration / 1000;
   ur_Hz = static_bi.updateRate_Hz;
   required_nSamples = floor(ur_Hz*stimOn_sec);
   current_nSamples = size(in_wave, 1);
   if (current_nSamples > required_nSamples) % Waveform needs to be truncated.
      out_wave = in_wave(1:required_nSamples);
      
      beep off
      nelwarn(['In function ''DALinloop_NI_wavfiles'': Input waveform has been truncated ' ...
                 'to fit requested duration.']);
           beep on
   elseif (current_nSamples < required_nSamples) % Waveform needs to be "repeated".
      nRepeats = ceil(required_nSamples / current_nSamples);
      out_wave = repmat(in_wave, nRepeats, 1);
      out_wave = out_wave(1:required_nSamples); % truncate any extra of the final repeat.
      nelwarn(['In function ''DALinloop_NI_wavfiles'': Input waveform has been repeated ' ...
                 'to fill requested duration.']);
   else
      out_wave = in_wave;
   end


%##########################################################################################
function out_wave = window_waveform(in_wave)
   global NelData 
   stimOn_sec = NelData.DAL.Gating.Duration / 1000;
   ur_Hz = NelData.DAL.Inloop.params.updateRate_Hz;
   rf_time_sec = NelData.DAL.Gating.Rise_fall_time / 1000;
   ramp_nSamples = ceil(ur_Hz * rf_time_sec);
   end_sample = size(in_wave, 1);
   windowing_vector = zeros(size(in_wave));
   windowing_vector(1:ramp_nSamples) = [0:(ramp_nSamples-1)]/ramp_nSamples;
   windowing_vector((ramp_nSamples+1):(end_sample-ramp_nSamples)) = 1;
   windowing_vector(((end_sample-ramp_nSamples)+1):end) = [(ramp_nSamples-1):-1:0]/ramp_nSamples;
   out_wave = windowing_vector .* in_wave;

%##########################################################################################
function out_wave = rescale_waveform(in_wave)
    scale_factor = 5;  % A .wav file with range [-1,1] is scaled to an output of 
                       %  +/-5V at the TDT analog output, before attenuation.
    % For the NI board implementation, multiply the .wav file data by this factor.
    %   The implementation of the NI board (d2a.c) assumes that the input data are
    %   already scaled to the exact voltage level.
    out_wave = scale_factor * in_wave;
   
   
%##########################################################################################
function rc_d2a = check_error_codes(statCodes, retVals, lineNum)
   rc_d2a = 0;
   if ( (~isempty(find(statCodes~=0))) | (~isempty(find(retVals~=0))) )
      rc_d2a = -1;
      disp(sprintf('\tDALInloop Location: %d',lineNum))
      disp([sprintf('\t\t') 'StatCodes: ' sprintf('%d ',statCodes')])
      disp([sprintf('\t\t') 'Return Vals: ' sprintf('%d ',retVals')])
      nelwarn(['non-zero status/return codes in ''DALinloop_NI_wavfiles'' ']);
   end
   