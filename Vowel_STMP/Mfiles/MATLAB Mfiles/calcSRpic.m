function [SR_sps,lineSRs_sps]=calcSRpic(PICnum,exclude_lines,x,SRwin_sec)
% File: calcSRpic.m
% Created: 010405 (M Heinz)
%
% For: General file to calculate the SR(sps) from the "OFF-time" from a single picture
%
% Usage: [SR_sps,lineSRs_sps]=calcSRpic(PICnum,exclude_lines,x,SRwin_sec)
%    - INPUTS:
%              PICnum: picture number to load from current directory, if x (PICTURE) is not passed
%       exclude_lines: vector of lines to exlude ([] to skip)
%                   x: if passed, this contains the picture data, e.g., from x=loadpic(PICnum) [used to avoid reloading pictures if possible] ([] to skip)
%           SRwin_sec: [start_sec, end_sec] is the window used to calculate Spont-rate
%    - OUTPUTS:
%              SR_sps: estimate of Spont-rate, taken as average of all lineSRs_sps
%         lineSRs_sps: [line#, SRest_sps] SR estimated from each line


%%%% Verify in passed parameters if needed
if ~exist('exclude_lines','var')
   exclude_lines=[];
end
if ~exist('x','var')
   x=[];
end

% Only load PICnum if needed
if isempty(x)
   x=loadPic(PICnum);
end

% Verify this picture has spike data
if ~isfield(x,'spikes')
   disp(sprintf('*****\n   ERROR:  NO SPIKES in PICnum=%d\n*****',PICnum))
   beep
   SR_sps=NaN;
   lineSRs_sps=[NaN NaN];
   return;
end

%% If SR template, processs whole line (0-LINEdur)
if strcmp(deblank(x.Stimuli.short_description),'SR')
	LINEdur_sec=(x.Hardware.Trigger.StmOn+x.Hardware.Trigger.StmOff)/1000;
	SRwin_sec(1)= 0;
	SRwin_sec(2)= LINEdur_sec;
end

% Compute spontaneous rate for each line
spikeTimes = x.spikes{1};
if ~exist('SRwin_sec','var')
	
	if x.Hardware.Trigger.StmOn+x.Hardware.Trigger.StmOff-200 > 250
		OFFSET=200;
	else
		OFFSET=50;
	end
	
	SRwin_sec(1)= (x.Hardware.Trigger.StmOn + OFFSET)/1000;   % Default is 200 msec beyond end of stimulus, OR 50 if shorter
	SRwin_sec(2)= (x.Hardware.Trigger.StmOn+x.Hardware.Trigger.StmOff)/1000;
end
spont_dur_sec = diff(SRwin_sec);

lineSRs_sps=[];
for line = 1:x.Stimuli.fully_presented_lines
	if ~(sum(ismember(x.Stimuli.bad_lines,line))+sum(ismember(exclude_lines,line)))
		spikeIndices = find( (spikeTimes(:,1) == line ) );
		Nsps_spont = length (find((spikeTimes(spikeIndices,2) >= SRwin_sec(1)) & (spikeTimes(spikeIndices,2) <= SRwin_sec(2))));
		lineSRs_sps(size(lineSRs_sps,1)+1,:)=[line Nsps_spont/spont_dur_sec];
	end
end
SR_sps=mean(lineSRs_sps(:,2));
