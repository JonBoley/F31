function [varargout] = sptr;
% PSTH - makes a PST histogram of files entered interactively. Handles all
% lab standard data types and includes code to process Hafter TDT click data.
% [absc, pstd] = psth returns the x-axis and the PSTs, as displayed, in spikes
% per second as column vectors (one abscissa, in ms, for all PSTs).

global PSTDATXYZ INDEXXYZ
str = input('Enter files, binwidth, filter (a4-6,7[b0.2][f5][h]): ','s');
[absc] = sptrain(str);
if nargout>=1
	varargout(1) = {absc};
end
if nargout>=2
	varargout(2) = {PSTDATXYZ};
end

return

clear PSTDATXYZ
clear INDEXXYZ
clear str
clear absc
