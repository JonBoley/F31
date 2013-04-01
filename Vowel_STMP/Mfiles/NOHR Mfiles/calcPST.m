function PIC=calcPST(PIC,binWidth_sec)
% M.Heinz 13Sep2004.  Taken from PSTview.m
% Calculates PST histogram from PIC.
%
% Usage:PIC=calcPST(PIC)
% Input PIC: has PIC.x with stimulus parameters and spikes
% Output PIC: Stores calcs in PIC.PST, with paramaters used

if nargin<2
    % binWidth_sec = 100e-6; % matches Miller et al 1999a,b
    % binWidth_sec = 20e-6; % matches Wong et al 1998, Miller et al 1997
    binWidth_sec = 500e-6; % for Ian
end

lastBin_sec = (PIC.x.Hardware.Trigger.StmOn + PIC.x.Hardware.Trigger.StmOff) / 1000;
pst_X_sec = [0:binWidth_sec:lastBin_sec];
pst_Y_sps=hist(PIC.x.spikes{1}(:,2), pst_X_sec);
%%% Convert PST to spikes per second
pst_Y_sps=pst_Y_sps/PIC.x.Stimuli.fully_presented_lines/binWidth_sec; % Convert to sp/sec

% Store calcs
PIC.PST.pst_X_sec=pst_X_sec;
PIC.PST.pst_Y_sps=pst_Y_sps;

% Store parameters used for calcs
PIC.PST.params.binWidth_sec=binWidth_sec;

return;
