function CD_us = findposnegCD_SCC(SCC,delays_usec)
% File:
% M. Heinz May 23, 2008
% NEW version - new
%
% OLD version: From: SCCdemo_EHrFF.m (2005)
%

% Find all peaks within 5% of true max, and pick one closest to zero delay
% posINDs=find(delays_usec>=0);
% delays_usec=delays_usec(posINDs);
% SCC=SCC(posINDs);

[pks,locs]=findpeaks(SCC);
BIGpeak_inds=find((pks-mean(SCC))>.95*(max(pks)-mean(SCC)));
BIGpeak_delays=delays_usec(locs(BIGpeak_inds));
[y,i]=min(abs(BIGpeak_delays));
CD_us=BIGpeak_delays(i);





% %%%% find CD, i.e., local max closest to zero delay
% diff_LtoR=diff(SCC);
% diff_RtoL=diff(fliplr(SCC));
% diffDelays_LtoR=delays_usec(1:end-1);
% diffDelays_RtoL=fliplr(delays_usec(2:end));
% LocalMaxDelays1=intersect(diffDelays_LtoR(find(diff_LtoR<0)),diffDelays_RtoL(find(diff_RtoL<0)));
% LocalMaxDelays=intersect(LocalMaxDelays1,delays_usec(find(SCC>=mean(SCC))));
% [y,i]=min(abs(LocalMaxDelays));
% CD_us=LocalMaxDelays(i);

return;
