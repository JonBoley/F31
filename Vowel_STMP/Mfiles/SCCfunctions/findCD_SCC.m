function CD_us = findCD_SCC(SCC,delays_usec)
% File:
% M. Heinz May 23, 2008
% NEW version - new
%
% OLD version: From: SCCdemo_EHrFF.m (2005)
%

% Find all peaks within 5% of true max, and pick one closest to zero delay
[pks,locs]=findpeaks(SCC);
BIGpeak_inds=find((pks-mean(SCC))>.95*(max(pks)-mean(SCC)));
BIGpeak_delays=delays_usec(locs(BIGpeak_inds));
[y,i]=min(abs(BIGpeak_delays));
CD_us=BIGpeak_delays(i);


    %%% Automatically select delays %%%
%     CDtfs_usec(k)=findposnegCD_SCC(SACSCCfunctions{end}.DIFCOR_AB,SACSCCfunctions{end}.delays_usec);
    % or...
    %%% Manually select delays %%%
%     scrsz = get(0,'ScreenSize');
%     figure('Position',[1 1 scrsz(3) scrsz(4)]);
%     plot(SACSCCfunctions{end}.delays_usec((end/4:3*end/4)),SACSCCfunctions{end}.SCC_AB_avg(end/4:3*end/4),'b'); hold on;
%     plot(SACSCCfunctions{end}.delays_usec((end/4:3*end/4)),SACSCCfunctions{end}.XpCC_AB_avg(end/4:3*end/4),'b');
% %     if exist('x'), plot(x,1,'.g'); end
%     plot(SACSCCfunctions{end}.delays_usec((end/4:3*end/4)),(SACSCCfunctions{end}.DIFCOR_AB(end/4:3*end/4)),'rx-');  hold off;
%     title(sprintf('Pick a peak\n%s\n%1.1f octaves',DATAfilename,octdiffs(k)));
%     [x,y] = ginput(1);
%     close(gcf);
%     CDtfs_usec(k)=x;%SACSCCfunctions{end}.delays_usec((round(x)-1)*2-1);


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
