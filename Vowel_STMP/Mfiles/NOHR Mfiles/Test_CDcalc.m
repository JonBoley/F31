% File: Test_CDcalc.m
% 
% 11Feb2005 M.Heinz
% used to test new, more robust way to calculate CD



%% SAVE data elsewhere so we can keep originals
F0per_us_XX=F0per_us;
delaysTEMP_usec_XX=delaysTEMP_usec;
nsccTEMP_XX=nsccTEMP;



[y_XX,i_XX]=max(nsccTEMP_XX);
NSCC_CDs_usec_XX=delaysTEMP_usec_XX(i_XX);
NSCC_peaks_XX=y_XX;




%%%%%%%%%%% LOCAL MAXIMA
TFiltWidth=5;
%%%% find local max within 1 period, based on smoothed SCCs
diff_LtoR=diff(trifilt(nsccTEMP_XX,TFiltWidth));
diff_RtoL=diff(fliplr(trifilt(nsccTEMP_XX,TFiltWidth)));
diffDelays_LtoR=delaysTEMP_usec_XX(1:end-1);
diffDelays_RtoL=fliplr(delaysTEMP_usec_XX(2:end));
LocalMaxDelays1=intersect(diffDelays_LtoR(find(diff_LtoR<0)),diffDelays_RtoL(find(diff_RtoL<0)));
LocalMaxDelays=intersect(LocalMaxDelays1,delaysTEMP_usec_XX(find(nsccTEMP_XX>=1)));

nsccLocalMax=zeros(size(LocalMaxDelays));
for i=1:length(LocalMaxDelays)
   nsccLocalMax(i)=nsccTEMP_XX(find(delaysTEMP_usec_XX==LocalMaxDelays(i)));
end

%% Restrict to local max within 10% of peak
PercCRIT=0.15;
RESTRICTinds=find(nsccLocalMax>max(nsccLocalMax)*(1-PercCRIT));
LocalMaxDelaysRESTRICT=LocalMaxDelays(RESTRICTinds);
nsccLocalMaxRESTRICT=nsccLocalMax(RESTRICTinds);

[y_XX,i_XX]=min(abs(LocalMaxDelaysRESTRICT));
CD_us_XX=LocalMaxDelaysRESTRICT(i_XX);
NSCC_CD_XX=nsccLocalMaxRESTRICT(i_XX);





figure
set(gcf,'pos',[115         540        1281         420])


plot(delaysTEMP_usec_XX/1000,nsccTEMP_XX)
hold on
plot([min(delaysTEMP_usec_XX) max(delaysTEMP_usec_XX)]/1000,ones(1,2),'k-')
plot(NSCC_CDs_usec_XX/1000,NSCC_peaks_XX,'go')
plot(CD_us_XX/1000,NSCC_CD_XX,'bs')  % NEW CALC
plot(LocalMaxDelays/1000,nsccLocalMax,'rx')


xlabel('Delay (ms)')





clear F0per_us_XX delaysTEMP_usec_XX nsccTEMP_XX
