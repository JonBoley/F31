function [CD,PEAK]=calcCD(Correl,Delays,F0per_us)
% File: calcCD
% 22Jun2009: R Sono - took out of UnitLook_EHIN_CoincDet2_MH2
%
%
% Computes Characteristic Delay (CD) of the correllogram Correl. PEAK is the
% value of Correl at CD.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Finding CD is TOUGH!!!! %%%%                                                %
%                                                                                %
% Need a robust way to do it                                                     %
% 1) Find Local Maxima within 1 period, based on 3-pt smoothed NSCC              %
% 2) Restrict to only those within 15% of max                                    %
% 3) Take local max closest to 0                                                 %
% *** THIS SOMETIMES CHOOSES CD close to 0, when the NSCC is fairly periodic     %
%     BUT, these are just tough ones to get right without a lot of "knowledge"   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if~exist('F0per_us','var')
    F0per_us=8.5e5;%I don't know if this is a good default. 120Hz gets you 8.33e5.
end

delaysTEMP_usec=Delays(find((Delays>=-F0per_us)&(Delays<=F0per_us)));
corrTEMP=Correl(find((Delays>=-F0per_us)&(Delays<=F0per_us)));

%%%%%%%%%%% LOCAL MAXIMA
TFiltWidth=3;
%%%% find local max within 1 period, based on smoothed SCCs
diff_LtoR=diff(trifilt(corrTEMP,TFiltWidth));
diff_RtoL=diff(fliplr(trifilt(corrTEMP,TFiltWidth)));
diffDelays_LtoR=delaysTEMP_usec(1:end-1);
diffDelays_RtoL=fliplr(delaysTEMP_usec(2:end));
LocalMaxDelays1=intersect(diffDelays_LtoR(find(diff_LtoR<0)),diffDelays_RtoL(find(diff_RtoL<0)));
LocalMaxDelays=intersect(LocalMaxDelays1,delaysTEMP_usec(find(corrTEMP>=1)));
corrLocalMax=zeros(size(LocalMaxDelays));
for i=1:length(LocalMaxDelays)
    corrLocalMax(i)=corrTEMP(find(delaysTEMP_usec==LocalMaxDelays(i)));
end

%% Restrict to local max within 15% of peak
PercCRIT=0.15;
RESTRICTinds=find(corrLocalMax>max(corrLocalMax)*(1-PercCRIT));
LocalMaxDelaysRESTRICT=LocalMaxDelays(RESTRICTinds);
corrLocalMaxRESTRICT=corrLocalMax(RESTRICTinds);

if ~isempty(LocalMaxDelaysRESTRICT)
    [y,i]=min(abs(LocalMaxDelaysRESTRICT));
    CD=LocalMaxDelaysRESTRICT(i);%NSCC_CDs_usec{ROWind,ATTind}{SCCind}
    PEAK=corrLocalMaxRESTRICT(i);  %NSCC_peaks{ROWind,ATTind}{SCCind} % This is not actually the peak, but the NSCC at the CD
else
    CD=NaN;
    PEAK=NaN;
end
return;