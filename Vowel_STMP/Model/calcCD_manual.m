function [CD,PEAK]=calcCD_manual(Correl,Delays,F0per_us,NSCC_CDs_usec,BFs_kHz,SCCind)
% File: calcCD
% 22Jun2009: R Sono - took out of UnitLook_EHIN_CoincDet2_MH2
% 01Oct2009: J Boley - added manual selection of CD
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Manually select delays %%%
scrsz = get(0,'ScreenSize');
figure('Position',[1 1 scrsz(3) scrsz(4)]);
subplot(2,1,1), plot(delaysTEMP_usec,corrTEMP,'b');
hold on; 

% % Plot Hilbert Env
% spec_temp = abs(fft(corrTEMP));
% peak_index = 1+find(spec_temp(2:end)==max(spec_temp(2:end)),1);
% % DF = ...; % dominant frequency
% hilbertenv = abs(hilbert(corrTEMP));
% n = round(1e6 / DF);
% win = hanning(n)/sum(hanning(n)); % LPF
% SCC_env = filter(win,1,hilbertenv);
% plot(delaysTEMP_usec,abs(hilbert(corrTEMP)),'g:');

subplot(2,1,2), hold on;
try
    subplot(2,1,1), plot(CD,corrTEMP(find(delaysTEMP_usec>=CD,1)),'rx-','MarkerSize',12,'LineWidth',2);
    octdiff = log2(BFs_kHz{SCCind}(1)/BFs_kHz{SCCind}(2));
    subplot(2,1,2), plot(octdiff,CD,'rx-','MarkerSize',12,'LineWidth',2);
catch
    subplot(2,1,1), text(0,corrTEMP(find(delaysTEMP_usec>=0,1)),'Could not find a peak!');
end
for i=1:size(NSCC_CDs_usec,2)
    if ~isempty(NSCC_CDs_usec{i})
        if ~isnan(CD)
            subplot(2,1,1), plot(NSCC_CDs_usec{i},corrTEMP(find(delaysTEMP_usec>=CD,1)),'b.');
        else
            subplot(2,1,1), plot(NSCC_CDs_usec{i},corrTEMP(find(corrTEMP==max(corrTEMP),1)),'b.');
        end
        octdiff = log2(BFs_kHz{i}(1)/BFs_kHz{i}(2));
        subplot(2,1,2), plot(octdiff,NSCC_CDs_usec{i},'b.');
    end
end

subplot(2,1,2), hold off;
title('Previous characteristic delays');
subplot(2,1,1), hold off;
octdiff = log2(BFs_kHz{SCCind}(1)/BFs_kHz{SCCind}(2));
title(sprintf('Pick a peak (or press Enter to accept)\nBF Difference = %1.2f octaves\n(BF1 = %1.2f; BF2 = %1.2f; Formant = %1.2f)',octdiff,BFs_kHz{SCCind}(1),BFs_kHz{SCCind}(2),1e6/F0per_us));

[x,y] = ginput(1);
close(gcf);
if ~isempty(x),
    CD=x;
    PEAK=corrTEMP(find(delaysTEMP_usec>=CD,1));
    if ((x<delaysTEMP_usec(1)) || (x>delaysTEMP_usec(end)))
        CD=NaN;
        PEAK=NaN;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


return;