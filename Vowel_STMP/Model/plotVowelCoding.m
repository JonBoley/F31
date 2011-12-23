function plotVowelCoding()
% plotVowelCoding()
% This function plots (for the first 3 formants):
%  - period histogram re BF
%  - strength of rate coding re BF
%  - strength of ALSR coding re BF
%  - cross-CF correlation
%  - cross-CF characteristic delay
%
% Input variables:
%    FeatINDs
%    BFs_kHz
%    Nattens_dB
%    PERhistXs_sec
%    PERhists_Smoothed
%    Rates
%    ALSRs
%    Synchs
%    data2plot

% Variables to set up: (also F2 & F3)
%    unit.EHvN_reBF_simFF.F1{1,1}.BFs_kHz(freqIND)
%    unit.EHvN_reBF_simFF.F1{1,1}.levels_dBSPL
%    unit.EHvN_reBF_simFF.F1{1,1}.Nattens_dB
%    unit.EHvN_reBF_simFF.F1{1,1}.picNums
%    unit.EHvN_reBF_simFF.F1{1,1}.excludeLines
%    unit.EHvN_reBF_simFF.F1{1,1}.rate(attenIND,freqIND)
%    unit.EHvN_reBF_simFF.F1{1,1}.synch{attenIND,freqIND}
%    unit.EHvN_reBF_simFF.F1{1,1}.phase{attenIND,freqIND}
%    unit.EHvN_reBF_simFF.F1{1,1}.RaySig{attenIND,freqIND}
%    unit.EHvN_reBF_simFF.F1{1,1}.TimeFact{attenIND,freqIND}
%    unit.EHvN_reBF_simFF.F1{1,1}.FeatureFreqs_Hz{freqIND}
%    unit.EHvN_reBF_simFF.F1{1,1}.FeatureLevels_dB


% vowel variables:
FeaturesText={'F0','T0','F1','T1','F2','T2','F3','T3','TN'};
FeatINDs=[3,5,7]; % just look at F1, F2, F3

% plotting variables:
ATTENcolors={'b','r','g','k','c','y'};
FeatureColors={'r','g','b'};

NUMcols=5;
NUMrows=3;
figure(501); clf
set(gcf,'Name','Comparison of Neural Coding Strategies');
set(gcf,'units','norm','pos',[0.2234    0.7119    0.4297    0.2344]);%,'Resize','off')
ROWind=0;
ALLlevelsTriFiltTONE=9;
ALLlevelsTriFilt=3;

ATTEN=max(Nattens_dB);

% Hardcode for now
HarmonicsIND=1;
PolarityIND=1;
    
for FeatIND=FeatINDs
    ROWind=ROWind+1;
    eval(['yTEMP=unit.EHvN_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
    if ~isempty(yTEMP)
        %%%% EH_reBF plots
        ATTind=find(yTEMP.Nattens_dB==ATTEN);
        
        %%%% Spatio-Temporal Plots
        PLOTnum=(ROWind-1)*NUMcols+1;
        eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
        LEGtext='';
        for BFind=1:length(BFs_kHz{ROWind,ATTind})
            if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,ATTind}/unit.Info.BF_kHz))<BFoctCRIT))
                LINEwidth=2;
            else
                LINEwidth=.5;
            end
            % This normalization plots each signal the same size on a log scale
            if ~isempty(PERhistXs_sec{ROWind,ATTind}{BFind})
                NormFact=(10^(PERhistGAIN*PERhists_logCHwidth)-1)*BFs_kHz{ROWind,ATTind}(BFind)/PERhistsMAX;
                semilogy(PERhistXs_sec{ROWind,ATTind}{BFind}*1000, ...
                    trifilt(PERhists_Smoothed{ROWind,ATTind}{BFind},ALLlevelsTriFilt)*NormFact+BFs_kHz{ROWind,ATTind}(BFind), ...
                    'LineWidth',LINEwidth,'Color',ATTENcolors{ATTind})
                hold on
                if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,ATTind}/unit.Info.BF_kHz))<BFoctCRIT))
                    LEGtext{length(LEGtext)+1}=sprintf('%.f dB',Nattens_dB(ATTind)+dBAtt_2_SNR);
                end
                for ATTind2=fliplr(find(Nattens_dB~=max(Nattens_dB)))
                    if ~isempty(PERhistXs_sec{ROWind,ATTind2}{BFind})
                        semilogy(PERhistXs_sec{ROWind,ATTind2}{BFind}*1000, ...
                            trifilt(PERhists_Smoothed{ROWind,ATTind2}{BFind},ALLlevelsTriFilt)*NormFact+BFs_kHz{ROWind,ATTind2}(BFind), ...
                            'LineWidth',LINEwidth,'Color',ATTENcolors{ATTind2})
                        if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,ATTind}/unit.Info.BF_kHz))<BFoctCRIT))
                            LEGtext{length(LEGtext)+1}=sprintf('%.f dB',Nattens_dB(ATTind2)+dBAtt_2_SNR);
                        end
                    end
                end
                if strcmp(FeaturesText{FeatIND},'F1')
                    hleg1000=legend(LEGtext,1);
                    set(hleg1000,'FontSize',10)
                    set(hleg1000,'pos',[0.356    0.8    0.07    0.17])
                end
            end
        end
        xlabel('Time (ms)')
        ylabel('Effective Best Frequency (kHz)')
        if ROWind==1
            title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n%s @ %.f dB SPL,   Harm: %d, Polarity: %d', ...
                ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10,FeaturesText{FeatIND}, ...
                yTEMP.levels_dBSPL,HarmonicsIND,PolarityIND),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
        else
            title(sprintf('%s @ %.f dB SPL,   Harm: %d, Polarity: %d',FeaturesText{FeatIND}, ...
                yTEMP.levels_dBSPL,HarmonicsIND,PolarityIND),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
        end
        xlim(XLIMITS_perhist)
        PLOThand=eval(['h' num2str(PLOTnum)]);
        set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
        ylim(YLIMITS)  % Same Ylimits for all plots
        %%%%%%%%%%%%%%%%%%%%%
        % Plot lines at all features
        for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
            if (yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000>=YLIMITS(1))&&(yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000<=YLIMITS(2))
                semilogy(XLIMITS_perhist,yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                text(XLIMITS_perhist(2)*1.005,yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000, ...
                    sprintf('%s (%.1f)',FeaturesText{FeatINDPlot},yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000), ...
                    'units','data','HorizontalAlignment','left','VerticalAlignment','middle','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
            end
            for BFind=1:length(BFs_kHz{ROWind,ATTind})
                if ~isempty(PERhistXs_sec{ROWind,ATTind}{BFind})
                    if (FeatINDPlot<=FeatIND)
                        if (FeatINDPlot>1)
                            text(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot),YLIMITS(1),sprintf('1/%s',FeaturesText{FeatINDPlot}),'units','data', ...
                                'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6,'Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                        else
                            text(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot),YLIMITS(1),sprintf('1/%s',FeaturesText{FeatINDPlot}),'units','data', ...
                                'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6,'Color','k')
                        end
                    end
                end
            end
        end
        hold off
        
        
        %%%% Rate Plot
        PLOTnum=(ROWind-1)*NUMcols+2;
        eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
        semilogy(Rates{ROWind,ATTind},BFs_kHz{ROWind,ATTind},'*-','Color',ATTENcolors{ATTind})
        hold on
        %                semilogy(Nsps{ROWind,ATTind}/10,BFs_kHz{ROWind,ATTind},'m+','MarkerSize',4,'Color',ATTENcolors{ATTind})
        semilogy(ALSRs{ROWind,ATTind},unit.Info.BF_kHz,'go','MarkerSize',6,'Color',ATTENcolors{ATTind})
        for ATTind2=fliplr(find(Nattens_dB~=max(Nattens_dB)))
            semilogy(Rates{ROWind,ATTind2},BFs_kHz{ROWind,ATTind2},'*-','Color',ATTENcolors{ATTind2})
            %                   semilogy(Nsps{ROWind,ATTind2}/10,BFs_kHz{ROWind,ATTind2},'m+','MarkerSize',4,'Color',ATTENcolors{ATTind2})
            semilogy(ALSRs{ROWind,ATTind2},unit.Info.BF_kHz,'go','MarkerSize',6,'Color',ATTENcolors{ATTind2})
        end
        semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
        %                xlabel(sprintf('Rate (sp/sec)\n[+: # of spikes/10]\nO: ALSR'),'FontSize',6)
        xlabel(sprintf('Rate (sp/sec)\nO: ALSR'),'FontSize',8)
        PLOThand=eval(['h' num2str(PLOTnum)]);
        xlim(XLIMITS_rate)
        set(PLOThand,'XDir','reverse')
        set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
        ylim(YLIMITS)  % Same Ylimits for all plots
        %%%%%%%%%%%%%%%%%%%%%
        % Plot lines at all features
        for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
            if (yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000>=YLIMITS(1))&&(yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000<=YLIMITS(2))
                semilogy(XLIMITS_rate,yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
            end
        end
        hold off
        
        %%%% Synch Plot
        PLOTnum=(ROWind-1)*NUMcols+3;
        eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
        semilogy(Synchs{ROWind,ATTind},BFs_kHz{ROWind,ATTind},'*-','Color',ATTENcolors{ATTind})
        hold on
        for ATTind2=fliplr(find(Nattens_dB~=max(Nattens_dB)))
            semilogy(Synchs{ROWind,ATTind2},BFs_kHz{ROWind,ATTind2},'*-','Color',ATTENcolors{ATTind2})
        end
        semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
        xlabel(sprintf('Synch Coef (to %s)',FeaturesText{FeatIND}))
        PLOThand=eval(['h' num2str(PLOTnum)]);
        xlim(XLIMITS_synch)
        set(PLOThand,'XDir','reverse')
        set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
        set(gca,'XTick',[0 .25 .5 .75 1],'XTickLabel',{'0','.25','.5','.75','1'})
        ylim(YLIMITS)  % Same Ylimits for all plots
        %%%%%%%%%%%%%%%%%%%%%
        % Plot lines at all features
        for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
            if (yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000>=YLIMITS(1))&&(yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000<=YLIMITS(2))
                semilogy(XLIMITS_synch,yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
            end
        end
        hold off
        
        %%%% Rho Plot
        PLOTnum=(ROWind-1)*NUMcols+4;
        eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);']);
        PLOThand=eval(['h' num2str(PLOTnum)]);
        for ATTind=1:length(Nattens_dB)
            %                 [SCCpos,SCCs_belowBF,SCCs_aboveBF,BFind,BFs(COLind)]=getSCCindsWithBF(unit.Info.BF_kHz,NSAC_BFs_kHz{COLind,ATTind},NSCC_BFs_kHz{COLind,ATTind});
            semilogy(data2plot{ROWind,ATTind}(:,4),unit.Info.BF_kHz*2.^(data2plot{ROWind,ATTind}(:,1)),'*-','Color',ATTENcolors{ATTind});
            hold on;
        end
        xlabel(sprintf('Corr Coef (to %s)',FeaturesText{FeatIND}))
        xlim(XLIMITS_synch)
        set(PLOThand,'XDir','reverse')
        set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
        set(gca,'XTick',[0 .25 .5 .75 1],'XTickLabel',{'0','.25','.5','.75','1'})
        ylim(YLIMITS)  % Same Ylimits for all plots
        %%%%%%%%%%%%%%%%%%%%%
        % Plot lines at all features
        for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
            if (yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000>=YLIMITS(1))&&(yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000<=YLIMITS(2))
                semilogy(XLIMITS_synch,yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
            end
        end
        hold off
        
        %%%% CD Plot
        PLOTnum=(ROWind-1)*NUMcols+5;
        eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);']);
        PLOThand=eval(['h' num2str(PLOTnum)]);
        XLIMITS_cd=[0 0];
        for ATTind=1:length(Nattens_dB)
            %                 [SCCpos,SCCs_belowBF,SCCs_aboveBF,BFind,BFs(COLind)]=getSCCindsWithBF(unit.Info.BF_kHz,NSAC_BFs_kHz{COLind,ATTind},NSCC_BFs_kHz{COLind,ATTind});
            semilogy(data2plot{ROWind,ATTind}(:,2),unit.Info.BF_kHz*2.^(data2plot{ROWind,ATTind}(:,1)),'*-','Color',ATTENcolors{ATTind});
            XLIMITS_cd = [min(XLIMITS_cd(1),min(data2plot{ROWind,ATTind}(:,2))) max(XLIMITS_cd(2),max(data2plot{ROWind,ATTind}(:,2)))];
            hold on;
        end
        xlabel([sprintf('Delay, re %s',FeaturesText{FeatIND}) ' (ms)'])
        xlim(XLIMITS_cd)
        set(PLOThand,'XDir','reverse')
        set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
        ylim(YLIMITS)  % Same Ylimits for all plots
        %%%%%%%%%%%%%%%%%%%%%
        % Plot lines at all features
        for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
            if (yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000>=YLIMITS(1))&&(yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000<=YLIMITS(2))
                semilogy(XLIMITS_cd,yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
            end
        end
        hold off
        
    end %End if data for this condition, plot
    
end % End Feature


Xcorner=0.05;
Xwidth1=.4;
Xshift1=0.05;
Xwidth2=.1;
Xshift2=0.03;

Ycorner=0.05;
Yshift=0.07;
Ywidth=(1-NUMrows*(Yshift+.01))/NUMrows;   %.26 for 3; .42 for 2

TICKlength=0.02;

if NUMrows>2 % these positions will need to be edited
    set(h11,'Position',[Xcorner Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
    set(h12,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
    set(h13,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
    set(h14,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
    set(h15,'Position',[Xcorner+Xwidth1+Xshift1+3*(Xwidth2+Xshift2) Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
end
if NUMrows>1
    set(h6,'Position',[Xcorner Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
    set(h7,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
    set(h8,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
    set(h9,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
    set(h10,'Position',[Xcorner+Xwidth1+Xshift1+3*(Xwidth2+Xshift2) Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
end
set(h1,'Position',[Xcorner Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
set(h2,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
set(h3,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
set(h4,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
set(h5,'Position',[Xcorner+Xwidth1+Xshift1+3*(Xwidth2+Xshift2) Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])

orient landscape


% EOF
