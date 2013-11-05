%% plotRateALSR
F1_index = find(strncmpi('F1',FeaturesText,2));
F2_index = find(strncmpi('F2',FeaturesText,2));
T1_index = find(strncmpi('T1',FeaturesText,2));
T2_index = find(strncmpi('T2',FeaturesText,2));
 
figure, hold on;
% ['k','r','m','g'] for [normal,impaired,linear,nonlinear]   
for thisColor = 'r'%'krmg'
    index=0;
    SNR_alsr = NaN*ones(numel(SNR)*3);
    slope_rate = NaN*ones(numel(SNR),3,2);
    slope_alsr = NaN*ones(numel(SNR),3,2);
    for ii=1:numel(SNR)
        if strcmpi(colors{ii},thisColor)
            for jj=1:numel(SNR{ii})
                index=index+1;
                SNR_alsr(index)=SNR{ii}(jj);
                slope_rate(ii,jj,1)=SMP_rate{ii}{jj}(F1_index)-SMP_rate{ii}{jj}(T1_index);
                slope_rate(ii,jj,2)=SMP_rate{ii}{jj}(F2_index)-SMP_rate{ii}{jj}(T2_index);
                slope_alsr(ii,jj,1)=SMP_alsr{ii}{jj}(F1_index)-SMP_alsr{ii}{jj}(T1_index);
                slope_alsr(ii,jj,2)=SMP_alsr{ii}{jj}(F2_index)-SMP_alsr{ii}{jj}(T2_index);
            end
        end
    end
    
    for ii=1:3
        temp = slope_rate(:,ii,:);
        temp = temp(~isnan(temp));
        mean_slope_rate(ii) = mean(temp);
        std_slope_rate(ii) = std(temp)/sqrt(numel(temp));
    end
    for ii=1:3
        temp = slope_alsr(:,ii,:);
        temp = temp(~isnan(temp));
        mean_slope_alsr(ii) = mean(temp);
        std_slope_alsr(ii) = std(temp)/sqrt(numel(temp));
    end
    errorbar((1:3)+0.05,mean_slope_rate,std_slope_rate,[thisColor 'o-'],'LineWidth',3,'MarkerSize',10);
    errorbar((1:3)+0.15,mean_slope_alsr,std_slope_alsr,[thisColor 'sq--'],'LineWidth',3,'MarkerSize',10);
end
xlim([0.5 3.5]);
set(gca,'XDir','reverse');
set(gca,'FontSize',14,'XTick',[1 2 3]);
set(gca,'XTickLabel',{'Equal SL','Equal SPL','In Quiet'});
title('Strength of Vowel Coding','FontSize',14);
ylabel('spikes/sec (formant-vs-trough)','FontSize',14);
if color=='k', legend('Rate','ALSR'); end
hold off;
