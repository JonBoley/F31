% Compare multiple CFs (STMP vs. actual CF)

%%%
% add ability to calculate error_rho as a function of CF 
%   - (correlate to phase locking?)
% also, add ability to calculate error in CD as a function of CF
%%%

% filesSTMP = dir('STMPvsCF_STMP*.mat');
% filesCF = dir('STMPvsCF_CF*.mat');
% numFiles = length(-2:0.5:2);
% filenamesSTMP = {filesSTMP(end-numFiles+1:end).name};
% filenamesCF = {filesCF(end-numFiles+1:end).name};

filenamesCF = {...
%     'STMPvsCF_CF_2012-07-27_125539.mat',... %425Hz, 100CFs, +/-1oct
%     'STMPvsCF_CF_2012-07-27_152137.mat',... %601Hz, 100CFs, +/-1oct
%     'STMPvsCF_CF_2012-07-27_174637.mat',... %850Hz, 100CFs, +/-1oct
%     'STMPvsCF_CF_2012-07-27_205755.mat',... %1202Hz, 100CFs, +/-1oct
%     'STMPvsCF_CF_2012-07-31_101231.mat',... %1700Hz, 100CFs, +/-1oct
%     'STMPvsCF_CF_2012-07-31_134705.mat',... %2404Hz, 100CFs, +/-1oct
%     'STMPvsCF_CF_2012-07-31_165605.mat',... %3400Hz, 100CFs, +/-1oct
%     'STMPvsCF_CF_2012-07-31_202826.mat',... %4808Hz, 100CFs, +/-1oct
%     'STMPvsCF_CF_2012-08-01_081428.mat',... %6800Hz, 100CFs, +/-1oct

%     'STMPvsCF_CF_2012-08-03_165823.mat',... %6800Hz, 100CFs, +/-1oct
%     'STMPvsCF_CF_2012-08-04_134726.mat',... %6800Hz, 100CFs, +/-1oct

%     'STMPvsCF_CF_2012-08-03_122331.mat',... %3400Hz, 100CFs, +/-1oct
%     'STMPvsCF_CF_2012-08-03_144449.mat',... %4808Hz, 100CFs, +/-1oct
%     'STMPvsCF_CF_2012-08-04_161201.mat',... %2404Hz, 100CFs, +/-1oct
%     'STMPvsCF_CF_2012-08-04_180859.mat',... %1700Hz, 100CFs, +/-1oct
%     'STMPvsCF_CF_2012-08-04_203253.mat',... %1202Hz, 100CFs, +/-1oct
%     'STMPvsCF_CF_2012-08-05_093428.mat',... %850Hz, 100CFs, +/-1oct
%     'STMPvsCF_CF_2012-08-05_121329.mat',... %601Hz, 100CFs, +/-1oct
%     'STMPvsCF_CF_2012-08-05_140905.mat',... %425Hz, 100CFs, +/-1oct

%     'STMPvsCF_CF_2012-08-06_095202.mat'... %1.7kHz, 20CFs
    'STMPvsCF_CF_2012-08-06_103432.mat'... %1.7kHz, 20CFs
    };

filenamesSTMP = {...
%     'STMPvsCF_STMP_2012-07-27_140934.mat',... %425Hz, 100CFs, +/-1oct, 2.6ms delay
%     'STMPvsCF_STMP_2012-07-27_163552.mat',... %601Hz, 100CFs, +/-1oct, 2.4ms delay
%     'STMPvsCF_STMP_2012-07-27_185954.mat',... %850Hz, 100CFs, +/-1oct, 2.2ms delay
%     'STMPvsCF_STMP_2012-07-27_221116.mat',... %1202Hz, 100CFs, +/-1oct, 2.0ms delay
%     'STMPvsCF_STMP_2012-07-31_112836.mat',... %1700Hz, 100CFs, +/-1oct, 1.8ms delay
%     'STMPvsCF_STMP_2012-07-31_150208.mat',... %2404Hz, 100CFs, +/-1oct, 1.6ms delay
%     'STMPvsCF_STMP_2012-07-31_181137.mat',... %3400Hz, 100CFs, +/-1oct, 1.4ms delay
%     'STMPvsCF_STMP_2012-07-31_212716.mat',... %4808Hz, 100CFs, +/-1oct, 1.3ms delay
%     'STMPvsCF_STMP_2012-08-01_091331.mat',... %6800Hz, 100CFs, +/-1oct, 1.2ms delay

%     'STMPvsCF_STMP_2012-08-03_175728.mat',... %6800Hz, 100CFs, +/-1oct, 1.2ms delay
%     'STMPvsCF_STMP_2012-08-04_144613.mat',... %6800Hz, 100CFs, +/-1oct, 1.2ms delay

%     'STMPvsCF_STMP_2012-08-03_132244.mat',... %3400Hz, 100CFs, +/-1oct, 1.4ms delay
%     'STMPvsCF_STMP_2012-08-03_155145.mat',... %4808Hz, 100CFs, +/-1oct, 1.3ms delay
%     'STMPvsCF_STMP_2012-08-04_171046.mat',... %2404Hz, 100CFs, +/-1oct, 1.6ms delay
%     'STMPvsCF_STMP_2012-08-04_190729.mat',... %1700Hz, 100CFs, +/-1oct, 1.8ms delay
%     'STMPvsCF_STMP_2012-08-04_213107.mat',... %1202Hz, 100CFs, +/-1oct, 2.0ms delay
%     'STMPvsCF_STMP_2012-08-05_103302.mat',... %850Hz, 100CFs, +/-1oct, 2.2ms delay
%     'STMPvsCF_STMP_2012-08-05_131152.mat',... %601Hz, 100CFs, +/-1oct, 2.4ms delay
%     'STMPvsCF_STMP_2012-08-05_150746.mat'... %425Hz, 100CFs, +/-1oct, 2.6ms delay
    
%     'STMPvsCF_STMP_2012-08-06_104625.mat'... %1.7kHz, 20CFs, 1ms neural delay
    'STMPvsCF_STMP_2012-08-06_100358.mat'... %1.7kHz, 20CFs, CF-dependent delay
    };

numFiles = length(filenamesSTMP);

%%
SCCscaleFactor=0.05;
colorList='rgbkrgbkrgbkrbgkrgbk';
figure(1), hold on;
set(gca,'XScale','log','XTick',[0.125 0.25 0.5 1 2 4 8],'XTickLabel',{'125','250','500','1k','2k','4k','8k'});
xlim([0.125 16]);
figure(2), hold on;
set(gca,'XScale','log','XTick',[0.125 0.25 0.5 1 2 4 8],'XTickLabel',{'125','250','500','1k','2k','4k','8k'});
xlim([0.125 16]);
for i=1:numFiles
    load(filenamesCF{i},'-regexp','[^i^hRho]');
    figure(1),
    hRho(1)=plot(midCF_kHz*2.^deltaCF(2:end),squeeze(Rho_vDeltaCF(:,:,2:end)),...
        [colorList(i) '-'],'LineWidth',3);  
    figure(2), 
    hCD(1)=plot(midCF_kHz*2.^deltaCF(2:end),medfilt1(squeeze(CD_vDeltaCF(:,:,2:end)),15),...
        [colorList(i) '-']);%,'LineWidth',3);  
    if exist('SCC','var')
        figure(3), hold on;
        for j=1:length(deltaCF)
            hSCC(1)=plot(SCCdelays_usec{j,:,:},...
                midCF_kHz*2.^deltaCF(j)+SCC{j,:,:}*SCCscaleFactor,...
                [colorList(i) '-']);
        end
        hold off;
    end
    
    load(filenamesSTMP{i},'-regexp','[^i^hRho]');
    figure(1),
    hRho(2)=plot(midCF_kHz*2.^deltaCF(2:end),squeeze(Rho_vDeltaCF(:,:,2:end)),...
        [colorList(i) '-.'],'LineWidth',3);
    figure(2),
    hCD(2)=plot(midCF_kHz*2.^deltaCF(2:end),medfilt1(squeeze(CD_vDeltaCF(:,:,2:end)),15),...
        [colorList(i) ':']);%,'LineWidth',3);  
    if exist('SCC','var')
        figure(3), hold on;
        for j=1:length(deltaCF)
            hSCC(2)=plot(SCCdelays_usec{j,:,:},...
                midCF_kHz*2.^deltaCF(j)+SCC{j,:,:}*SCCscaleFactor,...
                [colorList(i) ':']);
        end
        hold off;
    end
    NeuralDelays(i,:) = [NeuralDelay_sec midCF_kHz];
    pause(0.1);
end

figure(1), legend(hRho,{'Actual CF','STMP'});
xlabel('CF (Hz)'); ylabel('\rho (re F_2)');

h1=gca; h2=axes('Position',get(h1,'Position'));
plot(h2,NeuralDelays(:,2),1e3*NeuralDelays(:,1),'o','Color',[0.5 0.5 0.5],'LineWidth',3);
set(h2,'YAxisLocation','right','Color','none');
set(h2,'XScale','log','XTick',[],'XTickLabel',[]);
set(h2,'XLim',get(h1,'XLim'),'Layer','top');
ylim([1 3]); ylabel('Neural Latency (ms)');
axes(h1);

figure(3); 
set(gca,'YScale','log','YTick',[0.125 0.25 0.5 1 2 4 8],'YTickLabel',{'125','250','500','1k','2k','4k','8k'});
ylim([min(midCF_kHz*2.^deltaCF) max(midCF_kHz*2.^deltaCF)]);
legend(hSCC,{'Actual CF','STMP'});

% figure, plot(NeuralDelays(:,2),1e3*NeuralDelays(:,1));
% set(gca,'XScale','log','XTick',[0.25 0.5 1 2 4 8],'XTickLabel',{'250','500','1k','2k','4k','8k'});
% xlim([0.25 8]); ylim([0 3]);
% xlabel('CF (Hz)'); 


