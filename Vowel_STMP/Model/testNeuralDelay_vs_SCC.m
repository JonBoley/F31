%% evaluate the effects of neural delay
% CF, 1ms, 2ms, 3ms, 4ms
filename = {...
    'STMPvsCF_CF_2012-10-06_162912',...
    'STMPvsCF_STMP_2012-10-06_163854'};
%     'STMPvsCF_CF_2012-09-22_141530',...
%     'STMPvsCF_STMP_2012-09-22_141713',...
%     'STMPvsCF_STMP_2012-09-22_145132',...
%     'STMPvsCF_STMP_2012-09-22_150857',...
%     'STMPvsCF_STMP_2012-09-22_142549'};
 
figure(99); colors='rbgkm'; plotNum=0;
for fname=filename
    plotNum=plotNum+1;
    load(fname{1});
    disp(fname{1});

    scaleSize = 0.2;
    for i=1:length(SCC)
        figure(99);
        hSCC(plotNum)=plot(SCCdelays_usec{i}/1e6,CF_kHz(i)+scaleSize*SCC{i},colors(plotNum)); hold on;
    end
    if exist('NeuralDelay_sec'), hSCC_NeuralDelay(plotNum) = NeuralDelay_sec; end
end
% xlim([-2 2]/1e3);
legend(hSCC,num2str(hSCC_NeuralDelay'));
