%%
CFrange = [0.5 2];%4*[2^(-0.5) 2^(0.5)]; %kHz
FeatureNum = [1 2];
SNRindex = [3]; % quiet, spl, sl

NHindex = find(nhCFs>=CFrange(1) & nhCFs<=CFrange(end));
SNHLindex = find(snhlCFs>=CFrange(1) & snhlCFs<=CFrange(end));
NALindex = find(linCFs>=CFrange(1) & linCFs<=CFrange(end));
DSLindex = find(nonlinCFs>=CFrange(1) & nonlinCFs<=CFrange(end));

matNames = {'NH','SNHL','Linear','NonLinear'};
matrix{1}=[]; matrix{2}=[]; matrix{3}=[]; matrix{4}=[];
for m=FeatureNum
    for n=SNRindex
        matrix{1} = [matrix{1} array_cd.nh{m,n}(NHindex)];
        matrix{2} = [matrix{2} array_cd.snhl{m,n}(SNHLindex)];
        matrix{3} = [matrix{3} array_cd.lin{m,n}(NALindex)];
        matrix{4} = [matrix{4} array_cd.nonlin{m,n}(DSLindex)];
    end
end

% matrix{1}=[]; matrix{2}=[]; matrix{3}=[];
% for m=FeatureNum
%     for n=SNRindex
%         matrix{n} = [matrix{n} array_cd.snhl{m,n}(SNHLindex)];
%     end
% end

for refIndex = 1:4
    for testIndex = (refIndex+1):4
        [p_rs,h_rs,stats_rs] = ranksum(matrix{refIndex}(~isnan(matrix{refIndex})),...
            matrix{testIndex}(~isnan(matrix{testIndex})));
        
        fprintf('%s/%s...   ',matNames{refIndex},matNames{testIndex});
        if h_rs
            fprintf('Significant difference found! (p=%1.5f)\n',p_rs);
        else
            fprintf('No significant difference found (p=%1.5f)\n',p_rs);
        end
    end
end


figure, hold on;
% CD for {NH,SNHL,lin,nonlin}, in quiet then in noise
% errorbar([1:4]+0.05,[median(matrix{1}(~isnan(matrix{1}))),...
%     median(matrix{2}(~isnan(matrix{2}))),...
%     median(matrix{3}(~isnan(matrix{3}))),...
%     median(matrix{4}(~isnan(matrix{4})))],...
%     [prctile(matrix{1}(~isnan(matrix{1})),0.25),...
%     prctile(matrix{2}(~isnan(matrix{2})),0.25),...
%     prctile(matrix{3}(~isnan(matrix{3})),0.25),...
%     prctile(matrix{4}(~isnan(matrix{4})),0.25)],...
%     [prctile(matrix{1}(~isnan(matrix{1})),0.75),...
%     prctile(matrix{2}(~isnan(matrix{2})),0.75),...
%     prctile(matrix{3}(~isnan(matrix{3})),0.75),...
%     prctile(matrix{4}(~isnan(matrix{4})),0.75)],'bsq:');
errorbar([1:4]-0.05,[mean(matrix{1}(~isnan(matrix{1}))),...
    mean(matrix{2}(~isnan(matrix{2}))),...
    mean(matrix{3}(~isnan(matrix{3}))),...
    mean(matrix{4}(~isnan(matrix{4})))],...
    [stdev(matrix{1}(~isnan(matrix{1})))/sqrt(numel(matrix{1}(~isnan(matrix{1})))),...
    stdev(matrix{2}(~isnan(matrix{2})))/sqrt(numel(matrix{2}(~isnan(matrix{2})))),...
    stdev(matrix{3}(~isnan(matrix{3})))/sqrt(numel(matrix{3}(~isnan(matrix{3})))),...
    stdev(matrix{4}(~isnan(matrix{4})))/sqrt(numel(matrix{4}(~isnan(matrix{4}))))],'ko-');
set(gca,'XTick',[1 2 3 4]);
set(gca,'XTickLabel',{'Normal','Impaired','Linear Gain','Nonlinear Gain'});
ylabel('CD (CF cycles) [0.5 oct re F_{X}]');

%% rho for normal vs. impaired
figure, hold on;
errorbar(1,mean(mean(cell2mat(rhoCFgroups.nh))),(stdev(cell2mat(rhoCFgroups.nh))))
errorbar(2,mean(mean(cell2mat(rhoCFgroups.snhl))),(stdev(cell2mat(rhoCFgroups.snhl))))
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'Normal','Impaired'});
ylabel('\rho width (octaves)');
[p_rs,h_rs,stats_rs] = ranksum(reshape(cell2mat(rhoCFgroups.nh),numel(cell2mat(rhoCFgroups.nh)),1),...
    reshape(cell2mat(rhoCFgroups.snhl),numel(cell2mat(rhoCFgroups.snhl)),1));

%% rho for normal, quiet vs. noise
figure, hold on;
errorbar(0.9,mean(mean(cell2mat({rhoCFgroups.nh{:,1}}))),(stdev(cell2mat({rhoCFgroups.nh{:,1}}))),'k+')
errorbar(1.9,mean(mean(cell2mat({rhoCFgroups.nh{:,2:3}}))),(stdev(cell2mat({rhoCFgroups.nh{:,2:3}}))),'k+')
errorbar(1.1,mean(mean(cell2mat({rhoCFgroups.snhl{:,1}}))),(stdev(cell2mat({rhoCFgroups.snhl{:,1}}))),'ro')
errorbar(2.1,mean(mean(cell2mat({rhoCFgroups.snhl{:,2:3}}))),(stdev(cell2mat({rhoCFgroups.snhl{:,2:3}}))),'ro')
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'Quiet','In Noise'});
ylabel('\rho width (octaves)');
title('Normal Hearing');
[p_rs,h_rs,stats_rs] = ranksum(reshape(cell2mat({rhoCFgroups.nh{:,1}}),numel(cell2mat({rhoCFgroups.nh{:,1}})),1),...
    reshape(cell2mat({rhoCFgroups.nh{:,2:3}}),numel(cell2mat({rhoCFgroups.nh{:,2:3}})),1));

