% regressionTest

load('C:\Research\MATLAB\Vowel_STMP\ExpData\regressionTest.mat');

hlConfigs = {'nh','snhl'};
colors = ['k','r'];
FeatureNums = [1 3]; %F1,F2

figure(1); clf; hold on;

for i=1:length(hlConfigs)
    % concatenate all data
    eval(sprintf('%s.CF=[]; %s.CD=[]; %s.FeatureNum=[]; %s.SNRindex=[];',...
        hlConfigs{i},hlConfigs{i},hlConfigs{i},hlConfigs{i}));
    eval(sprintf('NumSNRs=size(arrayCDatHalfOct.%s,2);',hlConfigs{i}));
    for FeatureNum=FeatureNums
        for indx_snr=1:NumSNRs
            eval(sprintf('%s.CF = [%s.CF; arrayCDatHalfOct.%s{find(FeatureNum==FeatureNums),indx_snr}(:,1)];',...
                hlConfigs{i},hlConfigs{i},hlConfigs{i}));
            eval(sprintf('%s.CD = [%s.CD; arrayCDatHalfOct.%s{find(FeatureNum==FeatureNums),indx_snr}(:,2)];',...
                hlConfigs{i},hlConfigs{i},hlConfigs{i}));
            
            eval(sprintf('%s.FeatureNum = [%s.FeatureNum; FeatureNum*ones(size(arrayCDatHalfOct.%s{find(FeatureNum==FeatureNums),indx_snr}(:,2)))];',...
                hlConfigs{i},hlConfigs{i},hlConfigs{i}));
            eval(sprintf('%s.SNRindex = [%s.SNRindex; indx_snr*ones(size(arrayCDatHalfOct.%s{find(FeatureNum==FeatureNums),indx_snr}(:,2)))];',...
                hlConfigs{i},hlConfigs{i},hlConfigs{i}));
        end
    end
    
    
    eval(sprintf('temp=%s;',hlConfigs{i}));
    
    % sort the data
    [temp.CF,IX] = sort(temp.CF);
    temp.CD = temp.CD(IX);
    temp.FeatureNum = temp.FeatureNum(IX);
    temp.SNRindex = temp.SNRindex(IX);
    
    % throw out NaN's
    maxIX = min(find(isnan(temp.CF)))-1;
    temp.CF = temp.CF(1:maxIX);
    temp.CD = temp.CD(1:maxIX); 
    temp.FeatureNum = temp.FeatureNum(1:maxIX);
    temp.SNRindex = temp.SNRindex(1:maxIX); 
    
    % basic linear regression
    distr = 'normal'; % 'normal','poisson','gamma','binomial','inverse gaussian'
    link = 'log'; % 'identity','log','loglog','logit','probit','comploglog'
    [temp.b,temp.dev,temp.stats] = glmfit(temp.CF,temp.CD,distr,'link',link);
    [temp.yhat,temp.dylo,temp.dyhi] = glmval(temp.b,temp.CF,link,temp.stats);
    
    % multiple linear regression
    X = [ones(size(temp.CF)), temp.CF, temp.FeatureNum, temp.SNRindex, ...
        temp.CF.*temp.FeatureNum, temp.FeatureNum.*temp.SNRindex, ...
        temp.CF.*temp.SNRindex, temp.CF.*temp.FeatureNum.*temp.SNRindex];
    [temp.b_mlr,temp.bint,temp.r,temp.rint,temp.stats_mlr] = regress(temp.CD,X);
    % figure(2), rcoplot(temp.r,temp.rint); % plot residual (potential outliers in red)
    xfit2 = 0.1*2.^(0:0.1:6); % 100Hz + 6 octaves
    xfit3 = unique(X(:,3));
    xfit4 = unique(X(:,4));
    [XFIT2,XFIT3,XFIT4] = ndgrid(xfit2,xfit3,xfit4);
    YFIT = temp.b_mlr(1) ...
        + temp.b_mlr(2)*XFIT2 ...
        + temp.b_mlr(3)*XFIT3 ...
        + temp.b_mlr(4)*XFIT4 ...
        + temp.b_mlr(5)*XFIT2.*XFIT3 ...
        + temp.b_mlr(6)*XFIT3.*XFIT4 ...
        + temp.b_mlr(7)*XFIT2.*XFIT4 ...
        + temp.b_mlr(8)*XFIT2.*XFIT3.*XFIT4;
        
    
    figure(1), scatter(gca,temp.CF,temp.CD,[colors(i) '.']); 
%     set(gca,'XScale','log'); % log plots do not support OpenGL for alpha
    set(gca,'XTick',[1 2 4]);
    xlim([0.1 10]); ylim([0 10]);
                
    hold on;
    ha = area(temp.CF, [temp.yhat, (temp.dyhi+temp.dylo)]);
    set(ha(1), 'FaceColor', 'none'); % this makes the bottom area invisible
    set(ha, 'LineStyle', 'none');
    set(ha(2), 'FaceColor', colors(i));
    alpha(0.5);
    
%     for j=1:size(YFIT,2)
%         for k=1:size(YFIT,3)
%             plot(xfit2,YFIT(:,j,k),colors(i),'LineWidth',3); hold on;
%         end
%     end
    
    temp.figHandle = ha(2);
    eval(sprintf('%s=temp;',hlConfigs{i}));
end
xlabel('CF (kHz)'); ylabel('CD (CF cycles)');
legend([nh.figHandle snhl.figHandle],{'NH (95% CI)','SNHL (95% CI)'});
  
%% remove normal CF trend
nh.CD = nh.CD - nh.yhat;
snhl.CD = snhl.CD - interp1(nh.CF,nh.yhat,snhl.CF);
%% Non-parametric statistics on CD
for testNum=[5]
    KWmatrix{testNum}=[];
    RSmatrix{testNum,1}=[]; RSmatrix{testNum,2}=[];
    switch testNum
        case 1 % Normal vs. Impaired
            temp_nh=nh.CD(~isnan(nh.CD));
            temp_snhl=snhl.CD(~isnan(snhl.CD)); % get rid of extraneous NaNs
            minUnits = min(length(temp_nh),length(temp_nh));
            KWmatrix{testNum} = [[temp_nh;NaN*ones(max(0,length(temp_snhl)-length(temp_nh)),1)],...
                [temp_snhl;NaN*ones(max(0,length(temp_nh)-length(temp_snhl)),1)]];
            RSmatrix{testNum,1} = temp_nh; RSmatrix{testNum,2} = temp_snhl;
            groupNames = {'Normal','Impaired'};
            
        case 2 % Different CFs (normal)
            CFlist = unique(nh.CF); 
            numDataPoints=zeros(1,length(CFlist));
            KWmatrix{testNum}=NaN*ones(100,length(CFlist)); % assume <100 per CF
            for CFnum=1:length(CFlist)
                index = find(nh.CF==CFlist(CFnum));
                if ~isempty(index)
                    startIndex = numDataPoints(CFnum)+1;
                    endIndex = startIndex + length(index) - 1;
                    KWmatrix{testNum}(startIndex:endIndex,CFnum) = nh.CD(index);
                end
                numDataPoints(CFnum) = numDataPoints(CFnum) + length(index);
            end
            groupNames = CFlist;
            
        case 3 % Different CFs (impaired)
            CFlist = unique(snhl.CF);
            numDataPoints=zeros(1,length(CFlist));
            KWmatrix{testNum}=NaN*ones(100,length(CFlist));
            for CFnum=1:length(CFlist)
                index = find(snhl.CF==CFlist(CFnum));
                if ~isempty(index)
                    startIndex = numDataPoints(CFnum)+1;
                    endIndex = startIndex + length(index) - 1;
                    KWmatrix{testNum}(startIndex:endIndex,CFnum) = snhl.CD(index);
                end
                numDataPoints(CFnum) = numDataPoints(CFnum) + length(index);
            end
            groupNames = CFlist;
            
        case 4 % F1 vs F2
            FeatureNums = unique([nh.FeatureNum;snhl.FeatureNum]);
            numDataPoints=zeros(1,length(FeatureNums));
            KWmatrix{testNum}=NaN*ones(...
                max(length(nh.FeatureNum),length(snhl.FeatureNum)),length(FeatureNums));
            for i=1:length(FeatureNums)
                index = find(nh.FeatureNum==FeatureNums(i));
                if ~isempty(index)
                    startIndex = numDataPoints(i)+1;
                    endIndex = startIndex + length(index) - 1;
                    KWmatrix{testNum}(startIndex:endIndex,i) = nh.CD(index);
                end
                numDataPoints(i) = numDataPoints(i) + length(index);
                
                index = find(snhl.FeatureNum==FeatureNums(i));
                if ~isempty(index)
                    startIndex = numDataPoints(i)+1;
                    endIndex = startIndex + length(index) - 1;
                    KWmatrix{testNum}(startIndex:endIndex,i) = snhl.CD(index);
                end
                numDataPoints(i) = numDataPoints(i) + length(index);
                
                RSmatrix{testNum,1} = KWmatrix{testNum}(~isnan(KWmatrix{testNum}(:,1)),1); 
                RSmatrix{testNum,2} = KWmatrix{testNum}(~isnan(KWmatrix{testNum}(:,2)),2);
            end
            groupNames = {'F1','F2'};
            
        case 5 % Quiet / Equal SPL / Equal SL
            SNRnums = unique([nh.SNRindex;snhl.SNRindex]);
            numDataPoints=zeros(1,length(SNRnums));
            KWmatrix{testNum}=NaN*ones(...
                max(length(nh.SNRindex),length(snhl.SNRindex)),length(SNRnums));
            for i=1:length(SNRnums)
                index = find(nh.SNRindex==SNRnums(i));
                if ~isempty(index)
                    startIndex = numDataPoints(i)+1;
                    endIndex = startIndex + length(index) - 1;
                    KWmatrix{testNum}(startIndex:endIndex,i) = nh.CD(index);
                end
                numDataPoints(i) = numDataPoints(i) + length(index);
                
                index = find(snhl.SNRindex==SNRnums(i));
                if ~isempty(index)
                    startIndex = numDataPoints(i)+1;
                    endIndex = startIndex + length(index) - 1;
                    KWmatrix{testNum}(startIndex:endIndex,i) = snhl.CD(index);
                end
                numDataPoints(i) = numDataPoints(i) + length(index);
            end
            groupNames = {'Quiet','Equal SPL','Equal SL'};
            
            % Directly compare Quiet & Equal SPL
            RSmatrix{testNum,1} = KWmatrix{testNum}(~isnan(KWmatrix{testNum}(:,1)),1);
            RSmatrix{testNum,2} = KWmatrix{testNum}(~isnan(KWmatrix{testNum}(:,3)),3);
            
        case 6 % In CF Range: [Normal/Impaired] vs [Quiet/Equal SPL/Equal SL]
            
            CFrange = [0.1 2]; %kHz
            
            SNRnums = unique([nh.SNRindex;snhl.SNRindex]);
            numDataPoints=zeros(1,2*length(SNRnums));
            KWmatrix{testNum}=NaN*ones(...
                max(length(nh.SNRindex),length(snhl.SNRindex)),2*length(SNRnums));
            KWindex=0;
            for i=1:length(SNRnums)
                index = find(nh.SNRindex==SNRnums(i) ...
                    & nh.CF>CFrange(1) & nh.CF<CFrange(2) ...
                    );%& nh.FeatureNum==3);
                KWindex = KWindex + 1;
                if ~isempty(index)
                    startIndex = numDataPoints(KWindex)+1;
                    endIndex = startIndex + length(index) - 1;
                    KWmatrix{testNum}(startIndex:endIndex,KWindex) = nh.CD(index);
                end
                numDataPoints(KWindex) = numDataPoints(KWindex) + length(index);

                index = find(snhl.SNRindex==SNRnums(i) ...
                    & snhl.CF>CFrange(1) & snhl.CF<CFrange(2) ...
                    );%& snhl.FeatureNum==3);
                KWindex = KWindex + 1;
                if ~isempty(index)
                    startIndex = numDataPoints(KWindex)+1;
                    endIndex = startIndex + length(index) - 1;
                    KWmatrix{testNum}(startIndex:endIndex,KWindex) = snhl.CD(index);
                end
                numDataPoints(KWindex) = numDataPoints(KWindex) + length(index);
            end
            groupNames = {'Quiet (NH)','Quiet (SNHL)',...
                'Equal SPL (NH)','Equal SPL (SNHL)',...
                'Equal SL (NH)','Equal SL (SNHL)'};
            
            % Directly compare NH/SNHL @ SNRnum
            SNRnum = 3;
            RSmatrix{testNum,1} = KWmatrix{testNum}(~isnan(KWmatrix{testNum}(:,2*SNRnum-1)),2*SNRnum-1);
            RSmatrix{testNum,2} = KWmatrix{testNum}(~isnan(KWmatrix{testNum}(:,2*SNRnum)),2*SNRnum);

    end

    % Kruskal-Wallis Test (compares the median of each sample)
    % ... is each column from a different population?
    [p_kw{testNum},table_kw{testNum},stats_kw{testNum}] = kruskalwallis(KWmatrix{testNum},groupNames,'off');
%     figure;
%     [c_kw{testNum},m_kw{testNum},h_kw{testNum},gnames_kw{testNum}] = multcompare(stats_kw{testNum});
    
    % if only two groups, do Wilcoxon rank sum test (aka, Mann-Whitney U-test)
    if ~isempty(RSmatrix{testNum,1})
        [p_rs{testNum},h_rs{testNum},stats_rs{testNum}] = ranksum(RSmatrix{testNum,1},RSmatrix{testNum,2});
    end
end
p_rs

