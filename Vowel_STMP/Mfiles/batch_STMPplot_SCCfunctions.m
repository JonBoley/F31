%% batch_STMPplot_SCCfunctions
if ~exist('ExpList','var')
    setup_Vowel_STMP;
end

dates = {'062311','062311','062311','062311',...
    '072111','072111','072111','072111','072111','072111',...
    '080111','080111','080111','080111','080111',...
    '080911','080911','080911','080911','080911','080911','080911','080911'...
    };
unitNums = {1.01,1.09,1.11,1.14,...
    1.01, 1.06, 1.11, 1.12, 1.15, 1.16,...
    2.02, 4.01, 4.02, 4.03, 4.05,...
    1.04, 1.05, 1.06, 1.07, 1.09, 1.11, 1.15, 1.16...
    };

for i=1%:length(dates)
    FeatureIndices=3;
    AttenIndices=1;
    
    % just plot
%     STMPplot_SCCfunctions(dates{i},num2str(unitNums{i}),FeatureIndices,AttenIndices);
    
    % interactive plot (lets you manually fix peaks)
    if exist('h_calcCDreBF_manual','var'), uiwait(h_calcCDreBF_manual); end; 
    h_calcCDreBF_manual = calcCDreBF_manual(dates{i},num2str(unitNums{i}),FeatureIndices,AttenIndices);
    
%     count=fprintf('%d of %d',i,length(dates));
%     pause; fprintf(repmat('\b',1,count));
end
clear h_calcCDreBF_manual;
