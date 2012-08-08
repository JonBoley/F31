%% batch_STMPplot_SCCfunctions
if ~exist('ExpList','var')
    setup_Vowel_STMP;
end

% dates = {'062311','062311','062311','062311',...
%     '072111','072111','072111','072111','072111','072111',...
%     '080111','080111','080111','080111','080111',...
%     '080911','080911','080911','080911','080911','080911','080911','080911'...
%     };
% unitNums = {1.01,1.09,1.11,1.14,...
%     1.01, 1.06, 1.11, 1.12, 1.15, 1.16,...
%     2.02, 4.01, 4.02, 4.03, 4.05,...
%     1.04, 1.05, 1.06, 1.07, 1.09, 1.11, 1.15, 1.16...
%     };
date = '050412';
unitNums = {1.01, 1.02, 1.03, 1.05, 2.01, 2.03};
dates=repmat({date},size(unitNums));


for i=1:length(dates)
    count=fprintf('%d of %d',i,length(dates));
    for FeatureIndices=1:2%4
        for AttenIndices=1:3
            % just plot
%             STMPplot_SCCfunctions(dates{i},num2str(unitNums{i}),FeatureIndices,AttenIndices);

            % interactive plot (lets you manually fix peaks)
            calcCDreBF_manual({dates{i}},num2str(unitNums{i}),FeatureIndices,AttenIndices);
        end
    end

%     pause;
    fprintf(repmat('\b',1,count));
end
disp('Done!');
