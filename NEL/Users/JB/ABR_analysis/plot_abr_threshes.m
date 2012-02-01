% plot_abr_threshes
switch lower(getenv('computername'))
    case {'m4500'}
        abr_data_dir = 'C:\Research\MATLAB\Vowel_STMP\ExpData\Summary\';
    otherwise
        abr_data_dir = 'C:\NEL\ExpData\Summary\';
end
abr_files = dir([abr_data_dir 'chin*.mat']);

% pick earliest as normal (& latest as impaired)
normal_files={}; impaired_files={};
for i=1:length(abr_files)
    Chin1 = str2num(abr_files(i).name(5:8));
    Date1 = datenum(abr_files(i).name(end-12:end-4), 'ddmmmyyyy');
    for j=i+1:length(abr_files)
        Chin2 = str2num(abr_files(j).name(5:8));
        Date2 = datenum(abr_files(j).name(end-12:end-4), 'ddmmmyyyy');
        if Chin1==Chin2
            if length(normal_files) && Chin1==str2num(normal_files{end}(5:8)) % already in normal_files
                if Date2<datenum(normal_files{end}(end-12:end-4), 'ddmmmyyyy')
                    % if earlier than normal, replace normal
                    normal_files{end} = abr_files(j).name;
                elseif Date2>datenum(impaired_files{end}(end-12:end-4), 'ddmmmyyyy')
                    % if later than impaired, replace impaired
                    impaired_files{end} = abr_files(j).name;
                end
            else % we have more than 2 files for this chin
                if Date1<Date2
                    normal_files{end+1} = abr_files(i).name;
                    impaired_files{end+1} = abr_files(j).name;
                else %Date2<=Date1
                    normal_files{end+1} = abr_files(j).name;
                    impaired_files{end+1} = abr_files(i).name;
                end
            end
        end
    end
end

figure(124); hold on;
colors=['k';'r']; %normal, impaired
cd(abr_data_dir);

for impaired=0:1 %normal, impaired
    if impaired, files = impaired_files;
    else files = normal_files;
    end

    freqs_all=[];
    thresh_spl_all=[];
    for i=1:length(files)
        load(files{i});

        %%%
        freqs = abrs.thresholds(:,1);
        baseline = abrs.thresholds(:,2);
        thresh_spl = 150*ones(length(freqs),1);
        for j=1:length(abrs.thresholds(:,1))
            index = find(abrs.thresholds(j,1)==abrs.thresholds(:,1));
            if ~isempty(index)
                thresh_spl(j) = min(abrs.thresholds(j,2),thresh_spl(j));
            end
        end
        thresh_spl(thresh_spl>149)=NaN;
        %%%
        
        [freqs,indx]=sort(freqs);
        thresh_spl=thresh_spl(indx);
        if impaired
            plot(freqs,thresh_spl,['*:' colors(impaired+1)]);
        else
            plot(freqs,thresh_spl,['*-' colors(impaired+1)]);
        end

        freqs_all = union(freqs_all,freqs);
        for j=1:length(freqs_all)
            for k=1:length(freqs)
                index = find(freqs(k)==freqs_all);
                thresh_spl_all(index,i) = thresh_spl(k);
            end
        end
    end % for all files
    thresh_spl_all(~thresh_spl_all)=NaN; % replace zeros

    % calculate average shift
    for i=1:length(freqs_all)
        thresh_spl_avg(i) = mean(thresh_spl_all(i,~isnan(thresh_spl_all(i,:))));
    end
    h(impaired+1)=plot(freqs_all,thresh_spl_avg,['o-' colors(impaired+1)],'LineWidth',3);
    
    if ~impaired
        thresh_shift_all = thresh_spl_all;
    else
        thresh_shift_all = thresh_spl_all - thresh_shift_all;
    end

end %impaired?

figure(124); hold off;
title('ABR Thresholds');
xlabel('frequency (Hz)');
ylabel('threshold (dB SPL)');
set(gca,'XScale','log'); xlim([min(freqs_all) max(freqs_all)]);
legend(h,{'normal','impaired'});

figure(125), hold on;
plot(freqs_all,thresh_shift_all,'k*-'); 
plot(freqs_all,nanmean(thresh_shift_all,2),'k-','LineWidth',3); hold off;
title('ABR Threshold Shifts');
xlabel('frequency (Hz)');
ylabel('threshold shift (dB)');
set(gca,'XScale','log'); xlim([min(freqs_all) max(freqs_all)]);


