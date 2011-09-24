% plot_abr_shifts

abr_shift_files = dir('C:\Research\MATLAB\Vowel_STMP\ExpData\*_abr_shift.mat');
exclude = {'chin1103','chin1100'};

figure(123); hold on;
colorr=['k';'r';'b';'m';'g';'k';'r';'b';'m';'g';'k';'r';'b';'m';'g'];

freqs_all=[];
shift_db_all=[];

for i=1:length(abr_shift_files)
    load(abr_shift_files(i).name);
    plot(freqs,shift_db,['*-' colorr(i)]);
    
    ignore=0;
    for j=1:length(exclude)
        if strfind(abr_shift_files(i).name,exclude{j})
            ignore=1;
        end
    end

    if ~ignore
        if any(shift_db>0)
            fprintf('%s contains positive shifts!\n',abr_shift_files(i).name);
            shift_db(shift_db>0)=NaN;
        end
        
        if isempty(freqs_all)
            freqs_all=freqs;
            shift_db_all=shift_db;
        else
            for j=1:length(freqs_all)
                index = find(freqs==freqs_all);
                shift_db_all(index,i) = shift_db;
            end
        end
    end
end

% calculate average shift
for i=1:length(freqs_all)
    shift_db_avg(i) = mean(shift_db_all(i,~isnan(shift_db_all(i,:))));
end
plot(freqs_all,shift_db_avg,['o-' colorr(i)],'LineWidth',3);

figure(123); hold off;
xlabel('frequency (Hz)');
ylabel('threshold shift (-dB HL)');
set(gca,'XScale','log'); xlim([min(freqs_all) max(freqs_all)]);

