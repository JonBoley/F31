% calc_abr_shift.m
% (to be run from within comp_plot.m)

k=1; m=1;
switch numfiles
    case 2
        k = menu('Select the normal data',dd(1).name,dd(2).name);
        m = menu('Select the most recent impaired data',dd(1).name,dd(2).name);
    case 3
        k = menu('Select the normal data',dd(1).name,dd(2).name,dd(3).name);
        m = menu('Select the most recent impaired data',dd(1).name,dd(2).name,dd(3).name);
    case 4
        k = menu('Select the normal data',dd(1).name,dd(2).name,dd(3).name,dd(4).name);
        m = menu('Select the most recent impaired data',dd(1).name,dd(2).name,dd(3).name,dd(4).name);
    otherwise
        disp('Code is not written to handle more than 4 ABRs (calc_abr_shift.m)');
end

freqs = dd(k).data.abrs.thresholds(:,1);
baseline = dd(k).data.abrs.thresholds(:,2);
thresh_spl = 150*ones(length(freqs),1);
for i=m%1:numfiles
    if i~=k
        for j=1:length(dd(i).data.abrs.thresholds(:,1))
            index = find(dd(i).data.abrs.thresholds(j,1)==dd(k).data.abrs.thresholds(:,1));
            if ~isempty(index)
                thresh_spl(j) = min(dd(i).data.abrs.thresholds(j,2),thresh_spl(j));
            end
        end
    end
end
thresh_spl(thresh_spl>149)=NaN;
shift_db = baseline-thresh_spl;
filename = inputdlg('Please supply a *.mat file name','',1,{[dd(1).name(1:8) '_abr_shift']});
save(['C:\Research\MATLAB\Vowel_STMP\ExpData\' filename{1}], 'freqs', 'shift_db');
