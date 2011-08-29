function convertabr (prefix,first,last,chamber)

%call as convertabr(prefix,first,last,chamber)
%
%   prefix:  letter at beginning of old file
%   first:   first file number to be converted
%   last:    last file number to be converted
%   chamber: system where ABR was collected
%
%   example: convertabr('a',1,12,'north')
%
%   files must be in the present working directory (pwd)


if nargin < 2,
    warndlg('Must supply input filenumber.','Input Error');
    return
end
if nargin < 3
    last = first;
end
if nargin < 4
    chamber = 'ross';
end

for picnum = first:last
    
    if strcmp('north',chamber)
        ReadABRNorth;
    else
        command = sprintf('%s%s%d%s','eval(''',prefix,picnum,''');');
        eval(command);
        file_samp_rate = 12.2;
    end
    
    buf = zeros(length(file_data),4);
    buf(:,1) = file_data(:,1);
    buf(:,4) = file_data(:,2);
    x = struct('General', {struct('program_name', {file_prog } ...
            ,'picture_number', {picnum } ...
            ,'date', {file_date } ...
            ,'time', {'xx:xx:xx'} ...
            ,'comment', {'Converted from old format.'} ...
            )} ...
        ,'Stimuli', {struct('freq_hz', {file_stim_freq } ...
            ,'play_duration', {file_stim_dur } ...
            ,'record_duration', {file_rec_wind } ...
            ,'SampRate', {file_samp_rate } ...
            ,'pulses_per_sec', {file_stim_bursts } ...
            ,'rise_fall', {file_rise_fall } ...
            ,'naves', {file_num_samps } ...
            ,'db_atten', {file_stim_atten } ...
            ,'amp_gain', {file_samp_gain } ...
            )} ...
        ,'AverageData', {[buf]} ...
        ,'User', {'hwf' } ...
        ,'Hardware', {struct('amp_vlt', {file_amp_vlt } ...
            )} ...
        );
    
    if picnum < 10
        pic = ['00' num2str(picnum)];
    elseif picnum < 100
        pic = ['0' num2str(picnum)];
    else
        pic = [num2str(picnum)];
    end
    
    filename = ['EPavg_p' pic];
    fid = fopen([filename '.m'],'wt+');
    fprintf(fid,'function x = %s\n',filename);
    mat2text(x,fid);
    fclose(fid);
    
end