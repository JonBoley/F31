function convertcal (prefix,first,last,chamber)

%call as convertcal(prefix,first,last,chamber)
%
%   prefix:  letter at beginning of old file
%   first:   first file number to be converted
%   last:    last file number to be converted
%   chamber: system where cal was performed
%
%   example: convertcal('a',1,12,'north')
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
        ReadCalNorth;
    else
        command = sprintf('%s%s%d%s','eval(''',prefix,picnum,''');');
        eval(command);
    end
    x = struct('General', {struct('program_name', {file_prog} ...
            ,'picture_number', {picnum } ...
            ,'date', {file_date} ...
            ,'time', {'xx:xx:xx'} ...
            ,'comment', {'Converted from old format.'} ...
            )} ...
        ,'Stimuli', {struct('frqlo', {file_frqlo } ...
            ,'frqhi', {file_frqhi } ...
            ,'frqcal', {nan } ...
            ,'syslv', {file_isyslv } ...
            )} ...
        ,'CalibData', {[file_data(:,1:2)]} ...
        ,'User', {'hwf'} ...
        ,'Hardware', {struct('mic', {file_mic} ...
            ,'amp_vlt', {file_amp_vlt } ...
            )} ...
        );
    
    if picnum < 10
        pic = ['00' num2str(picnum)];
    elseif picnum < 100
        pic = ['0' num2str(picnum)];
    else
        pic = [num2str(picnum)];
    end
    filename = ['cal_p' pic];
    fid = fopen([filename '.m'],'wt+');
    fprintf(fid,'function x = %s\n',filename);
    mat2text(x,fid);
    fclose(fid);
    
end