fid = fopen([prefix num2str(picnum) '.dat']);
tline = fgetl(fid);
spc = find(isspace(tline));
file_prog = tline(1:spc(1)-1);
file_date = tline(spc(4)+1:length(tline));
tline = fgetl(fid);
tline = fgetl(fid);
file_frqlo = str2num(tline(1:6));
file_frqhi = str2num(tline(10:15));
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
spc = find(isspace(tline));
file_isyslv = str2num(tline(spc(2):spc(2)+3));
tline = fgetl(fid);
tline = fgetl(fid);
spc = find(isspace(tline));
file_mic = tline(spc(4):length(tline));
file_amp_vlt = NaN;
tline = fgetl(fid);
tline = fgetl(fid);

line = 0;

while 1
    tline = fgetl(fid);
    if tline ~= -1
        line = line + 1;
        file_data(line,1) = str2num(tline(2:10));
        file_data(line,2) = str2num(tline(12:20));
    else
        break
    end
end