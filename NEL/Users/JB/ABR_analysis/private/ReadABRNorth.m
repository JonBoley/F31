fid = fopen([prefix num2str(picnum) '.dat']);
file_prog = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
spc = find(isspace(tline));
file_date = tline(spc(3)+1:length(tline));
tline = fgetl(fid);
tline = fgetl(fid);
file_stim_freq = str2num(tline(length(tline)-6:length(tline)));
tline = fgetl(fid);
file_stim_dur = str2num(tline(length(tline)-6:length(tline)));
tline = fgetl(fid);
file_rec_wind = str2num(tline(length(tline)-6:length(tline)));
tline = fgetl(fid);
file_samp_rate = 1e3/str2num(tline(length(tline)-6:length(tline)));
tline = fgetl(fid);
tline = fgetl(fid);
file_stim_bursts = str2num(tline(length(tline)-6:length(tline)));
tline = fgetl(fid);
file_rise_fall = str2num(tline(length(tline)-6:length(tline)));
tline = fgetl(fid);
file_num_samps = str2num(tline(length(tline)-6:length(tline)));
tline = fgetl(fid);
file_stim_atten = str2num(tline(length(tline)-6:length(tline)));
tline = fgetl(fid);
file_samp_gain = str2num(tline(length(tline)-6:length(tline)));
tline = fgetl(fid);
file_amp_vlt = NaN;
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);

line = 0;

while 1
    tline = fgetl(fid);
    if tline ~= -1
        line = line + 1;
        file_data(line,1) = str2num(tline(5:10));
        file_data(line,2) = str2num(tline(15:20));
    else
        break
    end
end