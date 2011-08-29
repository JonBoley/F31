%ABR Analysis Instruction Block

get_analysis_ins_bak;
copyfile(fullfile(root_dir,'ABR_analysis','private','get_analysis_ins_bak.m'), ...
    fullfile(root_dir,'ABR_analysis','private','get_analysis_ins.m'),'f');


errstr = sprintf('%s\n%s','Corrupted instructions have been repaired.','Check your parameters.');
uiwait(errordlg(errstr,'Analysis Error'));

