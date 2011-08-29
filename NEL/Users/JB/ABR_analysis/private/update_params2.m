function update_params2;	
global abr_Stimuli abr_root_dir

filename = fullfile(abr_root_dir,'ABR_analysis','private','get_analysis_ins2.m');
fid = fopen(filename,'wt+');

fprintf(fid,'%s\n\n','%ABR Analysis Instruction Block');
fprintf(fid,'%s%d%s\n','abr_Stimuli = struct(''cal_pic'',''',str2num(abr_Stimuli.cal_pic),''', ...');
fprintf(fid,'\t%s%s%s\n','''abr_pic'',''',abr_Stimuli.abr_pic,''', ...');
fprintf(fid,'\t%s%5.2f%s\n','''start_resp'',',abr_Stimuli.start_resp,', ...');
fprintf(fid,'\t%s%5.2f%s\n','''end_resp'',',abr_Stimuli.end_resp,', ...');
fprintf(fid,'\t%s%5.2f%s\n','''start_back'',',abr_Stimuli.start_back,', ...');
fprintf(fid,'\t%s%5.2f%s\n','''end_back'',',abr_Stimuli.end_back,', ...');
fprintf(fid,'\t%s%5.2f%s\n','''scale'',',abr_Stimuli.scale,', ...');
fprintf(fid,'\t%s%s%s\n','''dir'',''',abr_Stimuli.dir,''');');

fclose(fid);
