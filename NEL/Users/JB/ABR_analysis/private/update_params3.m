function update_params3;	
global abr_Stimuli abr_root_dir

filename = fullfile(abr_root_dir,'ABR_analysis','private','get_analysis_ins3.m');
fid = fopen(filename,'wt+');

fprintf(fid,'%s\n\n','%ABR Analysis Instruction Block');
fprintf(fid,'%s%d%s\n','abr_Stimuli = struct(''cal_pic'',''',str2num(abr_Stimuli.cal_pic),''', ...');
fprintf(fid,'\t%s%s%s\n','''abr_pic'',''',abr_Stimuli.abr_pic,''', ...');
fprintf(fid,'\t%s%5.2f%s\n','''start'',',abr_Stimuli.start,', ...');
fprintf(fid,'\t%s%5.2f%s\n','''end'',',abr_Stimuli.end,', ...');
fprintf(fid,'\t%s%5.2f%s\n','''start_template'',',abr_Stimuli.start_template,', ...');
fprintf(fid,'\t%s%5.2f%s\n','''end_template'',',abr_Stimuli.end_template,', ...');
fprintf(fid,'\t%s%5.2f%s\n','''num_templates'',',abr_Stimuli.num_templates,', ...');
fprintf(fid,'\t%s%s%s\n','''dir'',''',abr_Stimuli.dir,''');');

fclose(fid);
