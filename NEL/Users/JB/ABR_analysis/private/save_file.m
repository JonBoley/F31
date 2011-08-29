function save_file

global num freq spl animal date freq_level data abr ABRmag w
	
default = strcat('chin',animal,'on',date);
filename = inputdlg('File name:','Save Data',1,default);

if ~isempty(filename)
	filename2 = char(strcat('C:\NEL\ExpData\Summary\',filename,'.mat'));

	freq2=ones(1,num)*freq; replaced=0;

	if exist(filename2,'file')
		load(filename2, 'abrs')
		if ismember(freq,abrs.thresholds(:,1))
			abrs.thresholds(abrs.thresholds(:,1)==freq,:)=[];
			abrs.z.par(abrs.z.par(:,1)==freq,:)=[];
			abrs.z.score(abrs.z.score(:,1)==freq,:)=[];
			abrs.amp(abrs.amp(:,1)==freq,:)=[];
			abrs.x(abrs.x(:,1)==freq,:)=[];
			abrs.y(abrs.y(:,1)==freq,:)=[];
			abrs.waves(abrs.waves(:,1)==freq,:)=[];
			replaced=1;
		end;
		abrs.thresholds = [abrs.thresholds; freq data.threshold data.amp_thresh -freq_level];
		abrs.z.par = [abrs.z.par; freq data.z.intercept data.z.slope];
		abrs.z.score = [abrs.z.score; freq2' spl' data.z.score' w'];
		abrs.amp = [abrs.amp; freq2' ABRmag];
		abrs.x = [abrs.x; freq2' spl' data.x'];
		abrs.y = [abrs.y; freq2' spl' data.y'];
		abrs.waves = [abrs.waves; freq2' spl' abr'];
		save(filename2, 'abrs'); clear abrs;
		if replaced==0
			msgbox('New data added')
		else
			msgbox('Data replaced')
		end
	
	else
		abrs.thresholds = [freq data.threshold data.amp_thresh -freq_level];
		abrs.z.par = [freq data.z.intercept data.z.slope];
		abrs.z.score = [freq2' spl' data.z.score' w'];
		abrs.amp = [freq2' ABRmag];
		abrs.x = [freq2' spl' data.x'];
		abrs.y = [freq2' spl' data.y'];
		abrs.waves = [freq2' spl' abr'];
		save(filename2, 'abrs'); clear abrs;
		msgbox('New file created')
	end
end