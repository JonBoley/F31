function save_dps

global dp_data_dir han_dp z zanimal zdate

default = strcat('chin',zanimal,'on',zdate);
filename = inputdlg('File name:','Save Data',1,default);


if ~isempty(filename)
	filename2 = char(strcat('C:\NEL\ExpData\Summary\',filename,'.mat'));
	dpoae.data=z.DpoaeData;

	if day_of(filename2) < 733976
		dpoae.data(:,4)=dpoae.data(:,4)-15; %Correction for flaw in early data (pre Nov 22, 2010)
	end
	
   if exist(filename2,'file')
	save(filename2,'dpoae','-append')
else
   save(filename2,'dpoae')
end
	clear dpoae;

	
end;

