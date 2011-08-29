function loaddp

global dp_data_dir han_dp z zanimal zdate

%clear 'global' 'd1' 'd2' 'd3'

d=dir(dp_data_dir);
d = d(find(strncmp('.',{d.name},1)==0)); % Only files which are not '.' nor '..'
str = {d.name};
[selection ok] = listdlg('Name', 'File Manager', ...
	'PromptString',   'Select data directory',...
	'SelectionMode',  'single',...
	'ListSize',       [300,300], ...
	'InitialValue',    1, ...
	'ListString',      str);
drawnow;
dir1=str{selection(1,1)};


cd([dp_data_dir dir1])
dp_pics=dir('*dpoae*');
str = {dp_pics.name};
[selection ok] = listdlg('Name', 'File Manager', ...
	'PromptString',   'Select dp file',...
	'SelectionMode',  'single',...
	'ListSize',       [300,300], ...
	'InitialValue',    1, ...
	'ListString',      str);
drawnow;

if ok ~= 0 & length(selection)==1
   currentDirectory=cd;
   if strcmp(currentDirectory(30:33),'chin') | strcmp(currentDirectory(30:33),'Chin')
         suggestedID={currentDirectory(34:37)};
   else
         suggestedID={['']};
   end
   zanimal = inputdlg('Animal ID number:','',1,suggestedID);
 	dp_pic=str{selection(1,1)};
	dp_pic=dp_pic(1:end-2);
	eval(['z=' dp_pic ';'])
	%zdate=regexprep(z.General.date,'-','');
	zdate=[z.General.date(1:2) z.General.date(4:6) z.General.date(8:11)];
   z.DpoaeData(:,3)=2*z.DpoaeData(:,1)-z.DpoaeData(:,2);
	axes(han_dp.dp_curve); plot(z.DpoaeData(:,3),z.DpoaeData(:,4),'*k-')
	title(['Chin' char(zanimal) ' on ' char(zdate)],'FontSize',14)
	ylabel('DP amplitude (dB SPL)','fontsize',14)
	xlabel('DP frequency (Hz)','fontsize',14)
	%set(han_dp.dp_curve,'Box','on','XGrid','on','YGrid','on','Xscale','log','Xlim',[300 12000],...
	%	'XTick',[300 400 500 600 800 1000 2000 3000 4000 5000 6000 8000 10000],'YLim',[0 100])
   set(han_dp.dp_curve,'Box','on','XGrid','on','YGrid','on','Xscale','log','Xlim',[300 10000],...
      'YLim',[0 100])
	
end
