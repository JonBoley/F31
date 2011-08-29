function abrcomp(command_str)

global abr_data_dir han_abrcomp dd colorr numfiles

%

abr_data_dir = 'C:\NEL\ExpData\'; % added by GE 04Mar2004.

if nargin < 1
	
	
	abrcomp_FIG.handle = figure('NumberTitle','off','Name','ABR/DP comparison','Units','normalized','Visible','on',...
		'position',[0 0.03 1 0.92]);
	
	colordef none;
	whitebg('w');
	
	abrcomp_FIG.push.data1 = uicontrol(abrcomp_FIG.handle,'callback','abrcomp(''data1'');','style','pushbutton',...
		'Units','normalized','position',[0.83 0.95 0.11 .03],'string','Load datasets');
	abrcomp_FIG.push.print = uicontrol(abrcomp_FIG.handle,'callback','abrcomp(''print'');','style','pushbutton',...
		'Units','normalized','position',[0.83 0.905 0.11 .03],'string','Print');
	abrcomp_FIG.push.tracecolors = uicontrol(abrcomp_FIG.handle,'callback','abrcomp(''tracecolors'');','style','pushbutton',...
		'Units','normalized','position',[0.83 0.86 0.11 .03],'string','Change colors');

	han_abrcomp.latlevel = uicontrol(abrcomp_FIG.handle,'callback','abrcomp(''latlevel'');','style','checkbox','Units','normalized',...
		'position',[0.9 0.43 0.05 .03],'String','SL','Value',1,'BackgroundColor','w');
	han_abrcomp.amplevel = uicontrol(abrcomp_FIG.handle,'callback','abrcomp(''amplevel'');','style','checkbox','Units','normalized',...
		'position',[0.9 0.03 0.05 .03],'String','SL','Value',0,'BackgroundColor','w');
	
	axes('Position',[0.05 0.85 0.75 0.14]); axis off; han_abrcomp.text=gca;
	
	axes('Position',[0.075 0.48 0.2 0.33]); han_abrcomp.thresh=gca;
	
	axes('Position',[0.075 0.08 0.2 0.33]); han_abrcomp.dpoae=gca;
	
	axes('Position',[0.83 0.48 0.11 0.33]); han_abrcomp.lat8k=gca;
	axes('Position',[0.71 0.48 0.11 0.33]); han_abrcomp.lat4k=gca;
	axes('Position',[0.59 0.48 0.11 0.33]); han_abrcomp.lat2k=gca;
	axes('Position',[0.47 0.48 0.11 0.33]); han_abrcomp.lat1k=gca;
	axes('Position',[0.35 0.48 0.11 0.33]); han_abrcomp.lat05k=gca;
	
	axes('Position',[0.83 0.08 0.11 0.33]); han_abrcomp.amp8k=gca;
	axes('Position',[0.71 0.08 0.11 0.33]); han_abrcomp.amp4k=gca;
	axes('Position',[0.59 0.08 0.11 0.33]); han_abrcomp.amp2k=gca;
	axes('Position',[0.47 0.08 0.11 0.33]); han_abrcomp.amp1k=gca;
	axes('Position',[0.35 0.08 0.11 0.33]); han_abrcomp.amp05k=gca;
	
	
elseif strcmp(command_str,'data1')
	clear 'dd'
	data1
elseif strcmp(command_str,'print')
	set(gcf,'PaperOrientation','Landscape','PaperPosition',[0 0 11 8.5]);
	if ispc
		print('-dwinc','-r300');
	else
		print('-PNeptune','-dpsc','-r200','-noui');
	end
elseif strcmp(command_str,'tracecolors') & strcmp(get(han_abrcomp.thresh,'Box'),'on') %%%The box is on if data are loaded
	tracecolor
	
	
elseif strcmp(command_str,'latlevel') & strcmp(get(han_abrcomp.thresh,'Box'),'on') %%%The box is on if data are loaded
	comp_plot
elseif strcmp(command_str,'amplevel') & strcmp(get(han_abrcomp.thresh,'Box'),'on')
	comp_plot
end
