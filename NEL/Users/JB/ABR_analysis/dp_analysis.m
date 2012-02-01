function dp_analysis(command_str)

global dp_data_dir han_dp z zanimal zdate

%
switch lower(getenv('computername'))
    case {'m4500'}
        dp_data_dir = 'C:\Research\MATLAB\Vowel_STMP\ExpData\'; 
    otherwise
        dp_data_dir = 'C:\NEL\ExpData\'; 
end

if nargin < 1
	
	
	dp_FIG.handle = figure('NumberTitle','off','Name','DP analysis','Units','normalized','Visible','on',...
		'position',[0 0.03 1 0.92]);
	
	colordef none;
	whitebg('w');
	
	dp_FIG.push.loaddp = uicontrol(dp_FIG.handle,'callback','dp_analysis(''loaddp'');','style','pushbutton',...
		'Units','normalized','position',[0.05 0.05 0.1 .03],'string','Load DP data');
	dp_FIG.change_spectrum = uicontrol(dp_FIG.handle,'callback','dp_analysis(''change_spectrum'');','style','pushbutton',...
		'Units','normalized','position',[0.07 0.86 0.07 0.02],'string','Edit DPs');
	dp_FIG.push.save_dps = uicontrol(dp_FIG.handle,'callback','dp_analysis(''save_dps'');','style','pushbutton',...
		'Units','normalized','position',[0.18 0.05 0.1 .03],'string','Save DP data');
	dp_FIG.push.print_dps = uicontrol(dp_FIG.handle,'callback','dp_analysis(''print_dps'');','style','pushbutton',...
		'Units','normalized','position',[0.31 0.05 0.1 .03],'string','Print');
	
	axes('Position',[0 0 1 1]); han_dp.text=gca; axis off;
	
	axes('Position',[0.075 0.20 0.8 0.65]); han_dp.dp_curve=gca;
	
elseif strcmp(command_str,'loaddp')
	loaddp
elseif strcmp(command_str,'change_spectrum') & strcmp(get(han_dp.dp_curve,'Box'),'on')
	change_spectrum
elseif strcmp(command_str,'save_dps') & strcmp(get(han_dp.dp_curve,'Box'),'on')
	save_dps
elseif strcmp(command_str,'print_dps') & strcmp(get(han_dp.dp_curve,'Box'),'on')
    set(gcf,'PaperOrientation','Landscape','PaperPosition',[0 0 11 8.5]);
    if ispc
        print('-dwinc','-r200');
    else
        print('-PNeptune','-dpsc','-r200','-noui');
    end
end