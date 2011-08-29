function data1

global abr_data_dir dd han_abrcomp colorr numfiles
% symbolandcolor=[...
% 	'\color[rgb]{0 0 0}* ';...
% 	'\color[rgb]{1 0 0}* ';...
% 	'\color[rgb]{0 0 1}* ';...
% 	'\color[rgb]{1 0 1}* ';...
% 	'\color[rgb]{0 1 0}* ';...
% 	'\color[rgb]{0 0 0}o ';...
% 	'\color[rgb]{1 0 0}o ';...
% 	'\color[rgb]{0 0 1}o ';...
% 	'\color[rgb]{1 0 1}o ';...
% 	'\color[rgb]{0 1 0}o ';...
% 	'\color[rgb]{0 0 0}x ';...
% 	'\color[rgb]{1 0 0}x ';...
% 	'\color[rgb]{0 0 1}x ';...
% 	'\color[rgb]{1 0 1}x ';...
% 	'\color[rgb]{0 1 0}x '];

d=dir([abr_data_dir 'Summary\']);
d = d(find(strncmp('.',{d.name},1)==0)); % Only files which are not '.' nor '..'
str = {d.name};
[selections ok] = listdlg('Name', 'File Manager', ...
	'PromptString',   'Select two datasets',...
	'SelectionMode',  'multiple',...
	'ListSize',       [300,300], ...
	'InitialValue',    1, ...
	'ListString',      str);
drawnow;

%clear 'titlename'

if ok ~= 0 & length(selections) >= 1 & length(selections) <=15
	numfiles=length(selections);
	for i=1:numfiles %%% Sort the files by date
		unsorted(i,1)=i; unsorted(i,2)=day_of(str{selections(1,i)});
	end
	[xx,ii]=sort(unsorted(:,2));
	sorted=unsorted(ii,:);
	
	
	for i=1:numfiles %Load the selected data files in order by date
		file(i,:)=str{selections(1,sorted(i))};
		dd(i).data=load([abr_data_dir 'Summary\' file(i,:)]);
		dd(i).name=file(i,:);
%		titlename(i,:)=[file(i,1:8) ' on '  file(i,13:15) ' ' file(i,11:12) ', ' file(i,16:19)];
		
		if isfield(dd(i).data,'abrs')
			
			dd(i).data.abrs.y(:,11)=dd(i).data.abrs.y(:,3)-dd(i).data.abrs.y(:,4); %Compute Wave I amplitude
			
			if day_of(dd(i).name) > 733818 %i.e. if after June18 2010
				dd(i).data.abrs.x(:,3)=dd(i).data.abrs.x(:,3) - 6.635; %Correction for delay of TDT + ER2
			elseif day_of(dd(i).name) <= 733818
				dd(i).data.abrs.x(:,3)=dd(i).data.abrs.x(:,3) - 6.635 + 0.147; %Correction for delay of TDT + free field speakers
			end
			
			[qq,ww]=sort(dd(i).data.abrs.thresholds(:,1)); %sort ABR thresholds by frequency
			dd(i).data.abrs.thresholds=dd(i).data.abrs.thresholds(ww,:);
			for ii=1:length(dd(i).data.abrs.y(:,1)) %Compute Sensation levels
				
				if dd(i).data.abrs.y(ii,1)==8000
					dd(i).data.abrs.y(ii,12)=dd(i).data.abrs.y(ii,2)-dd(i).data.abrs.thresholds(dd(i).data.abrs.thresholds(:,1)==8000,2);
				elseif dd(i).data.abrs.y(ii,1)==4000
					dd(i).data.abrs.y(ii,12)=dd(i).data.abrs.y(ii,2)-dd(i).data.abrs.thresholds(dd(i).data.abrs.thresholds(:,1)==4000,2);
				elseif dd(i).data.abrs.y(ii,1)==2000
					dd(i).data.abrs.y(ii,12)=dd(i).data.abrs.y(ii,2)-dd(i).data.abrs.thresholds(dd(i).data.abrs.thresholds(:,1)==2000,2);
				elseif dd(i).data.abrs.y(ii,1)==1000
					dd(i).data.abrs.y(ii,12)=dd(i).data.abrs.y(ii,2)-dd(i).data.abrs.thresholds(dd(i).data.abrs.thresholds(:,1)==1000,2);
				elseif dd(i).data.abrs.y(ii,1)==500
					dd(i).data.abrs.y(ii,12)=dd(i).data.abrs.y(ii,2)-dd(i).data.abrs.thresholds(dd(i).data.abrs.thresholds(:,1)==500,2);
					
				end
			end
			
			
		elseif ~isfield(dd(i).data,'abrs') %%% Dummy matrices for plots
			dd(i).data.abrs.thresholds=NaN*ones(1,4);
			dd(i).data.abrs.x=NaN*ones(1,10);
			dd(i).data.abrs.y=NaN*ones(1,12);
			
		end
		
		if ~isfield(dd(i).data,'dpoae') %%% Dummy matrix for plot
			dd(i).data.dpoae.data=NaN*ones(1,4);
			
		end
		
	end
	if length(selections) < 15
		for i=length(selections)+1:15
			dd(i).data.abrs.thresholds=NaN*ones(1,4);
			dd(i).data.abrs.x=NaN*ones(1,10);
			dd(i).data.abrs.y=NaN*ones(1,12);
			dd(i).data.dpoae.data=NaN*ones(1,4);
			dd(i).name='                       ';
%			titlename(i,:)='                        ';
		end
	end
	
	colorr=['k';'r';'b';'m';'g';'k';'r';'b';'m';'g';'k';'r';'b';'m';'g'];
	

	% 	if length(titlename(:,1)) <=5
	% 		text(0.05,0.98,titlename,'FontSize',14,'HorizontalAlignment','left','VerticalAlignment','top');
	% 	elseif length(titlename(:,1)) >5 & length(titlename(:,1)) <= 10
	% 		text(0.05,0.98,titlename(1:5,:),'FontSize',14,'HorizontalAlignment','left','VerticalAlignment','top');
	% 		text(0.28,0.98,titlename(6:end,:),'FontSize',14,'HorizontalAlignment','left','VerticalAlignment','top');
	% 	elseif length(titlename(:,1)) >10 & length(titlename(:,1)) <= 15
	% 		text(0.05,0.98,titlename(1:5,:),'FontSize',14,'HorizontalAlignment','left','VerticalAlignment','top');
	% 		text(0.28,0.98,titlename(6:10,:),'FontSize',14,'HorizontalAlignment','left','VerticalAlignment','top');
	% 		text(0.51,0.98,titlename(11:end,:),'FontSize',14,'HorizontalAlignment','left','VerticalAlignment','top');
	% 	end
	% 	for i=6:10
	% 		text(0.5,1.1-(i-5)*0.2,['* ' titlename(i,:)],'Fontsize',14,'HorizontalAlignment','Center','VerticalAlignment','Middle')
	% 	end
	% 	for i=11:15
	% 			text(5/6,1.1-(i-10)*0.2,['* ' titlename(i,:)],'Fontsize',14,'HorizontalAlignment','Center','VerticalAlignment','Middle')
	% 	end
	comp_plot
	
elseif length(selections)>15
	msgbox('15 files max')
end

