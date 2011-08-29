function refdata

global abr_data_dir freq reff

d=dir([abr_data_dir 'Summary\']);
d = d(find(strncmp('.',{d.name},1)==0)); % Only files which are not '.' nor '..'
str = {d.name};
[selection ok] = listdlg('Name', 'File Manager', ...
    'PromptString',   'Select a reference dataset',...
    'SelectionMode',  'single',...
    'ListSize',       [300,300], ...
    'InitialValue',    1, ...
    'ListString',      str);
drawnow;

if (ok==0 | isempty(selection))
else
	reffile=str{selection};
	reff = load([abr_data_dir 'Summary\' reffile], 'abrs');
	
	if ismember(freq,reff.abrs.thresholds(:,1))
 		reff.abrs.thresholds(reff.abrs.thresholds(:,1)~=freq,:)=[];    
		reff.abrs.z.par(reff.abrs.z.par(:,1)~=freq,:)=[];              
		reff.abrs.z.score(reff.abrs.z.score(:,1)~=freq,:)=[];        
		reff.abrs.amp(reff.abrs.amp(:,1)~=freq,:)=[];                  
		reff.abrs.x(reff.abrs.x(:,1)~=freq,:)=[];                      
		reff.abrs.y(reff.abrs.y(:,1)~=freq,:)=[];                     
		reff.abrs.waves(reff.abrs.waves(:,1)~=freq,:)=[];              
		
		plot_data
		
	else
		msgbox('Reference data does not exist')
	end;
	
end;