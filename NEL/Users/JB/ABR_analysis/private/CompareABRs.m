global abr_root_dir abr_data_dir

mFilesDir=fullfile(abr_root_dir,'ABR_analysis','private');
Datadir=fullfile(abr_data_dir,'Summary');

ChinName=input('Please enter Chin Name: ');

filename=sprintf('chin%d*.mat',ChinName);
cd(Datadir);
ChinFiles=dir(filename);
%There will always be multiple files with a given chin Name
%From each file we need Threshold and ABR amplitudes
months={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
freq_vect=[500,1000,2000,4000,8000];
freq_vect=freq_vect';
Date=zeros(1,length(ChinFiles));
ABR_thr=NaN*ones(5,length(ChinFiles));
Delta_ABR_thr=NaN*ones(5,length(ChinFiles));
for FileIND=1:length(ChinFiles)
	DataFile=ChinFiles(FileIND,1).name;
	ABRdateText=ChinFiles(FileIND,1).name(11:15); %to be used later more eff.
	Day=str2double(ABRdateText(1:2)); Month=find(ismember(months,ABRdateText(3:5)));
	if (Day<10)
		DayText=strcat('0',num2str(Day));
	else
		DayText=num2str(Day);
	end
	if (Month<10)
		MonthText=strcat('0',num2str(Month));
	else
		MonthText=num2str(Month);
	end
	DateText=strcat(MonthText,DayText);
	Date(1,FileIND)=str2double(DateText);
	%Date(1,FileIND)=(Month*1000)+str2double(Day);
	cd(Datadir);
	load(DataFile);
	Thrs=abrs.thresholds;
	for FreqIND=1:size(Thrs,1)
		Findex=find(ismember(freq_vect,Thrs(FreqIND,1)));
		ABR_thr(Findex,FileIND)=Thrs(FreqIND,2);
		%Now find the difference in ABR amplitudes
		ampINDs=find(abrs.y(:,1)==Thrs(FreqIND,1));
		ABRamp(:,2)=abs(abrs.y(ampINDs,4)-abrs.y(ampINDs,3));
		ABRamp(:,1)=abrs.y(ampINDs,2);
		ABRamplitude{FileIND,Findex}=ABRamp;
		clear ABRamp;
	end
	% 	for FreqIND=1:length(freq_vect)
	% 		Findex=find(Thrs(:,1)==freq_vect(FreqIND));
	% 		ABR_thr(FreqIND,FileIND)=Thrs(Findex,2);
	% 	end
end
minDateIND=find(Date==min(Date));
for FileIND=1:length(ChinFiles)
	Delta_ABR_thr(:,FileIND)=ABR_thr(:,FileIND)-ABR_thr(:,minDateIND);
	neg_inds=find(Delta_ABR_thr(:,FileIND)<0);
	Delta_ABR_thr(neg_inds,FileIND)=0;
end
figure(99); clf;
plot(freq_vect,Delta_ABR_thr,'-*')
ylabel('ABRthrReNormal'); xlabel('Freq (Hz)');
xlim([0 10000]);
legend(num2str(Date'),2)

cd(mFilesDir);