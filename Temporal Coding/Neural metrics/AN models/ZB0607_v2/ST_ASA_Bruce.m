function [STAplus,STAminus,dur]=ST_ASA_Bruce(filename1,cf,nrep,stimdb);
%
%
%%% this function is to create auditory nerve spike trains from Zilany and Bruce model..
% Program structure: [SpikeTrainx]=spiketrain(filename,cf,nreps,dur);
% SpikeTrainx is a cell array with spike times. 
%filename is the input .wav file. Please note that the SR should be
%100000Hz. cf is the characteristic frequency in Hz; nreps is the number of
%repetitions over which the model has to be run; dur is the duration of the
%stimuli over which the spike train has to be computed

tic
% clear all;clc
% warning('off','MATLAB:dispatcher:InexactMatch')
% 

%%%%%%%%%%To read the Chimaera files%%%%%%%%%%%%%%
filename3=strcat(filename1,'.wav');
[x1,f1]=wavread(filename3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%input parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x1,f1]=wavread('BBN_B_100k.wav');
% nrep=5; stimdb=50;
% cf=750;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%to set the model and stimulus parameters%%%%%%%%%%%%%%%%%%%%%%
dur=length(x1)/f1;
Npt=round(dur*f1);
RFT=10/1000;Nrft=round(RFT*f1);
fsres=100000;
tdres=1/fsres;
spont=50; %Zilany model has been tesed only for SR = 50
cohc = 1; cihc=1;
nrepm = 10;            % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.5e-3; % binwidth in seconds;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Chk for model compatibility%%%%%%%%%%%%%%%%%%%%%%%%%
if f1~=fsres
    warning('Wav file incorrectly sampled; Please use a sampling frequency of 100K')
else
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%%To scale the stimuli in Pascal%%%%%%%%%%%%%%%%%%%
a=20*sqrt(2)*10^(-6)*10^(stimdb/20);
x1=((a*x1)/sqrt(mean(x1.^2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%To add 10 msec ramp%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1(1:Nrft) = x1(1:Nrft).*(0:(Nrft-1)).'/(Nrft-1); % add triangular ramps
x1(end-Nrft+1:end) = x1(end-Nrft+1:end).*((Nrft-1):-1:0).'/(Nrft-1); % add triangular ramps
acx1=-x1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%Break it down; STAplus%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth] ...
    = zbcatmodel(x1.',cf,nrepm,1/fsres,dur,cohc,cihc,spont);

[sptime1,nspikes1] = sgmodel([tdres, nrep],sout); 
pause(1)
[sptime2,nspikes2] = sgmodel([tdres, nrep],sout);
pause(1)
[sptime3,nspikes3] = sgmodel([tdres, nrep],sout);
pause(1)
[sptime4,nspikes4] = sgmodel([tdres, nrep],sout);
pause(1)

nrepTEMP=nrep*4;
sptime=[sptime1(1:nspikes1); sptime2(1:nspikes2); sptime3(1:nspikes3); sptime4(1:nspikes4)];
nspikes=length(sptime);
LASTspikeINDS=[find(diff(sptime)<0)' nspikes];
FIRSTspikeINDs=[1 LASTspikeINDS(1:end-1)+1];
STAplus=cell(1,nrepTEMP);
for i=1:nrepTEMP
    STAplus{i}=sptime(FIRSTspikeINDs(i):LASTspikeINDS(i));
end

% filename1='BBN_B_100k';
filename=sprintf('STAplus_%s.mat',filename1);
eval(['save ''' filename ''' STAplus']);
clear sptime1 sptime2 sptime3 sptime4 sout
clear nspikes1 nspikes2 nspikes3 nspikes4 %STAplus
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%%%%%%%%%%%%%%Break it down; STAminus%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth] ...
    = zbcatmodel(acx1.',cf,nrepm,1/fsres,dur,cohc,cihc,spont);

[sptime1,nspikes1] = sgmodel([tdres, nrep],sout); 
pause(1)
[sptime2,nspikes2] = sgmodel([tdres, nrep],sout);
pause(1)
[sptime3,nspikes3] = sgmodel([tdres, nrep],sout);
pause(1)
[sptime4,nspikes4] = sgmodel([tdres, nrep],sout);
pause(1)

nrepTEMP=nrep*4;
sptime=[sptime1(1:nspikes1); sptime2(1:nspikes2); sptime3(1:nspikes3); sptime4(1:nspikes4)];
nspikes=length(sptime);
LASTspikeINDS=[find(diff(sptime)<0)' nspikes];
FIRSTspikeINDs=[1 LASTspikeINDS(1:end-1)+1];
STAminus=cell(1,nrepTEMP);
for i=1:nrepTEMP
    STAminus{i}=sptime(FIRSTspikeINDs(i):LASTspikeINDS(i));
end

filename=sprintf('STAminus_%s.mat',filename1);
eval(['save ''' filename ''' STAminus']);
clear sptime1 sptime2 sptime3 sptime4 sout
clear nspikes1 nspikes2 nspikes3 nspikes4 %STAminus
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
