%% calculate neurograms
clear input;
[input,Fs]=wavread('heed_100.wav');
% Fs=16000; input=chirp((1:.5*Fs)/Fs,4000,.5,1000,'logarithmic')';
OALevel_dBSPL = 70;
lowCF=80; %Hz
highCF=8000; %Hz
numCFs=80;  % 80 = 1 semitone steps from 80Hz to 8kHz
spread=log2(highCF/lowCF); %number of octaves to span
CFs=(lowCF*2.^(0:spread/(numCFs-1):spread));

% Normal = neurogram(input,Fs,OALevel_dBSPL,lowCF,highCF,numCFs);
% Impaired = neurogram(input,Fs,OALevel_dBSPL,lowCF,highCF,numCFs,'mixed');
Aided = neurogram(input,Fs,OALevel_dBSPL,lowCF,highCF,numCFs,'mixed',true);
% save test

%% Phase-Correction
GD=0.050*Fs; % group delay (samples)
k0=(GD-2)/(GD+2); 
k1=-cos(2*pi*2250/Fs); 
k2=1; 
B=[k0 k1*(1+k0*k2) k2]; 
A=fliplr(B);
% figure(999), plot(1:Fs/2,grpdelay(B,A,Fs/2,Fs));
input2 = filter(B,A,input);
Corrected = neurogram(input2,Fs,OALevel_dBSPL,lowCF,highCF,numCFs,'mixed',true);
save test

%% plot neurograms
clim = [min(min(Normal)) max(max(Normal))];
% clim = [0 0.4*(max(max(Normal))-min(min(Normal)))] + min(min(Normal));
zoomarea = 20e3:35e3;
figure(1); colormap hot;%(flipud(gray(256))); %([1 1 1;0 0 0]);
subplot(4,1,1), imagesc(Normal(:,zoomarea),clim); axis off; title('Normal /i/');
subplot(4,1,2), imagesc(Impaired(:,zoomarea),clim); axis off; title('Impaired /i/');
subplot(4,1,3), imagesc(Aided(:,zoomarea),clim); axis off; title('Aided (NAL) /i/');
subplot(4,1,4), imagesc((zoomarea(1):zoomarea(end))/100e3,CFs/1000,Corrected(:,zoomarea),clim); 
set(gca,'YTickLabel',round(CFs(floor((1:4)*end/4))/10)/100);
title('Phase-Corrected (NAL) /i/'); xlabel('time (sec)'); ylabel('CF (kHz)');

%% plot neural "spectrum"
figure(2), semilogx(CFs,sum(Normal,2),'b'); hold on;
semilogx(CFs,sum(Impaired,2),'r'); hold on;
semilogx(CFs,sum(Aided,2),'g'); hold on;
legend({'Normal','Impaired','Aided'},'location','nw'); 
axis tight; title('Average Rate'); xlabel('Frequency (Hz)');

%% plot avg PSTH (averaged across CF)
figure(3), plot(sum(Normal,1),'b'); hold on;
plot(sum(Impaired,1),'r'); hold on;
plot(sum(Aided,1),'g'); hold on;
legend({'Normal','Impaired','Aided'},'location','nw'); 
axis tight; title('Average Rate'); xlabel('Time');

%% plot ACF
CorrSeq = zeros(size(Aided,1),size(xcorr(Aided(1,:),Aided(1,:)),2));
for i=1:size(Aided,1)
   CorrSeq(i,:) =  xcorr(Aided(i,:),Aided(i,:));
end
figure, imagesc(CorrSeq);
set(gca,'YTickLabel',round(CFs(floor((1:8)*end/8))/10)/100);
ylabel('CF (kHz)');
