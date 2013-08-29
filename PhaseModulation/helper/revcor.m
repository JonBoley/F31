function [H1spec,freqs] = revcor(spikeTrains,stimWav,sFreq)

spikeTrain = cell2mat(spikeTrains')';
Nreps = numel(spikeTrains);

stimWav=resample(stimWav,1,2);
stimPower=mean(stimWav.^2);
stimRMS=sqrt(stimPower);
stimPSD=stimPower/16500;

timeVector=[0:length(stimWav)-1]/sFreq;
endOnset=0.030;
stimEnd=max(timeVector);

spikeTrain=spikeTrain(1,spikeTrain(1,:)>endOnset);
spikeTrain=spikeTrain(1,spikeTrain(1,:)<stimEnd);

numSpikes=length(spikeTrain);
H0=numSpikes/(stimEnd-endOnset)/Nreps;

Npoints=2048;
delay=[0:Npoints-1]/sFreq;
R1=zeros(Npoints,1);
R1null=R1;
R2=zeros(Npoints,Npoints);
stimAutoCorrMx=R2;
waveXcorr=xcorr(stimWav,Npoints-1,'unbiased');
for ii=1:Npoints
	stimAutoCorrMx(:,ii)=waveXcorr(Npoints-ii+1:end-ii+1);
end


for i=1:numSpikes
	
	t=spikeTrain(i);
	bins=[round(t*sFreq):-1:round(t*sFreq)-Npoints+1];
    bins = max(1,bins);
		
	R1=R1+stimWav(bins);  %revcor - avg of waveform preceeding each spike
	
	R2=R2+stimWav(bins)*stimWav(bins)'-stimAutoCorrMx; %autocorrelation of waveform preceeding each spike
	
	tRan=random('Uniform',endOnset,stimEnd);
	bins=[round(tRan*sFreq):-1:round(tRan*sFreq)-Npoints+1];
    bins = max(1,bins);
		
	R1null=R1null+stimWav(bins);  %NF? - stim waveform preceeding random time to estimate noise floor?
	
end


R1=R1/numSpikes;
H1=R1*stimRMS/stimPower;

R2=R2/numSpikes;
H2=R2*stimPower/(2*stimPower^2);

R1null=R1null/numSpikes;
H1null=R1null*stimRMS/stimPower;


% figure(100); plot(delay, H1);
% figure(101); imagesc(delay, delay,H2);


freqs = (0:Npoints-1)/(Npoints-1)*sFreq;
H1spec = fft(H1);
H1spec = H1spec(2:floor(end/2));
freqs = freqs(2:floor(end/2));
if nargout<1
    figure, semilogx(freqs,abs(H1spec),'b');
    hold on; semilogx(freqs,unwrap(angle(H1spec)),'k:'); hold off;
end

