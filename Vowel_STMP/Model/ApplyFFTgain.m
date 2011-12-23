function AmplifiedSignal=ApplyFFTgain(input,Fs,freqs,gains)

SpecOrig = fft(input);
MagSpecOrig = abs(SpecOrig);

f = (1:length(MagSpecOrig)/2)'/(Fs/2)*Fs;
maxLevel = max(MagSpecOrig(2:end/2));


FFTgain = interp1(freqs,gains,f,'linear','extrap');
FFTgain(FFTgain<1)=1;

% % Ramp gain back down at high freqs
% index = find(f>5e3,1,'first'); % ramp down after 5kHz
% n = length(FFTgain)-index;
% FFTgain(end-n+1:end)=FFTgain(end-n+1:end)-(1:n)'/n.*FFTgain(end-n+1:end);

MagSpecNew = MagSpecOrig.*[FFTgain; flipud(FFTgain)];
AmplifiedSignal = real(ifft(MagSpecNew.*exp(angle(SpecOrig))));


% % plot results
% figure, subplot(2,1,1),
% semilogx(f(2:end),20*log10(MagSpecNew(2:end/2)/maxLevel),'r');
% axis([20 Fs/2 -80 max(gains)]);
% hold on;
% semilogx(f(2:end),20*log10(MagSpecOrig(2:end/2)/maxLevel),'b');
% semilogx(f,FFTgain,'r:');
% 
% OrigSignal = real(ifft(MagSpecOrig.*exp(angle(SpecOrig))));
% subplot(2,1,2), plot(AmplifiedSignal,'r'); hold on;
% plot(OrigSignal,'b'); hold off;


