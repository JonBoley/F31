% Generate LTASS noise
function Z=GenLTASS(sec,Fs,SaveWAV)
addpath octave;
Z = IntLTASS(randn(sec*Fs,1),Fs);
Z = Z/(max(abs(Z))+eps);

if 0 % plot FFT?
    BufferLen = 2048;
    Z2 = buffer(Z,BufferLen,BufferLen/2); % buffer w/ 50% overlap
    Z2 = repmat(hamming(BufferLen),1,size(Z2,2)).*Z2;
    
    LTASSspec = mean(abs(fft(Z2)),2);
    LTASSspec = LTASSspec/max(LTASSspec);
    freqs = (Fs/length(Z2)):(Fs/length(Z2)):Fs;
    semilogx(freqs(1:end/2),20*log10(LTASSspec(1:end/2)));
    
    % plot 1/3-octave spectrum
    [p,f] = oct3bank(Z,Fs);
    figure, semilogx(f,p); xlabel('Frequency (Hz)'); ylabel('dBFS');
end

if nargin>2 % save WAV file?
    [FileName,PathName,FilterIndex] = ...
        uiputfile('*.wav','Save WAV File...','LTASS_noise.wav');
    wavwrite(Z,Fs,16,[PathName FileName]);
end
