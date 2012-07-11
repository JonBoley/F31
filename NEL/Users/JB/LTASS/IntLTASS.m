% Byrne, et al. (1994), "An international comparison of long-term average speech spectra," J. Acoust. Soc. Am. 96, 2108-2120.

function Z = IntLTASS(x,Fs)

band = [100,125,160,200,250,315,400,500,630,800,1000,...
    1250,1600,2000,2500,3150,4000,5000,6300,8000,10000];
noiseLevel = [54.4,57.7,56.8,60.2,60.3,59,62.1,62.1,60.5,56.8,53.7,...
    53,52,48.7,48.1,46.8,45.6,44.5,44.3,43.7,43.4];

noisePower = 10.^(noiseLevel./20);
nF = length(band);
N = 3; 					% Order of analysis filters.

npts = length(x);
z = zeros(1,npts);

ff = (1000).*((2^(1/3)).^[-10:nF-11]); 	% Exact center freq.

fnyq=0.5*Fs; %Nyquist frequency
edge=0.5*(ff(1:end-1) + ff(2:end)); %Band edges in Hz
tfir=8; %Length of the FIR filter impulse response in msec
nfir=round(0.001*tfir*Fs); %Length of the FIR filters in samples
nfir=2*floor(nfir/2); %Force filter length to be even
ft=175; %Half the width of the filter transition region

for i = 1:2
    % Set up storage for the output
    z=zeros(npts,nF); % Do not include filter transients
    
    % First band is a low-pass filter
    y=zeros(1,npts+nfir); % Include filter transients
    w=zeros(1,npts); %Do not Include filter transients
    gain=[1 1 0 0];
    f=[0,(edge(1)/1.1),(edge(1)*1.1),fnyq];
    b=fir2(nfir,f/fnyq,gain); %FIR filter design
    y=conv(x,b);
    w = y(nfir+1:end);
    z(:,1)=w.*noisePower(1)./rms2(w);
    
    % Last band is a high-pass filter
    y=zeros(1,npts+nfir); % Include filter transients
    w=zeros(1,npts); %Do not Include filter transients
    gain=[0 0 1 1];
    f=[0,(edge(nF-1)/1.1),(edge(nF-1)*1.1),fnyq];
    b=fir2(nfir,f/fnyq,gain); %FIR filter design
    y=conv(x,b);
    w = y(nfir+1:end);
    z(:,nF)=w.*noisePower(nF)./rms2(w);
    
    % Remaining bands are bandpass filters
    y=zeros(1,npts+nfir); % Include filter transients
    w=zeros(1,npts); %Do not Include filter transients
    gain=[0 0 1 1 0 0];
    for n=2:nF-1
        f=[0,(edge(n-1)/1.1),(edge(n-1)*1.1),(edge(n)/1.1),(edge(n)*1.1),fnyq];
        b=fir2(nfir,f/fnyq,gain); %FIR filter design
        y=conv(x,b);
        w = y(nfir+1:end);
        z(:,n)=w.*noisePower(n)./rms2(w);
    end
    
    Z = sum(z');
    
    [p,f] = oct3bank2(-99,-99,Z,Fs);
    p = p+(noiseLevel(1)-p(1));
    p(end+1) = noiseLevel(end);
    if length(p) < length(noiseLevel)
        noisePower = 10.^((noiseLevel+(noiseLevel-p))./20);
    end
end;