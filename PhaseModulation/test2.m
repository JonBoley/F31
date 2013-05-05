%% test2
% apply AP filter to subset of harmonics
% do we hear 2 pitches?

f0 = 100; %Hz
numHarm = 12; % number of harmonics
Fs = 24e3; %sampling frequency
t=(1:Fs)/Fs;

input = zeros(1,length(t));
for i=1:numHarm
    input = input + sin(2*pi*i*f0*t)/numHarm;
end

%filters
fPeak = 3:3:numHarm;
numStages = 5;
Hcas=dfilt.scalar;
GD=0.001*Fs; % group delay (samples)
k0=(GD-2)/(GD+2);
k2=1;
for i=1:numel(fPeak)
    % note that phase delay (~0.2ms/section) is independent of group delay
    % also note that max phase delay is not at center freq
    k1=-cos(2*pi*fPeak(i)/Fs);
    B=[k0 k1*(1+k0*k2) k2];
    A=fliplr(B);
    % figure(999), plot(1:Fs/2,grpdelay(B,A,Fs/2,Fs));
    
    H(i)=dfilt.df2t(B,A);
    for j=1:numStages
        Hcas=dfilt.cascade(Hcas,H(i));
    end
end

output = filter(Hcas,input);

figure, plot(t,input,'b'); hold on;
plot(t,output,'r'); hold off;

rampdur = round(0.020*Fs);
win = [linspace(0,1,rampdur),...
    ones(1,Fs*max(t)-2*rampdur),...
    linspace(1,0,rampdur)];

soundsc(win.*input,Fs);
pause(max(t)+0.5);
soundsc(win.*output,Fs);
