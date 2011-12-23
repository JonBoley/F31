% GenEh.m
% From Miller '97:
% F1(0.5/0.2 kHz)
% F2(1.7/0.2 kHz)
% F3(2.5/0.2 kHz)
% F4(3.3/0.25 kHz)
% F5(3.7/0.2 kHz)
formants =  [500 1700 2500 3300 3700];
BWs= [200 200 200 250 200];
FeaturesText={'F1','F2','F3','F4','F5'};

dur = 0.5;
F0=100;
Fs=24414.062500;

[time, vowel] = dovowl(formants,BWs,F0,dur,Fs);
vowel=vowel./max(abs(vowel))*0.99; % normalize
