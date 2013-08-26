% File: GenVowels_1.m
% M. Heinz
% 18 May 2007


heed_forms =  [250 2250 3050 3350 3850];
hod_forms =   [750 1050 2950 3350 3850];
whod_forms =  [250  850 2250 3350 3850];
had_forms =   [750 1450 2450 3350 3850];
heard_forms = [450 1150 1250 3350 3850];

BWs= [ 90 110 170 250 300];
dur = 0.5;
F0s=[100 150 200];
Fs=24414.062500;

[time, heed_100_WAV] = dovowl(heed_forms,BWs,F0s(1),dur,Fs);
heed_100_WAV=heed_100_WAV./max(abs(heed_100_WAV))*0.99;
wavplay(heed_100_WAV,Fs);
wavwrite(heed_100_WAV,Fs,'heed_100.wav');

[time, heed_150_WAV] = dovowl(heed_forms,BWs,F0s(2),dur,Fs);
heed_150_WAV=heed_150_WAV./max(abs(heed_150_WAV))*0.99;
wavplay(heed_150_WAV,Fs);
wavwrite(heed_150_WAV,Fs,'heed_150.wav');


[time, hod_100_WAV] = dovowl(hod_forms,BWs,F0s(1),dur,Fs);
hod_100_WAV=hod_100_WAV./max(abs(hod_100_WAV))*0.99;
wavplay(hod_100_WAV,Fs);
wavwrite(hod_100_WAV,Fs,'hod_100.wav');

[time, hod_150_WAV] = dovowl(hod_forms,BWs,F0s(2),dur,Fs);
hod_150_WAV=hod_150_WAV./max(abs(hod_150_WAV))*0.99;
wavplay(hod_150_WAV,Fs);
wavwrite(hod_150_WAV,Fs,'hod_150.wav');


[time, whod_100_WAV] = dovowl(whod_forms,BWs,F0s(1),dur,Fs);
whod_100_WAV=whod_100_WAV./max(abs(whod_100_WAV))*0.99;
wavplay(whod_100_WAV,Fs);
wavwrite(whod_100_WAV,Fs,'whod_100.wav');

[time, whod_150_WAV] = dovowl(whod_forms,BWs,F0s(2),dur,Fs);
whod_150_WAV=whod_150_WAV./max(abs(whod_150_WAV))*0.99;
wavplay(whod_150_WAV,Fs);
wavwrite(whod_150_WAV,Fs,'whod_150.wav');


[time, had_100_WAV] = dovowl(had_forms,BWs,F0s(1),dur,Fs);
had_100_WAV=had_100_WAV./max(abs(had_100_WAV))*0.99;
wavplay(had_100_WAV,Fs);
wavwrite(had_100_WAV,Fs,'had_100.wav');

[time, had_150_WAV] = dovowl(had_forms,BWs,F0s(2),dur,Fs);
had_150_WAV=had_150_WAV./max(abs(had_150_WAV))*0.99;
wavplay(had_150_WAV,Fs);
wavwrite(had_150_WAV,Fs,'had_150.wav');


[time, heard_100_WAV] = dovowl(heard_forms,BWs,F0s(1),dur,Fs);
heard_100_WAV=heard_100_WAV./max(abs(heard_100_WAV))*0.99;
wavplay(heard_100_WAV,Fs);
wavwrite(heard_100_WAV,Fs,'heard_100.wav');

[time, heard_150_WAV] = dovowl(heard_forms,BWs,F0s(2),dur,Fs);
heard_150_WAV=heard_150_WAV./max(abs(heard_150_WAV))*0.99;
wavplay(heard_150_WAV,Fs);
wavwrite(heard_150_WAV,Fs,'heard_150.wav');



