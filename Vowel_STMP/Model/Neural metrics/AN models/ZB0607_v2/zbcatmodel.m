% zbcatmodel - Zilany and Bruce (JASA 2006, 2007) Auditory Nerve Model
%
% [timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,synout,psth] ...
%   = zbcatmodel(pin,CF,nrep,binwidth,stimtime,cohc,cihc,spont);
%
% timeout is an array of times in seconds
% meout is the output of the middle-ear filter
% c1filterout is the output of the C1 (signal path) BM filter
% c2filterout is the output of the C2 (parallel path) BM filter
% c1vihc is the output of the C1 IHC transduction function
% c2vihc is the output of the C2 IHC transduction function
% vihc is the IHC potential
% synout is the synapse output in spikes/s
% psth is the peri-stimulus time histogram
%
% pin is the input sound wave in Pa sampled at the appropriate sampling rate (see instructions below)
% cf is the characteristic frequency of the fiber in Hz
% nrep is the number of repetitions for the psth
% binwidth is the binsize in seconds, i.e., the reciprocal of the sampling rate (see instructions below)
% reptime is the time between stimulus repetitions in seconds - NOTE should be equal to or longer than the duration of pin
% cohc is the ohc scaling factor: 1 is normal OHC function; 0 is complete OHC dysfunction
% cihc is the ihc scaling factor: 1 is normal IHC function; 0 is complete IHC dysfunction
% spont is the spontaneous rate of the fiber in spikes/s - NOTE a value of 50 was used in Zilany and Bruce (2006)
% 
% Fore example,
%
%   [timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,synout,psth] ...
%    = zbcatmodel(pin,1e3,10,1/500e3,0.200,1,1,50);
%
% models a normal fiber of spontaneous rate 50 spikes/sec (normal OHC & IHC function) with a CF of 1 kHz, 
% for 10 repititions and a sampling rate of 500kHz, for a repetition duration of 200 ms.
%
% NOTE ON SAMPLING RATE:-
% Since version 2 of the code, it is possible to run the model at a range
% of sampling rates between 100 kHz and 500 kHz.
% It is recommended to run the model at 100 kHz for CFs up to 20 kHz, and
% at 200 kHz for CFs up to 40 kHz.