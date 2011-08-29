%function manual_NI(STIMind)
% File: manual_NI.m
%
% Written by MH 24oct2003.
%
% to debug NI board by hand

global signals_dir
global root_dir % Used to find correct control 'rco' file.
global RP

if ~exist('STIMind','var')
   STIMind=1;
end

uRate=80000;
atten=40;
Gating.Duration=300;
Gating.Period=1000;
MAXtrigs=1;

rconame = [root_dir 'stimulate\object\control.rco'];
rc = RPload_rco(rconame);   %Creates click: MH 10/20/03
% use L3 as source for Chan0data:
dev_description = nel_devices_vector('L3','list');
if (rc == 0)
   return;
end

% Load the waveform from the file and place in 'Chan0data':
Nstim=10;
test_tone_list=cell(1,Nstim);
for ifiles = 1:Nstim
%   test_tone_list{ifiles} = sprintf('%sMH\\test_tone%d.wav', signals_dir, ifiles);
%   test_tone_list{ifiles} = sprintf('%sGE\\RSS_stimA_set01\\RSS_%03d.wav', signals_dir, ifiles);
   test_tone_list{ifiles} = sprintf('%sGE\\RSS_stimA_set01\\RSS_003.wav', signals_dir);
end

[data,sr,rc] = nel_wavread(test_tone_list{STIMind});
if (rc == 0)
   nelerror(['Can''t read wavfile ''' static_bi.list{min(ind,end)} '''']);
else
   Chan0data = data;
end
% Rescape to 5 volts
Chan0data=Chan0data*5.0;

[statCodes retVals ireturnWAV] = d2a(1, uRate, Chan0data);  % Loads data to NI board
%%% ERROR check return codes from d2a.c
lineNum=56;
rc_d2a = 0;
if ( (~isempty(find(statCodes~=0))) | (~isempty(find(retVals~=0))) )
   rc_d2a = -1;
   disp(sprintf('\tManual_NI Location: %d',lineNum))
   disp([sprintf('\t\t') 'StatCodes: ' sprintf('%d ',statCodes')])
   disp([sprintf('\t\t') 'Return Vals: ' sprintf('%d ',retVals')])
   nelwarn(['non-zero status/return codes in ''DALinloop_NI_wavfiles'' ']);
end
%%%

%%%% Setting up switch box: select/connect/Attens
left_dev =  atten * nel_devices_vector('L3');  % Using L3 as device for Chan0data.
right_dev = nel_devices_vector([]);
attens_devices = [left_dev right_dev];     
[select,connect,PAattns] = find_mix_settings(attens_devices);
if (isempty(select) | isempty(connect))
   nelerror('''data_acquisition_loop'': Can''t find appropriate select and connect parameters. Aborting...');
   return;
end
rc = SBset(select,connect);
rc = PAset(PAattns) & (rc==1);

%%%% Setting up triggering on RP2
rc = (Trigcheck & (rc==1));
if (rc)
   rc = TRIGset(Gating.Duration,Gating.Period,MAXtrigs) & (rc==1);
   trig_off_time = (Gating.Period - Gating.Duration) / 1000;
   trig_period   = Gating.Period / 1000;
   [rc_set,RP] = RPset_params(RP);
   rc = rc_set & (rc==1);
end
if rc~=1
   nelerror('''STM'': Error(s) detected within manual_NI setup');
end


%% Setup to re-create BUG, by performing "find" MATLAB function right after starting trigger.
%% Seems to be a memory allocation issue 

%temp=input('Press Enter to start');
clear dummy_vec dummy_ans

TRIGstart;
dummy_vec=round(rand(1,10000)*10);  
%%% *** MH: 11/3/03: 10000 cause bug (most of the time, if not, switch 10000 to 1000000 and back!), 
% 1000000 does not cause bug???
dummy_ans=find(dummy_vec==3);

pause(Gating.Period*MAXtrigs/1000+.5)  % Wait for trigs to run out on RP2


%Cleanup TDT & NI systems to end.
[statCodes retVals] = d2a(2);  % Terminate/Clean-up NI 6052 board (see file "d2a.c")
lineNum=125;
rc_d2a = 0;
if ( (~isempty(find(statCodes~=0))) | (~isempty(find(retVals~=0))) )
   rc_d2a = -1;
   disp(sprintf('\tManual_NI Location: %d',lineNum))
   disp([sprintf('\t\t') 'StatCodes: ' sprintf('%d ',statCodes')])
   disp([sprintf('\t\t') 'Return Vals: ' sprintf('%d ',retVals')])
   nelwarn(['non-zero status/return codes in ''DALinloop_NI_wavfiles'' ']);
end

PAset(PAattns-PAattns+120);
SBset([7 7],[0 0]);
RPhalt(RP);
RPclear(RP);
SBset([],[]);
if (rc ~= 1)
   nelerror('''STM'': Error(s) detected within stimulus presentation loop');
end



