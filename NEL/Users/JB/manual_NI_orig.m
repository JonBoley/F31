function manual_NI(STIMind)
% File: manual_NI.m
%
% Written by MH 24oct2003.
%
% to debug NI board by hand

global signals_dir
global root_dir % Used to find correct control 'rco' file.

if ~exist('STIMind','var')
   STIMind=1;
end

load RP_eg.mat

disp(sprintf('\n\n***************'))

uRate=100000;
atten=40;
Gating.Duration=300;
Gating.Period=1000;
MAXtrigs=30;

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
   test_tone_list{ifiles} = sprintf('%sMH\\test_tone%d.wav', signals_dir, ifiles);
%   test_tone_list{ifiles} = sprintf('%sGE\\RSS_stimA_set01\\RSS_%03d.wav', signals_dir, ifiles);
end


[data,sr,rc] = nel_wavread(test_tone_list{STIMind});
% [data2,sr2,rc2] = nel_wavread(test_tone_list{10});
% data=[data' data2']';
% data(30001:30011)=ones(11,1);


if (rc == 0)
   nelerror(['Can''t read wavfile ''' static_bi.list{min(ind,end)} '''']);
else
   Chan0data = data;
end
% Rescape to 5 volts
Chan0data=Chan0data*5.0;

[statCodes retVals ireturnWAV] = d2a(1, uRate, Chan0data);
lineNum=56;
rc_d2a = 0;
if ( (~isempty(find(statCodes~=0))) | (~isempty(find(retVals~=0))) )
   rc_d2a = -1;
   disp(sprintf('\tManual_NI Location: %d',lineNum))
   disp([sprintf('\t\t') 'StatCodes: ' sprintf('%d ',statCodes')])
   disp([sprintf('\t\t') 'Return Vals: ' sprintf('%d ',retVals')])
   nelwarn(['non-zero status/return codes in ''DALinloop_NI_wavfiles'' ']);
end

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

%temp=input('Press Enter to start');
TRIGstart;


%temp=input('Press Enter to stop');
%pause(3)

for i=1:(MAXtrigs-20)
   pause(Gating.Duration/1000+.1);
%   pause(.8)
   disp(sprintf('Index=%d',i+1))   
   [statCodes retVals] = d2a(2);  % Terminate/Clean-up NI 6052 board (see file "d2a.c")
   lineNum=94;
   rc_d2a = 0;
   if ( (~isempty(find(statCodes~=0))) | (~isempty(find(retVals~=0))) )
      rc_d2a = -1;
      disp(sprintf('\tManual_NI Location: %d',lineNum))
      disp([sprintf('\t\t') 'StatCodes: ' sprintf('%d ',statCodes')])
      disp([sprintf('\t\t') 'Return Vals: ' sprintf('%d ',retVals')])
      nelwarn(['non-zero status/return codes in ''DALinloop_NI_wavfiles'' ']);
   end
   
   [statCodes retVals ireturnWAV] = d2a(1, uRate, Chan0data);
   lineNum=39;
   rc_d2a = 0;
   if ( (~isempty(find(statCodes~=0))) | (~isempty(find(retVals~=0))) )
      rc_d2a = -1;
      disp(sprintf('\tManual_NI Location: %d',lineNum))
      disp([sprintf('\t\t') 'StatCodes: ' sprintf('%d ',statCodes')])
      disp([sprintf('\t\t') 'Return Vals: ' sprintf('%d ',retVals')])
      nelwarn(['non-zero status/return codes in ''DALinloop_NI_wavfiles'' ']);
   end 
end


[statCodes retVals] = d2a(2);  % Terminate/Clean-up NI 6052 board (see file "d2a.c")
lineNum=94;
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



