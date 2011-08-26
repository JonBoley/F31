% RATV - makes rate level plots from lab-standard data, including PRE_GNYSYS
% data. This is a calling interface to ratepV(str).

global INDEXXYZ

str = input('Enter filenames (''a4-7,15[s][f]'', Superimpose, Filter: ','s');
[drat,srat,abval,npts,nralvs] = ratepV(str);

clear str
clear INDEXXYZ
