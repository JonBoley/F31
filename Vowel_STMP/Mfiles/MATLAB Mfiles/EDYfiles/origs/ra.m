% RATE - makes rate level plots from lab-standard data. This is a calling
% interface to ratep(str).

global INDEXXYZ

str = input('Enter filenames (''a4-7,15[s][f]'', Superimpose, Filter: ','s');
[drat,srat,abval,npts,nralvs] = ratep(str);

clear str
clear INDEXXYZ
