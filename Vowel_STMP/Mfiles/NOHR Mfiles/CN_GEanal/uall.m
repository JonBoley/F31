function rc = uall(trackNum, unitNum, data_folder)

% Written by GE 19Mar2004.

if ~exist('data_folder')
   data_folder = [cd '\'];
end

figure(1); clf;
figure(2); clf;
figure(3); clf;
figure(5); clf;

urss(trackNum, unitNum, data_folder);
usmp(trackNum, unitNum, data_folder);
uchar(trackNum, unitNum, data_folder);
url(trackNum, unitNum, data_folder);

