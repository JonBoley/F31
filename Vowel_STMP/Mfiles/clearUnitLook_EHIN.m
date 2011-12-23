function clearUnitLook_EHIN()
% clears out all variables in a UnitLookEHIN .mat file
% *** but keeps characteristic delays ***

[FileName,PathName,FilterIndex] = uigetfile('UnitLook_EHIN*.mat');
FileName = [PathName FileName];

load(FileName);
clearvars -except NSCC_CDs_usec NSCC_peaks FileName
save(FileName,'NSCC_CDs_usec','NSCC_peaks');

disp([FileName ' cleared!']);
