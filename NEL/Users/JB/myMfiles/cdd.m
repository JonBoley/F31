% CDD - Change to the data directory
global NelData SaveDirXYZ

SaveDirXYZ = cd;
cd(NelData.File_Manager.dirname);
% fprintf('\nRemember to call RDD to reset directory.\n')
