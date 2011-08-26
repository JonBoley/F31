function rm(picNum, figNum)

params = picview;
params.picArray = {picNum};
params.plotType ='rm';
params.checkPicType = 0;

if(exist('figNum', 'var'))
   params.figNum = figNum;
end

params = picview(params);
