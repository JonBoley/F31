function raster(picNum, figNum)

params = picview;
params.picArray = {picNum};
params.plotType ='raster';
params.checkPicType = 0;

if(exist('figNum', 'var'))
   params.figNum = figNum;
end

params.figNum = 101;

params = picview(params);