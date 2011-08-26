function doSCCmex(dist,debug)


if (exist('dist','var') ~= 1)
   dist = 0;
end
if (exist('debug','var') ~=1)
   debug = 0;
end

mex_files = {'InnerSCCmex.c' ...
   };
doNOHRmex(mex_files,dist,debug);
