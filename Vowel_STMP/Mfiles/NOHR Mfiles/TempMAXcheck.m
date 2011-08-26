COUNTmh=0;
TEMPmaxs=zeros(1,4*3*11)+NaN;
for iROW=1:4
   for iCOL=1:3
      for jj=1:11
         COUNTmh=COUNTmh+1;
         TEMP=max(PERhists{iROW,iCOL}{jj});
         if ~isnan(TEMP)
            TEMPmaxs(COUNTmh)=TEMP;
         end
      end
   end
end
