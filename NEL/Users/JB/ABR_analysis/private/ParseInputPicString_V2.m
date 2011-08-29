function [picnums] = ParseInputPicString_V2(picst)

% Takes the input number string (eg, '5-7,9') and turns it into an array
% of picture numbers, picnums=[5,6,7,9]

c='0';
i=0;j=1;numpics=1;dashflag=0;
while i<length(picst)
   while c~='-' & c~=',' & i+j~=length(picst)+1
      b(j)=picst(i+j);
      c=b(j);
      j=j+1;
   end
   if c=='-' | c==','
      b=b(1:end-1);
   end
   if dashflag==1
      try
         upto=str2num(b);
      catch
         error('Can''t parse picture numbers.');
      end
      numdash=upto-picnums(numpics-1);
      for k=1:numdash
         picnums(k+numpics-1)=picnums(numpics-1)+k;
      end
      numpics=length(picnums);
   else  % if dashflag==1
      try
         picnums(numpics)=str2num(b);
      catch
         error('Can''t parse picture numbers!\n');
      end
   end
   clear b;
   i=i+j-1;
   j=1;
   if c=='-'
      dashflag=1;
   else
      dashflag=0;
   end
   c='0';
   numpics=numpics+1;
end  % while i<length(picst)