	function ratef(str)
% RATEF('a4') - computes and displays  rate versus level for driven and
% spontaneous rate from the file 'a4.DAT', which should be in standard
% lab format.
% RATEF('a4F') or 'a4F' filters data once (filterwidth = 3)
% Uses mex files (fast), but does not process file types other than rate-
% levels correctly.

plcolor=['k',  'r',  'g',  'm',  'b'];

% thatsall is a loop-breakout flag
% nplots counts plots superimposed, controls plot appearance
thatsall = 0;
nplots = 1;
hold off

while thatsall==0
	
% Parse for 'F' or 'f' at end of str
	if str(size(str,2))=='F' | str(size(str,2))=='f'
		str = str(1:size(str,2)-1);
		autofilt = 1;
	else
		autofilt = 0;
	end
   fname = sprintf('%s.DAT',str);
%
%fname1=fname;
%parm = get_parm_block(fname1);                % get parameter block
%dur = parm(47);
%if parm(39)==11
%   dur = parm(51);
%else if parm(41) ~= 11
%   dur = max([parm(47),parm(51)]);
%end
%tim1 = 10;
%tim2 = dur+10;
%tim4 = parm(94);
%tim3 = (tim2+tim4)/2;
   tim1 = 10;
   tim2 = 210;
   tim3 = 600;
   tim4 = 1000;
%
% read_rate() returns row vectors
%
   fname1=fname;
   [rate,atten,stim]=read_rate(fname1, tim1, tim2);
   fname1=fname;
   [spont,satten,stim]=read_rate(fname1, tim3, tim4);
   
   legend = sprintf('File %s, %s.', fname, stim);
%
% If autofiltering,
	if autofilt==1
		rate = trifilt(rate, 3);
		spont = trifilt(spont, 3);
		legend = strcat(legend,', F3');
	end
%
   displayloop = 0;
   while displayloop==0
      plot(atten, rate, strcat(plcolor(nplots),'-'), satten, spont, strcat(plcolor(nplots),'--'))
	  if nplots==1
         xlabel('Attenuation, dB')
         ylabel('Rate, spikes/s')
		 axlim = axis;
		 xleg = axlim(1)+5;
		 yleg = 0.95*axlim(4);
		 ydlt = 0.1*axlim(4);
	  end
	  text(xleg, yleg, legend)
	  yleg = yleg - ydlt;

      what = input('\nFn-filter, Sa5-superimpose, Na5-next: ','s');
      if what(1)=='F' | what(1)=='f'
	     nfw = sscanf(what(2:size(what,2)), '%g');
	     nfw = 2*floor(nfw/2) + 1;
	     legend = strcat(legend, ', F', sprintf('%g',nfw));
		 rate = trifilt(rate, nfw);
		 spont = trifilt(spont, nfw);
		 yleg = yleg + ydlt;
	  elseif what(1)=='N' | what(1)=='n'
		 str = what(2:size(what,2));
		 nplots = 1;
		 hold off
		 displayloop = 1;
	  elseif what(1)=='S' | what(1)=='s'
		 str = what(2:size(what,2));
		 nplots = nplots + 1;
		 if nplots>5
			 nplots = 2;
	     end
		 displayloop = 1;
		 hold on
	  else
		 displayloop = 1;
	     thatsall = 1;
	  end
   end
end
