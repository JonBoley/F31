function change_spectrum

global z han_dp zanimal zdate


x=300;
while x >= 300 & x <= 12000
	[x,y,marker] = ginput(1);
	for i=1:length(z.DpoaeData(:,1))
		fdiff(i,1)=abs(x-z.DpoaeData(i,3));
	end
	[toss,index]=min(fdiff);
	
	
	if x >= 300 & x <= 12000
		
		if marker==3 & x >= 300 & x <= 12000
			if isnan(z.DpoaeData(index,4))
				[toss, index2]=min(abs(z.Dpoaefreqs(1,:)-z.DpoaeData(index,3))); %index of expected DP freq
				z.DpoaeData(index,4)=max(z.DpoaeSpectra(index,index2-20:index2+20));
			elseif ~isnan(z.DpoaeData(index,4))
				z.DpoaeData(index,4)=NaN;
			end
		end
		axes(han_dp.dp_curve); plot(z.DpoaeData(:,3),z.DpoaeData(:,4),'*k-',z.Dpoaefreqs(1,:),z.DpoaeSpectra(index,:),'r-',...
			[z.DpoaeData(index,3) z.DpoaeData(index,3)],[0 100],'y:')
%	z.DpoaeData(index,3),z.DpoaeData(index,4),'ro')
		%set(han_dp.dp_curve,'Box','on','XGrid','on','YGrid','on','Xscale','log','Xlim',[300 12000],...
			%'XTick',[300 400 500 600 800 1000 2000 3000 4000 5000 6000 8000 10000],'YLim',[0 100])
      set(han_dp.dp_curve,'Box','on','XGrid','on','YGrid','on','Xscale','log','Xlim',[300 10000],...
         'YLim',[0 100])
      axes(han_dp.dp_curve); title(['Chin' char(zanimal) ' on ' char(zdate)],'FontSize',14)
		ylabel('DP amplitude (dB SPL)','fontsize',14)
		xlabel('DP frequency (Hz)','fontsize',14)
		
	end
	
end;


