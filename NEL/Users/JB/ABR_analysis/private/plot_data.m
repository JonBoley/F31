function plot_data

global 	data han spl attn line_width abr_Stimuli abr_time abr w...
	upper_y_bound lower_y_bound y_shift padvoltage num thresh_mag reff

%clear out contents of all axes
set([han.abr_panel han.amp_panel han.lat_panel han.z_panel han.peak_panel],'NextPlot','replacechildren')
axes(han.abr_panel); plot(0,0,'-w'); axes(han.amp_panel); plot(0,0,'-w');
axes(han.lat_panel); plot(0,0,'-w');
axes(han.z_panel); plot(0,0,'-w'); axes(han.peak_panel); plot(0,0,'-w');


%abr panel
axes(han.abr_panel); set(han.abr_panel,'NextPlot','Add');
plot([abr_Stimuli.start_template abr_Stimuli.start_template abr_Stimuli.end_template abr_Stimuli.end_template...
	abr_Stimuli.start_template],[upper_y_bound(1,1)-padvoltage lower_y_bound(1,1)+padvoltage lower_y_bound(1,1)+padvoltage...
	upper_y_bound(1,1)-padvoltage upper_y_bound(1,1)-padvoltage],'-r','LineWidth',line_width)
for i=1:num
	plot(abr_time,abr(:,i)+y_shift(1,i),'-k',[abr_Stimuli.start abr_Stimuli.end],[upper_y_bound(1,i) upper_y_bound(1,i)],'-k',...
		'LineWidth',line_width)
	text(abr_Stimuli.start+(abr_Stimuli.end-abr_Stimuli.start)*0.01,upper_y_bound(1,i)-0.5*padvoltage,num2str(spl(i),'%10.1f'),...
		'fontsize',10,'horizontalalignment','left','VerticalAlignment','middle','color','b')
	text(abr_Stimuli.start+(abr_Stimuli.end-abr_Stimuli.start)*0.15,upper_y_bound(1,i)-0.5*padvoltage,...
		['(' num2str(-attn(i),'%10.1f') ')'],'fontsize',10,'horizontalalignment','left','VerticalAlignment','middle','color','k')
	if get(han.cbh,'Value')==get(han.cbh,'Max') & ~isempty(reff)
		for ii=1:length(reff.abrs.waves(:,1))
			reffattn=reff.abrs.waves(ii,2)+reff.abrs.thresholds(1,4);
			if reffattn==attn(1,i)
				plot(abr_time,reff.abrs.waves(ii,3:end)+y_shift(1,i),'-m','LineWidth',line_width)
			end;
		end;
	end;
end;




%peak panel
axes(han.peak_panel);
plot(data.x(1,:),data.y_forfig(1,:),'r+',data.x(3,:),data.y_forfig(3,:),'b+',...
	data.x(5,:),data.y_forfig(5,:),'m+',data.x(7,:),data.y_forfig(7,:),'k+','LineWidth',line_width)
text(data.x(1,:),data.y_forfig(1,:),'1','FontSize',10,'Color','r','horizontalalignment','left','VerticalAlignment','bottom')
text(data.x(2,:),data.y_forfig(2,:),'----','FontSize',10,'Color','r','horizontalalignment','center','VerticalAlignment','middle')
text(data.x(3,:),data.y_forfig(3,:),'2','FontSize',10,'Color','b','horizontalalignment','left','VerticalAlignment','bottom')
text(data.x(4,:),data.y_forfig(4,:),'----','FontSize',10,'Color','b','horizontalalignment','center','VerticalAlignment','middle')
text(data.x(5,:),data.y_forfig(5,:),'3','FontSize',10,'Color','m','horizontalalignment','left','VerticalAlignment','bottom')
text(data.x(6,:),data.y_forfig(6,:),'----','FontSize',10,'Color','m','horizontalalignment','center','VerticalAlignment','middle')
text(data.x(7,:),data.y_forfig(7,:),'4','FontSize',10,'Color','k','horizontalalignment','left','VerticalAlignment','bottom')
text(data.x(8,:),data.y_forfig(8,:),'----','FontSize',10,'Color','k','horizontalalignment','center','VerticalAlignment','middle')

%latency panel
set(han.lat_panel,'Box','on','XLim',[min(spl)-10 max(spl)+10],'YLim',[0.8 2.8],'XGrid','on','YGrid','on');
if isempty(reff)
	axes(han.lat_panel); plot(...%%% -6.635 is correction for delay of TDT/ER2
		spl,data.x(1,:)-6.635,'-r*',spl,data.x(3,:)-6.635,'-b*',...
		spl,data.x(5,:)-6.635,'-m*',spl,data.x(7,:)-6.635,'-k*','LineWidth',line_width)
else
	axes(han.lat_panel); plot(...%%% -6.635 is correction for delay of TDT/ER2
		spl,data.x(1,:)-6.635,'-r*',spl,data.x(3,:)-6.635,'-b*',...
		spl,data.x(5,:)-6.635,'-m*',spl,data.x(7,:)-6.635,'-k*',...
		reff.abrs.x(:,2),reff.abrs.x(:,3)-6.635,':r.',...
		reff.abrs.x(:,2),reff.abrs.x(:,5)-6.635,':b.',...
		reff.abrs.x(:,2),reff.abrs.x(:,7)-6.635,':m.',...
		reff.abrs.x(:,2),reff.abrs.x(:,9)-6.635,':k.','LineWidth',line_width)
end;

%zscore panel
axes(han.z_panel)
if isempty(reff)
	set(han.z_panel,'Box','on','XLim',[min(spl)-10 max(spl)+10],'YLim',[0 max(data.z.score)*1.1],...
		'XTickLabel',[],'XGrid','on','YGrid','on','NextPlot','Add');
	text(min(spl)-10,max(data.z.score),[' \theta      ' num2str(data.threshold,'%10.1f') ' dB SPL'],...
		'FontSize',10,'Color','r','horizontalalignment','left','VerticalAlignment','bottom')
	for i=1:num
		if w(1,i) ~= 0
			plot(spl(1,i),data.z.score(1,i),'r*')
		else
			plot(spl(1,i),data.z.score(1,i),'k*')
		end;
	end;
	plot([0 100],[3 3],'-r',[0 100],[data.z.intercept data.z.intercept+100*data.z.slope],'-r','LineWidth',line_width)
else
	set(han.z_panel,'Box','on','XLim',[min(spl)-10 max(spl)+10],'YLim',[0 max([data.z.score reff.abrs.z.score(:,3)'])*1.1],...
		'XTickLabel',[],'XGrid','on','YGrid','on','NextPlot','Add');
	text(min(spl)-10,max([data.z.score reff.abrs.z.score(:,3)']),[' \theta      ' num2str(data.threshold,'%10.1f') ' dB SPL'],...
		'FontSize',10,'Color','r','horizontalalignment','left','VerticalAlignment','bottom')
	text(min(spl)-10,0.9*max([data.z.score reff.abrs.z.score(:,3)']),[' \theta ref  ' num2str(reff.abrs.thresholds(1,2),'%10.1f') ' dB SPL'],...
		'FontSize',10,'Color','r','horizontalalignment','left','VerticalAlignment','bottom')
	for i=1:num
		if w(1,i) ~= 0
			plot(spl(1,i),data.z.score(1,i),'r*')
		else
			plot(spl(1,i),data.z.score(1,i),'k*')
		end;
	end;
	for i=1:length(reff.abrs.z.score')
		if reff.abrs.z.score(i,4) ~= 0
			plot(reff.abrs.z.score(i,2),reff.abrs.z.score(i,3),'r.')
		else
			plot(reff.abrs.z.score(i,2),reff.abrs.z.score(i,3),'k.')
		end;
	end;
	plot([0 100],[3 3],'-r',[0 100],[data.z.intercept data.z.intercept+100*data.z.slope],'-r','LineWidth',line_width)
	plot([0 100],[reff.abrs.z.par(1,2) reff.abrs.z.par(1,2)+100*reff.abrs.z.par(1,3)],':r','LineWidth',line_width)
end;

%amp panel
axes(han.amp_panel)
if isempty(reff)
	set(han.amp_panel,'Box','on','XLim',[min(spl)-10 max(spl)+10],'YLim',[0 max(data.amp)*1.1],'XGrid','on','YGrid','on')
	axes(han.amp_panel); plot(spl,data.amp,'-k*',...
		spl,data.amp_null,'k.',[0 100],[thresh_mag thresh_mag],'k-','LineWidth',line_width)
	text(min(spl)-10,max(data.amp),[' \theta      ' num2str(data.amp_thresh,'%10.1f') ' dB SPL'],...
		'FontSize',10,'horizontalalignment','left','VerticalAlignment','bottom')
else
	set(han.amp_panel,'Box','on','XLim',[min(spl)-10 max(spl)+10],'YLim',[0 max([data.amp reff.abrs.amp(:,3)'])*1.1],...
		'XGrid','on','YGrid','on')
	axes(han.amp_panel); plot(spl,data.amp,'-k*',spl,data.amp_null,'k.',[0 100],[thresh_mag thresh_mag],'k-',...
		reff.abrs.amp(:,2),reff.abrs.amp(:,3),':k.','LineWidth',line_width)
	text(min(spl)-10,max([data.amp reff.abrs.amp(:,3)']),[' \theta      ' num2str(data.amp_thresh,'%10.1f') ' dB SPL'],...
		'FontSize',10,'horizontalalignment','left','VerticalAlignment','bottom')
	text(min(spl)-10,0.9*max([data.amp reff.abrs.amp(:,3)']),[' \theta ref  ' num2str(reff.abrs.thresholds(1,3),'%10.1f') ' dB SPL'],...
		'FontSize',10,'horizontalalignment','left','VerticalAlignment','bottom')
end;