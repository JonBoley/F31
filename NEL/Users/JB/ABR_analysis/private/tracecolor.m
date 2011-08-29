function tracecolor

global han_abrcomp colorr

axes(han_abrcomp.text) 
xxx=0.5; yyy=0.5;
while xxx >=0.01 & xxx <= 0.99 & yyy >=0.01 & yyy <=0.99
	[xxx,yyy,marker] = ginput(1);
	x_contrib=floor(xxx*3)*5;
	y_contrib=ceil((1-yyy)*5);
	ind=x_contrib+y_contrib;
	if marker == 1 & xxx >=0.01 & xxx <= 0.99 & yyy >=0.01 & yyy <=0.99
		if strcmp(colorr(ind,:),'k')==1
			colorr(ind,:)='r';

		
		elseif strcmp(colorr(ind,:),'r')==1
			colorr(ind,:)='b';
		elseif strcmp(colorr(ind,:),'b')==1
			colorr(ind,:)='m';
		elseif strcmp(colorr(ind,:),'m')==1
			colorr(ind,:)='g';
		elseif strcmp(colorr(ind,:),'g')==1
			colorr(ind,:)='k';

		end
		comp_plot
	elseif marker == 3 & xxx >=0.01 & xxx <= 0.99 & yyy >=0.01 & yyy <=0.99
		colorr=['k';'r';'b';'m';'g';'k';'r';'b';'m';'g';'k';'r';'b';'m';'g'];
		comp_plot
	end
	

	
	
	
end