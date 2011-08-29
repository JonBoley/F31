function peak(type,id)

global abr_Stimuli num abr upper_y_bound lower_y_bound y_shift data


x=abr_Stimuli.start;
while x >= abr_Stimuli.start & x <= abr_Stimuli.end
    [x,y,marker] = ginput(1);
	for i=1:num
		if y>lower_y_bound(1,i) & y<upper_y_bound(1,i)
			index=i;
		end;
	end;
	if type == 'p' & x >= abr_Stimuli.start & x <= abr_Stimuli.end
		[yy, xx]=max(abr(bin_of_time(x-0.15):bin_of_time(x+0.15),index));
		xxx=time_of_bin(xx)+time_of_bin(bin_of_time(x-0.15));
		if marker==1
			data.x(id*2-1,index)=xxx; data.y(id*2-1,index)=yy; data.y_forfig(id*2-1,index)=yy+y_shift(1,index);
		elseif marker==3
			data.x(id*2-1,index)=NaN; data.y(id*2-1,index)=NaN; data.y_forfig(id*2-1,index)=NaN;
		end;
	elseif type == 'n' & x >= abr_Stimuli.start & x <= abr_Stimuli.end
		[yy, xx]=min(abr(bin_of_time(x-0.15):bin_of_time(x+0.15),index));
		xxx=time_of_bin(xx)+time_of_bin(bin_of_time(x-0.15));
		if marker==1
			data.x(id*2,index)=xxx; data.y(id*2,index)=yy; data.y_forfig(id*2,index)=yy+y_shift(1,index);
		elseif marker==3
			data.x(id*2,index)=NaN; data.y(id*2,index)=NaN; data.y_forfig(id*2,index)=NaN;
		end;

	end;

plot_data

end;




