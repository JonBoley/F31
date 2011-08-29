function peak1af

global num abr attn y_shift data

for pk=1:4
	if ~isnan(data.x(pk*2-1,1))
		data.x(pk*2,1)=NaN; data.y(pk*2,1)=NaN; data.y_forfig(pk*2,1)=NaN; old_val=100;
		for bin=bin_of_time(data.x(pk*2-1,1)):bin_of_time(data.x(pk*2-1,1)+1)
			if abr(bin,1) < abr(bin-1,1) & abr(bin,1) < abr(bin+1,1) & abr(bin,1) < old_val
				data.x(pk*2,1)=time_of_bin(bin); data.y(pk*2,1)=abr(bin,1); data.y_forfig(pk*2,1)=abr(bin,1)+y_shift(1,1);
				old_val=abr(bin,1);
			end
		end
	end
end
for i=2:num
	add_attn=attn(1,i-1)-attn(1,i); add_attn2=attn(1,1)-attn(1,i); 
	for pk=1:4
		if data.z.score(1,i)>3 & ~isnan(data.x(pk*2-1,i-1))
			exp_lat= 0.5*((data.x(pk*2-1,i-1)+add_attn/40)+(data.x(pk*2-1,1)+add_attn2/40));
			data.x(pk*2,i)=NaN; data.y(pk*2,i)=NaN; data.y_forfig(pk*2,i)=NaN; old_val=-100;
			for bin=bin_of_time(exp_lat-0.5):bin_of_time(exp_lat+0.5)
				if abr(bin,i) > abr(bin-1,i) & abr(bin,i) > abr(bin+1,i) & abr(bin,i) > old_val
					data.x(pk*2-1,i)=time_of_bin(bin); data.y(pk*2-1,i)=abr(bin,i); data.y_forfig(pk*2-1,i)=abr(bin,i)+y_shift(1,i);
					old_val=abr(bin,i);
				end
			end
		end
		
		if ~isnan(data.x(pk*2-1,i))
			data.x(pk*2,i)=NaN; data.y(pk*2,i)=NaN; data.y_forfig(pk*2,i)=NaN; old_val=100;
			for bin=bin_of_time(data.x(pk*2-1,i)):bin_of_time(data.x(pk*2-1,i)+1)
				if abr(bin,i) < abr(bin-1,i) & abr(bin,i) < abr(bin+1,i) & abr(bin,i) < old_val
					data.x(pk*2,i)=time_of_bin(bin); data.y(pk*2,i)=abr(bin,i); data.y_forfig(pk*2,i)=abr(bin,i)+y_shift(1,i);
					old_val=abr(bin,i);
				end
			end
		end
	end
end

plot_data


