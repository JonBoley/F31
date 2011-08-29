function change_weights

global data spl w num

x=min(spl)-10; distance=NaN*ones(1,num);
while x >= min(spl)-10 & x <= max(spl)+10
    [x,y,marker] = ginput(1);
	for i=1:num
		distance(1,i)=abs(x-spl(1,i))/(max(spl)+10-min(spl)-10)+abs(y-data.z.score(1,i))/max(data.z.score)/1.1;
	end
	[toss,index]=min(distance);
	if marker == 1 & x >= min(spl)-10 & x <= max(spl)+10
		w(1,index)=1;
	elseif marker ==3 & x >= min(spl)-10 & x <= max(spl)+10
		w(1,index)=0;
	end;

	X=ones(num,2);
	X(:,2)=spl';
	b=lscov3(X,data.z.score',w');
	data.z.slope=b(2,1);data.z.intercept=b(1,1);
	data.threshold=(3-b(1,1))/b(2,1);

	plot_data
	
end;