% TESTP - test pointer functions
%

	fighan = gcf;
	figure(fighan);
	
	perd = 20*rand;
	plot([1:100],sin(2*pi*[1:100]/perd),'g-')
	
	doclick=['INDEX=1;'];
	set(fighan, 'WindowButtonUpFcn', doclick)
	set(fighan, 'Pointer', 'crosshair')
	x = 0;
	
	INDEX = 0;
	while INDEX==0
		x = x + 1;
		drawnow;
	end
	mouselocw = get(fighan,'CurrentPoint')
	mouselocs = get(0, 'PointerLocation')
	
	unitv = get(gca,'Units');
	set(gca,'Units', 'Pixels')
	axlimpix = get(gca,'Position')
	set(gca, 'Units', unitv)
	
	axlim = axis
	fprintf(1,'Mouse posn: x=%g, y=%g.\n', ...
	 axlim(1) + (axlim(2)-axlim(1))*(mouselocw(1)-axlimpix(1))/axlimpix(3), ...
	 axlim(3) + (axlim(4)-axlim(3))*(mouselocw(2)-axlimpix(2))/axlimpix(4));

	
	set(fighan, 'WindowButtonUpFcn', '')
	set(fighan, 'Pointer', 'arrow')
	
