function [xval, yval] = cursor(fighan)

% CURSOR launches a crosshair cursor in the figure whose handle is
% fighan. Waits for user to position cursor and click the mouse. Returns
% the position of cursor in the plot in the figure, in units of the
% axis scaling of the figure.
%	 global INDEXXYZ        <---NOTE this must be in calling script
%    [xval,yval] = cursor(fighan);
% NOTE: if there is more than one plot in the figure, I don't know what
% this will return.

global INDEXXYZ

	figure(fighan)
	doclick=['INDEXXYZ=1;'];
	set(fighan, 'WindowButtonUpFcn', doclick)
	set(fighan, 'Pointer', 'crosshair')
	drawnow
	
% Get axis limits in units of pixels (same units as mouse posn.)
% axlimpix = [left bottom width, height].
	unitv = get(gca,'Units');
	set(gca,'Units', 'Pixels')
	axlimpix = get(gca,'Position');
	set(gca, 'Units', unitv)
	
% Wait for the mouse-click:
INDEXXYZ = 1;
INDEXXYZ = waitforbuttonpress;
   while INDEXXYZ==1
      drawnow
     % disp(INDEXXYZ); modified for PC BY M.A.

	end

% Get mouse position and convert to axis units. Recall axlim is
% [left right bottom top].
	mouselocw = get(fighan,'CurrentPoint');
	axlim = axis;
	xval = axlim(1) + (axlim(2)-axlim(1))*(mouselocw(1)-axlimpix(1))/axlimpix(3);
	yval = axlim(3) + (axlim(4)-axlim(3))*(mouselocw(2)-axlimpix(2))/axlimpix(4);

% Clean up the mess
	set(fighan, 'WindowButtonUpFcn', '')
	set(fighan, 'Pointer', 'arrow')

return
