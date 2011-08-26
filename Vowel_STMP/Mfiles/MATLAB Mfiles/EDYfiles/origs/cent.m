% CENT: Compute percentile scores for a data vector da.
%
if exist('da')==0
	fprintf('*** Got to give me a data vector. ***\n')
else
	mv = mean(da);
	sd = std(da);
	fprintf('Mean=%g,  stdev=%g.\n',mv,sd);
	zv = (da-mv)/sd;
	fprintf(' Data   z-score \n')
	for j=1:max(size(da))
		fprintf('%g  %g\n',da(j),zv(j))
	end
end
