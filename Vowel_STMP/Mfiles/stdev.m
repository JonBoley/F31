function y=stdev(x)
%STDEV standard deviation

y = sqrt(sum((x(~isnan(x))-mean(x(~isnan(x)))).^2/numel(x)))/sqrt(numel(x));
