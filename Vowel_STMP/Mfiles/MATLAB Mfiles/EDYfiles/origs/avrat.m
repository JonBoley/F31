% AVRAT - compute rate at fixed abval across a set of curves.

loval = input('Minimum dB value for window: ');
hival = input('Maximum dB value for window: ');

arate = zeros(nralvs,1);

for jc=1:nralvs
	ind = abval(:,jc)<0 & abval(:,jc)>=loval & abval(:,jc)<=hival;
	arate(jc) = mean(drat(ind,jc));
end

clear jc
clear ind
