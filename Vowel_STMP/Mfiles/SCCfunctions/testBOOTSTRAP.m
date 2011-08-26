% TESTBOOTSTRAP
% - work out indices


X=6
inds=1:X;
percent=.8;
Nreps=length(inds);
Ninclude=ceil(percent*Nreps)
Nexclude=max([2 Nreps-Ninclude])
Navgs=floor(Nreps/Nexclude)

BOOTinds=cell(1,Navgs);

for BOOTrep=1:Navgs
	BOOTinds=setdiff(inds,(Nreps-Nexclude+1:Nreps)-(BOOTrep-1)*Nexclude)
end
