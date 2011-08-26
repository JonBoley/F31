function syncs = compsync(data, nbn, npsts, binw, T);

% Each column of data(nbn,npsts) contains a PST histogram. Computes
% |synchronized rates| at frequencies 1/T, where T is a row vector containing
% the value of T for each column in data; returns syncs as a row vector.

	sines = zeros(npsts, nbn);
	coses = zeros(npsts, nbn);
	tm = zeros(1,nbn);
	for j=2:nbn
		tm(j) = tm(j-1) + binw;
	end
	for j=1:npsts
		sines(j,:) = sin(2*pi*tm/T(j));
		coses(j,:) = cos(2*pi*tm/T(j));
	end
	ipt = diag(sines*data(1:nbn, :))/nbn;
	rpt = diag(coses*data(1:nbn, :))/nbn;
	syncs = sqrt(ipt.^2 + rpt.^2)';
	
% OK, now that that's done, set to -1 the syncs above the Nyquist freq
	nyqT = 2*binw;
	syncs(T<nyqT) = -1;
return
