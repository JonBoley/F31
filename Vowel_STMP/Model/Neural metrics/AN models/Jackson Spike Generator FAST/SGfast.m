%SGFAST Generates spike times from an inhomogeneous Poisson process with refractoriness.
%
%	SGFAST is a faster (and better) replacement for SGMODEL.  Its parameters are entirely
%	backwards compatible with SGMODEL.
%	
%	SGFAST([DT, NREP], RATE) returns spike times in seconds for an inhomogeneous Poisson 
%	process with a driving rate specified by the vector of RATE at sampling period of DT 
%	seconds and with refractoriness.  Any negative elements in RATE are treated as zeros.
%	DT is also used as the the width of discrete time bins in the algorithm, and spike times 
%	will be multiples of DT.  The input RATE is repeated NREP contiguous times, but the spike 
%	times are given relative to the beginning of the repetition (not from the beginning of 
%	the first repetition).  The absolute refractory period (i.e. deadtime) is 0.00075 sec 
%	long and the relative refractory function is:
%		(1 - c0*exp(-(t-td)/s0) - c1*exp(-(t-td)/s1), 
%	where 't' is the time since the last spike, 'td' is the length of the absolute refractory 
%	period, c0 = 0.5, s0 = 0.001 sec, c1 = 0.5, and s1 = 0.0125 sec.
%
%	SGFAST([DT, NREP], RATE, DEADTIME) uses an absolute refractory period of length DEADTIME
%	seconds.
%
%	SGFAST([DT, NREP], RATE, DEADTIME, [C0, S0, C1, S1]) also uses the relative refractory 
%	parameters C0, S0, C1, and S1.  S0 and S1 should be in seconds.  DEADTIME may be an empty 
%	matrix, in which case the default value of 0.00075 sec is used.  
%
%	[SPKTIMES, NSPIKES] = SGFAST(. . .) also returns the number of spikes in NSPIKES.  This is
%	not necessary since it will always be equal to the length of SPKTIMES.  However, this makes
%	SGFAST backwards compatible with SGMODEL, which outputs the number of spikes in its second 
%	output parameter.
%
%	NOTE: SGFAST uses the MATLAB pseudo-random number generator RAND without seeding it.  Thus, 
%	the calling function can and should handle seeding of the random number generator, if 
%	necessary.
%
%   See also SGMODEL.

%   Copyright © 2003 by B. Scott Jackson
%   Revision: 1.0    Date: 7/28/03
%
%	SGFAST is a MEX-function.  Its code is contained in the file 'SGfast.c'.
%
%	This function, SGfast, is much faster than the SGmodel function for two reasons:
%		1) SGmodel uses a Bernoulli approximation to the Poisson process in each time bin, while SGfast
%			uses the time transformation method.  This has the added benefit that it requires significantly
%			fewer pseudo-random numbers (one per spike versus one per time bin for SGmodel).
%		2) SGmodel calculates the relative refractory ratio from scratch at each time bin.  SGfast uses
%			running approximations to the differential equations of which the exponentials in the relative
%			refractory equation are solutions.
%			
%	Also, SGfast uses the MATLAB 'rand' function, whereas SGmodel uses the C 'rand' function compiled along
%	with it.  However, the default 'rand' functions that come with many C compilers are often not very good
%	ones, and MATLAB's is known to be a good one (and won't differ from compiler to compiler).  While SGmodel
%	seeded its 'rand' function, SGfast does not seed the MATLAB 'rand' function.  This is beneficial in that
%	the user now has more control over the pseudo-random number sequence used by the function, but could be
%	a drawback if for some reason the results are dependent on the 'rand' function being randomly seeded and
%	the user does not do so.
