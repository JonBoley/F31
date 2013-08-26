#include <stdlib.h>
#include <math.h>

#include "mex.h"

/*  Generates spike times from an inhomogeneous Poisson process with absolute and relative refractoriness.
	
	This function, SGfast, is much faster than the SGmodel function for two reasons:
		1) SGmodel uses a Bernoulli approximation to the Poisson process in each time bin, while SGfast
			uses the time transformation method.  This has the added benefit that it requires significantly
			fewer pseudo-random numbers (one per spike versus one per time bin for SGmodel).
		2) SGmodel calculates the relative refractory ratio from scratch at each time bin.  SGfast uses
			running approximations to the differential equations of which the exponentials in the relative
			refractory equation are solutions.
			
	Also, SGfast uses the MATLAB 'rand' function, whereas SGmodel uses the C 'rand' function compiled along
	with it.  However, the default 'rand' functions that come with many C compilers are often not very good
	ones, and MATLAB's is known to be a good one (and won't differ from compiler to compiler).  While SGmodel
	seeded its 'rand' function, SGfast does not seed the MATLAB 'rand' function.  This is beneficial in that
	the user now has more control over the pseudo-random number sequence used by the function, but could be
	a drawback if for some reason the results are dependent on the 'rand' function being randomly seeded and
	the user does not do so.
	
	Negative "rates" are acceptable, but this algorithm treats them as zero, i.e. it (half-wave) rectifies 
	the rate waveform.
	
	Calling syntax:
		[spktimes, {nspikes}] = SGfast([dt, nrep], rate, {deadtime, refracparams})
	
	where		spktimes	 is a vector containing the times at which spikes occurred (sec).
				nspikes		 is the number of spikes that occurred.
				dt			 is the sample period of the rate function (sec).  This is also used as
								the width of discrete time bins in the algorithm, and spike times will
								be multiples of 'dt'.
				nrep		 is the number of repetitions of the rate function.
				rate		 is the rate function vector (spikes/sec).
				deadtime	 is the dead time or absolute refractory period (sec).  DEFAULT: 0.00075
				refracparams is a vector [c0 s0 c1 s1] containing the parameters for the relative
									refractory period (dimensionless, sec, dimensionless, sec).  
									DEFAULT: [0.5 0.001 0.5 0.0125]
	
	
				
File: SGfast.c
Content: MEX-file for SGfast 
Date: 07/28/2003
Revision: 1.0

Copyright © 2003 B. Scott Jackson

*/


/*------------------------------------------------------------------------------------------------------
	The MEX-function for SGfast.
------------------------------------------------------------------------------------------------------*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
double		*rate, *vectorArg, T, dt, deadtime =  0.75e-3;
double		c0 = 0.5, c1 = 0.5, s0 = 0.001, s1 = 0.0125;
long		j, k, numRate, nrep;
long		NoutMax, Nout, randBufLen, randBufIndex;
mxArray		*randInputArray[1], *randOutputArray[1];
double		*spktimes, nspikes, meanRate, *tempOutputPtr, *randNums, *randDims;
double		deadtimeRnd, endOfLastDeadtime, refracMult0, refracMult1, refracValue0, refracValue1;
double		Xsum, unitRateIntrvl, time;
long		deadtimeIndex;

/* Check number of input arguments */
if ( (nrhs < 2) || (nrhs > 3) )  
	mexErrMsgTxt("Requires two or three input arguments.");

/* Check and get each input argument */
	/* First input argument */
if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || (mxGetNumberOfElements(prhs[0]) != 2) )
	mexErrMsgTxt("The first argument must be a two-element, real, double-precision vector.");
vectorArg = mxGetPr(prhs[0]);
dt = vectorArg[0];
if ( dt <= 0 )
	mexErrMsgTxt("The rate function sampling period must be greater than zero.");

	/* Second input argument */
nrep = floor(vectorArg[1]+0.5);
if ( (nrep < 0) || (nrep != vectorArg[1]) )
	mexErrMsgTxt("The number of repetitions must be a postive integer.");

if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) 
								|| ( (mxGetM(prhs[1]) != 1) && (mxGetN(prhs[1]) != 1) ) )
	mexErrMsgTxt("The second argument must be a real, double-precision vector.");
rate = mxGetPr(prhs[1]);
numRate = mxGetNumberOfElements(prhs[1]);

	/* Third input argument */
if (nrhs > 2) 
	if ( !mxIsEmpty(prhs[2]) )
		{
		if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || (mxGetNumberOfElements(prhs[2]) != 1) )
			mexErrMsgTxt("The third argument must be a real, double-precision scalar.");
		deadtime = mxGetScalar(prhs[2]);
		if (deadtime < 0)
			mexErrMsgTxt("The dead time must be non-negative.");
		}

	/* Fourth input argument */
if (nrhs > 3) 
	if ( !mxIsEmpty(prhs[3]) )
		{
		if ( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || (mxGetNumberOfElements(prhs[3]) != 4) 
											|| ( (mxGetM(prhs[3]) != 1) && (mxGetN(prhs[3]) != 1) ) )
			mexErrMsgTxt("The fourth argument must be a real, double-precision, four-element vector.");
		vectorArg = mxGetPr(prhs[3]);
		c0 = vectorArg[0];
		s0 = vectorArg[1];
		c1 = vectorArg[2];
		s1 = vectorArg[3];
		
		if ( (c0 < 0) || (c1 < 0) || (s0 <= 0) || (s1 <= 0) )
			mexErrMsgTxt("c0 & c1 must be non-negative and s0 & s1 must be positive.");
		if ( c0 + c1 != 1.0 )
			mexWarnMsgTxt("c0 + c1 is not equal to one.");
		}

/* Determine the mean of the rate vector */
meanRate = 0.0;
for (k=0; k<numRate; ++k)
	if (rate[k] > 0.0)  meanRate += rate[k];
meanRate /= numRate;

T = numRate * dt;  /* Total duration of rate function */

/* Create output buffer for spike times */
NoutMax = (long) ceil(meanRate * T * nrep);
if (NoutMax < 100)  NoutMax = 100;
spktimes = (double *)mxCalloc(NoutMax, sizeof(double));
Nout = 0;

/* Get a vector of pseudo-random numbers from MATLAB.  Calling MATLAB from a MEX-function is slow, 
	 so calling the MATLAB 'rand' function for each individual pseudo-random number needed greatly 
	 increases the computation time.																*/
if ( NoutMax < 1 )  mexErrMsgTxt("Random number buffer length must be positive.");
randBufLen = NoutMax+1;  /* Need 1 "extra" pseudo-random number to initialize the process */
randInputArray[0] = mxCreateDoubleMatrix(1, 2, mxREAL);
randDims = mxGetPr(randInputArray[0]);
randDims[0] = 1;
randDims[1] = randBufLen;
mexCallMATLAB(1, randOutputArray, 1, randInputArray, "rand");
randNums = mxGetPr(randOutputArray[0]);
randBufIndex = 0;

/* Calculate useful constants */
deadtimeIndex = (long) floor(deadtime/dt);	/* Integer number of discrete time bins within deadtime */
deadtimeRnd = deadtimeIndex*dt;				/* Deadtime rounded down to length of an integer number
													   						  of discrete time bins */

refracMult0 = 1 - dt/s0;  /* If y0(t) = c0*exp(-t/s0), then y0(t+dt) = y0(t)*refracMult0 */
refracMult1 = 1 - dt/s1;  /* If y1(t) = c1*exp(-t/s1), then y1(t+dt) = y1(t)*refracMult1 */

/* Calculate effects of a random spike before t=0 on refractoriness and the time-warping sum at t=0 */
endOfLastDeadtime = log( randNums[randBufIndex++] ) / rate[0] + deadtime;  /* End of last deadtime before t=0 */
refracValue0 = c0*exp(endOfLastDeadtime/s0);  /* Value of first exponential in refractory function */
refracValue1 = c1*exp(endOfLastDeadtime/s1);  /* Value of second exponential in refractory function */
Xsum = rate[0] * ( -endOfLastDeadtime + c0*s0*(exp(endOfLastDeadtime/s0)-1) 
											+ c1*s1*(exp(endOfLastDeadtime/s1)-1) );  /* Value of time-warping sum */
					/*  ^^^^ This is the "integral" of the refractory function ^^^^ (normalized by 'dt') */

/* Calculate first interspike interval in a homogeneous, unit-rate Poisson process (normalized by 'dt') */
unitRateIntrvl = -log( randNums[randBufIndex++] ) / dt;

/* NOTE: Both 'unitRateInterval' and 'Xsum' are divided (or normalized) by 'dt' in order to reduce calculation time.  
		This way we only need to divide by 'dt' once per spike (when calculating 'unitRateInterval'), instead of 
		multiplying by 'dt' once per time bin (when calculating the new value of 'Xsum').                            */

time = dt;
k=0;
for (j=0; j<nrep; ++j)  /* Loop through "stimulus" repetitions */
	{
	for (; (k<numRate) && (time<T); ++k, time+=dt, 
									refracValue0*=refracMult0, refracValue1*=refracMult1)  /* Loop through rate vector */
		{
		if (rate[k] > 0.0)  /* Nothing to do for non-positive rates, i.e. Xsum += 0 for non-positive rates. */
			{
			Xsum += rate[k] * (1 - refracValue0 - refracValue1);  /* Add rate*(refractory value) to time-warping sum */
			
			if ( Xsum >= unitRateIntrvl )  /* Spike occurs when time-warping sum exceeds interspike "time" in unit-rate process */
				{
				spktimes[Nout] = time;
				if (++Nout >= NoutMax)  /* If the output buffer for the spike times is too small . . . */
					{
					NoutMax += (long)ceil(meanRate * ((T-time) + (nrep-j-1)*T));
					if (NoutMax - Nout < 100)  NoutMax += 100;
					spktimes = (double *)mxRealloc(spktimes, NoutMax*sizeof(double));  /* . . . allocate additional memory . . . */
					if (spktimes == NULL)  mexErrMsgTxt("Out of Memory");
					
					/* . . . and get more pseudo-random numbers since we will have used all of them up. */
					randBufLen = NoutMax;
					randDims[1] = randBufLen;
					mxDestroyArray(randOutputArray[0]);
					mexCallMATLAB(1, randOutputArray, 1, randInputArray, "rand");
					randNums = mxGetPr(randOutputArray[0]);
					randBufIndex = 0;
					}
				
				unitRateIntrvl = -log( randNums[randBufIndex++] ) / dt;  /* Next interspike "time" in unit-rate process */
				Xsum = 0.0;
				
				/* Increase index and time to the last time bin in the deadtime, and reset (relative) refractory function */
				k += deadtimeIndex;
				time += deadtimeRnd;
				refracValue0 = c0;
				refracValue1 = c1;
				}
			}
		} /* End of rate vector loop */
		
		/* Reset index and time to begin a new repetion.  Don't just set to zero, since new repetition may start within the 
				deadtime following the last spike in the previous repetition.                                               */
		time -= T;
		k -= numRate;
	} /* End of "stimulus" repetition loop */

/* Delete spike(s) that occur after the last repetition of the rate function ends */
for (; (Nout>0)&&(spktimes[Nout-1]>T); --Nout)  spktimes[Nout-1] = mxGetNaN();

/* Finish up */
mxDestroyArray(randInputArray[0]);
mxDestroyArray(randOutputArray[0]);

nspikes = Nout;  /* Number of spikes that occurred. */
if (Nout > 0)
	{
	/* Make the output (spike) buffer size equal to the number of spikes. */
	tempOutputPtr = (double *)mxRealloc(spktimes, Nout*sizeof(double));
	if (tempOutputPtr == NULL)
		{
		mexWarnMsgTxt("Unable to deallocate the unused portion of memory in the output array.");
		for (k=Nout; k<NoutMax; ++k)
			spktimes[k] = mxGetNaN();  /* Fill unused part of the output (spike) buffer with NaN's. */
		Nout = NoutMax;
		}
	else
		spktimes = tempOutputPtr;  /* Able to reduce the buffer size. */

	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	mxFree(mxGetPr(plhs[0]));
	mxSetN(plhs[0], Nout);
	mxSetPr(plhs[0], spktimes);
	
	/* Output the number of spikes. */
	if (nlhs > 1)
		{
		plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(plhs[1]) = nspikes;
		}
	}
else  /* Return empty matrix since no spikes occurred. */
	{
	mxFree(spktimes);
	plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);

	/* Output the number of spikes, which in this case is zero. */
	if (nlhs > 1)
		{
		plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
		*mxGetPr(plhs[1]) = 0;
		}
	}
}
