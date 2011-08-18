/*********************************************************************
* InnerSACmex.c 
* Created by: M. Heinz 03Sept2003
* 
* This is a MEX-file for MATLAB.
*
* Computes Shuffled Auto Correlogram from a Spike Matrix
* 
* Call:  [intsMEX,TOTALints] = InnerSACmex(SpikeMAT1,NUMspikes,TOTALspikes);
* 		SpikesMAT1: (Kmax x NUMreps) matrix with spike trains (NaNs fill end of each column for shorter spike trains)
*		NUMspikes: (1x NUMreps) vector with number of spikes in each train
*		TOTALspikes: needed to allocate intsMEX
*
*		intsMEX: vector with all-order intervals
*		TOTALints: number of actual ints (non NaN)
*		
*********************************************************************/

#include "mex.h"

//InnerShufAutoCorrMEX(SpikeMAT1,NUMspikes,NUMspikeREPS,Kmax,ints,TOTALints);
void InnerSACmex(double *SpikeMAT1, double *NUMspikes, long NUMspikeREPS, long Kmax, double *ints, double TOTALints[])
{
   long REPindREF,SpikeIND,REPindCOMP,COMPSpikeIND,intIND=0;
   
   for (REPindREF=0; REPindREF<NUMspikeREPS; REPindREF++) {
	   for (SpikeIND=0; SpikeIND<NUMspikes[REPindREF]; SpikeIND++) {
		   for (REPindCOMP=0; REPindCOMP<NUMspikeREPS; REPindCOMP++) {
			   if (REPindCOMP!=REPindREF) {
				   for (COMPSpikeIND=0; COMPSpikeIND<NUMspikes[REPindCOMP]; COMPSpikeIND++) {
					   ints[intIND]=SpikeMAT1[REPindCOMP*Kmax+COMPSpikeIND]-SpikeMAT1[REPindREF*Kmax+SpikeIND];
						intIND++;
      }  }  }	}  }
	TOTALints[0]=intIND;
}


// the gateway function
// [intsMEX,TOTALints] = InnerSACmex(SpikeMAT1,NUMspikes,TOTALspikes);
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
   double   *SpikeMAT1,*intsMEX,*NUMspikes; // input-output arrays
	long     TOTALspikes; // input scalar
   double   *TOTALints;  // output scalar
   long     NUMspikeREPS,Kmax;

	int      i,j;

//	__asm int 3h ; // Setting a breakpoint for debugging (compile with the -g flag!)
							// For "Nel_matlab", run 'doallmex' with '(0,1)' flags.
							// Run Nel software as usual, and Matlab will automatically
							//  break at this line and open this file in Visual Studio
							//  debugger environment.

   /*  check for proper number of arguments */
   if(nrhs!=3)
    mexErrMsgTxt("Three inputs required.");
   if(nlhs!=2)     mexErrMsgTxt("One output required.");
  
  	/*  create a pointer to the input matrix SpikeMAT1 */
   // SpikeMAT1 is stored in MATLAB as: each row is one spike train
	//   in MEX/C, linear indexing goes down columns first, so we pass the transpose to 
	//   allow easier indexing
   SpikeMAT1 = mxGetPr(prhs[0]);  // each column is one spike train
	/*  get the dimensions of the matrix input SpikeMAT1 */
   NUMspikeREPS = mxGetN(prhs[0]);
   Kmax = mxGetM(prhs[0]);
  
	/*  create a pointer to the input vector NUMspikes */
   NUMspikes = mxGetPr(prhs[1]);

   /*  get the scalar input TOTALspikes */
   TOTALspikes = mxGetScalar(prhs[2]);

   /* check to make sure the last input argument is a scalar */
   if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
      mxGetN(prhs[2])*mxGetM(prhs[2])!=1 ) {
      mexErrMsgTxt("Input TOTALspikes must be a scalar.");
   }

   
   /*  set the output pointer to the output matrix intsMEX */
	plhs[0] = mxCreateDoubleMatrix(1,TOTALspikes*(NUMspikeREPS-1)*Kmax, mxREAL);
   /*  create a C pointer to a copy of the output matrix intsMEX */
   intsMEX = mxGetPr(plhs[0]);
  
   /*  set the output pointer to the output scalar TOTALints */
   plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
   /*  create a C pointer to a copy of the output variable TOALints */
   TOTALints = mxGetPr(plhs[1]);
   

   /*  call the C subroutine */
   InnerSACmex(SpikeMAT1,NUMspikes,NUMspikeREPS,Kmax,intsMEX,TOTALints);
	
	}
