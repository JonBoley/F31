%%% readme.txt for zbcatmodel %%%

This is Version 2 of the public distribution of the code for the cat auditory
periphery model of:

        Zilany, M. S. A. and Bruce, I. C. (2006). "Modeling auditory-nerve
        responses for high sound pressure levels in the normal and impaired
        auditory periphery," Journal of the Acoustical Society of
        America 120(3):1446�1466,

   and
 
        Zilany, M. S. A. and Bruce, I. C. (2007). "Representation of the vowel
        /eh/ in normal and impaired auditory nerve fibers: Model predictions of
        responses in cats," Journal of the Acoustical Society of America
        122(1):402�417.

Please cite these papers if you publish any research results obtained with this
code or any modified versions of this code.

*** Change History ***

Version 2.0:-

-  Model now runs at sampling rates of 100, 200 and 500 kHz:
   The major improvement in this version of the code is that it is able run at
   100 kHz for model fibers with characteristic frequencies (CFs) up to 20 kHz.
   To model fibers up to CFs of 40 kHz, a minimum sampling rate of 200 kHz
   should be used.  The flexibility in sampling rate is made possible by the
   changes to the middle-ear filter and spike generator code described below.

-  New middle-ear filter code:
   A new formulation of the middle-ear filter (by Rasha Ibrahim) allows for an
   accurate and stable middle-ear transfer function for sampling rates between
   100 and 500 kHz.

-  New spike generator code:
   We have incorporated the spike generator code written by B. Scott Jackson (bsj22@cornell.edu) 
   (Scott's original code is available from Laurel Carney's web site at:
    http://www.urmc.rochester.edu/smd/Nanat/faculty-research/lab-pages/LaurelCarney/auditory-models.cfm)
   This code uses a renewal process (i.e., inter-spike interval) approach that
   is more efficient than the old Bernoulli (i.e., spike per time bin) method.

-  Improved speed and reduced memory consumption:
   Running the model at a lower sampling rate produces a major speed increase
   and reduction in memory usage.  In addition, the new spike generator code is
   substantially faster than the old method.  We have also found some savings that can
   be made in memory allocation.

Version 1.1:-

-  Switched to new CA gain versus CF function:
   The only difference between the models presented in the two papers is
   the function for cochlear amplifier (CA) gain versus characteristic frequency (CF),
   given by Eq. (6) in Zilany and Bruce (2006) and by Eq. (1) in Zilany and Bruce (2007).
   The new CA gain versus CF function of Zilany and Bruce (2007) is now used by default.
   If you wish to use the old function, then you will need to uncomment line 288 of
   zbcatmodel.c and comment out line 289, and then recompile the code.

-  Added unpublished arbitrary spont rate feature:
   This version of the code has an "unpublished feature" of being able to model
   AN fibers with different spontaneous rates (from 0 to 150 spikes/s) and the
   corresponding changes in the threshold and the rate-level function. A spontaneous
   rate of 50 spikes/s (before refractory effects) was used in the Zilany and Bruce
   papers.

Version 1.0:-

-  Original Public Release


The Matlab and C code included with this distribution is designed to be
compiled as a Matlab MEX file, i.e., the compiled model MEX function will run
as if it were a Matlab function.  The code can be compiled within Matlab using
the function:

    mexzbcatmodel_lcc.m     if you are using the lcc compiler supplied with
                            Matlab on a Win32 system; or

    mexzbcatmodel_unix.m    if you are using a compiler such as gcc on a Unix
                            system.

The reason for the difference is that Unix compilers typically include the
drand48() random number generator, which has the precision required for spike
generation in this model.  The lcc compiler does not include drand48(), so I
have made use of the Gnu Scientific Library (GSL; http://www.gnu.org/software/gsl/)
equivalent.  The required files from the GSL are included with this distribution,
as is a copy of the GNU General Public License (gpl.txt), under which the GSL
is released.

Note that since version 1.1 of the code, Microsoft C compilers are no longer supported,
because there are issues with recent versions of the  Microsoft C compilers and
recent Matlab releases.  Please switch to using the lcc compiler that comes with
Matlab if you are using MS Windows.

Once you have compiled the MEX file in Matlab, type:

    help zbcatmodel

for instructions on how to call the MEX function.

I have also included:-

1. a sample Matlab script "testzbcatmodel.m" for setting up an acoustic stimulus
   and the model parameters and running the model, and

2. a function "fitaudiogram.m" for estimating the parameters for outer and
   inner hair cell impairment, Cohc and Cihc, respectively, for a given
   audiogram.


ACKNOWLEDGMENTS

Some of this code is based on code written by Xuedong (Frank) Zhang, Michael
Heinz, Ian Bruce and Laurel Carney for the model of:

    Zhang, X., Heinz, M. G., Bruce, I. C., and Carney, L. H. (2001). "A
    phenomenological model for the responses of auditory-nerve fibers: I.
    Nonlinear tuning with compression and suppression," J. Acoust. Soc. Am. 109,
    648�670,

and by Qing Tan and Laurel Carney for the model of:

    Tan, Q. and Carney, L. H. (2003). �A phenomenological model for the
    responses of auditory nerve fibers. II. Nonlinear tuning with a
    frequency glide,� J. Acoust. Soc. Am. 114, 2007�2020.

%%% � Ian C. Bruce (ibruce@ieee.org), M. S. Arefeen Zilany and Rasha Ibrahim, June 2006 - December 2007 %%%
