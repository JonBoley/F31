#include "filters.hpp"
#include "complex.hpp"
/**
    The low-pass filter is in the form of cascade of first order lowpass filter\\
    Each of them is implemented as y(i) = c1*y(i-1)+c2*(x(i)+x(i-1))

    This function initializes the filter:\\
    1. initialize/reset the filter\\
    2. calculation c1 c2 thru bilinear transformation\\
    3. set the gain and order of the filter\\

    @author Xuedong Zhang
 */
void TLowPass::init(double _tdres,double _Fc,double _gain,int _LPorder)
{
  tdres = _tdres;
  double c = 2.0/tdres;
  Fc = _Fc;
  LPorder = _LPorder;
  c1LP = ( c - TWOPI*Fc ) / ( c + TWOPI*Fc );
  c2LP = TWOPI*Fc / (TWOPI*Fc + c);
  for(int i=0; i<=LPorder; i++) hc[i] = hcl[i] = 0.0;
  gain = _gain;
};
/**
  This function runs the low-pass filter
   @author Xuedong Zhang
 */
double TLowPass::run(double x)
{
  hc[0] = x*gain;
  for(int i=0; i<LPorder;i++)
    hc[i+1] = c1LP*hcl[i+1] + c2LP*(hc[i]+hcl[i]);
  for(int i=0; i<=LPorder;i++) hcl[i] = hc[i];
  return(hc[LPorder]);
};
/**

    The high-pass filter is in the form of a cascade of first order highpass filter\\
    Each of them is implemented as y(i) = c1*y(i-1)+c2*(x(i)-x(i-1))

    This function does several things:\\
    1. initialize/reset the filter\\
    2. calculation c1 c2 thru bilinear transformation\\
    3. set the gain and order of the filter\\

    @author Xuedong Zhang
*/
void THighPass::init(double _tdres,double _Fc,double _gain,int _HPorder)
{
  tdres = _tdres;
  double c = 2.0/tdres;
  Fc = _Fc;
  HPorder = _HPorder;
  c1HP = ( c - TWOPI*Fc ) / ( c + TWOPI*Fc );
  c2HP = c / (TWOPI*Fc + c);
  for(int i=0; i<=HPorder; i++) hc[i] = hcl[i] = 0.0;
  gain = _gain;
};
/**

   This function get the filtering output of the low-pass filter
   @author Xuedong Zhang
 */
double THighPass::run(double x)
{
  hc[0] = x*gain;
  for(int i=0; i<HPorder;i++)
    hc[i+1] = c1HP*hcl[i+1] + c2HP*(hc[i]-hcl[i]);
  for(int i=0; i<=HPorder;i++) hcl[i] = hc[i];
  return(hc[HPorder]);
};
/**
    The gammatone filter is in the form of cascade of first order lowpass filter
    shifted by the center frequency (from Carney,1993)\\
    Each of the low pass filter is implemented as y(i) = c1*y(i-1)+c2*(x(i)-x(i-1))

    This function does several things:\\
    1. initialize/reset the filter\\
    2. calculation c1 c2 thru bilinear transformation from tau\\
    3. set the gain and order of the filter\\

    @author Xuedong Zhang
*/
void TGammaTone::init(double _tdres,double _Fshift,double _tau,double _gain,int _order)
{
  tdres = _tdres;
  F_shift = _Fshift;
  delta_phase = -TWOPI*F_shift*tdres;
  phase = 0;
  tau = _tau;
  double c = 2.0/tdres; /* for bilinear transformation */
  c1LP = (tau*c-1)/(tau*c+1);
  c2LP = 1/(tau*c+1);
  gain = _gain;
  gtf_order = _order;
  for(int i=0; i<=gtf_order; i++) gtf[i].x = gtfl[i].x = gtf[i].y = gtfl[i].y = 0.0;
};
/**

   Reset the tau of the gammatone filter\\
   it recalculate the c1 c2 used by the filtering function
 */
void TGammaTone::settau(double _tau)
{
  tau = _tau;
  double dtmp = tau*2.0/tdres;
  c1LP = (dtmp-1)/(dtmp+1);
  c2LP = 1.0/(dtmp+1);
};
/**

   Pass the signal through the gammatone filter\\
   1. shift the signal by centeral frequency of the gamma tone filter\\
   2. low pass the shifted signal \\
   3. shift back the signal \\
   4. take the real part of the signal as output
   @author Xuedong Zhang
 */
double TGammaTone::run(double x)
{
  x = gain*x;
  phase += delta_phase;
  gtf[0] = compmult(x,compexp(phase));     // FREQUENCY SHIFT
  for(int j = 1; j <= gtf_order; j++)      // IIR Bilinear transformation LPF
    gtf[j] = comp2sum(compmult(c2LP,comp2sum(gtf[j-1],gtfl[j-1])),
		      compmult(c1LP,gtfl[j]));
  double out = REAL(compprod(compexp(-phase), gtf[gtf_order])); // FREQ SHIFT BACK UP
  for(int i=0; i<=gtf_order;i++) gtfl[i] = gtf[i];
  return(out);
};

