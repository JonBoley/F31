#ifndef _FILTER_HPP
/**
   This file defines several filters that is used in the program
 */
#define _FILTER_HPP
#include <math.h>
#include <stdlib.h>
#include "complex.hpp"

#ifndef TWOPI
#define TWOPI 6.2831853
#endif
/// Maximum order of the filter
#define MAX_ORDER 10

// ----------------------------------------------------------------------------
/**
   Base Class of the filter
 */
class TFilter{
public:
  ///input one signal to the filter and get the output
  virtual double run(double x){return x;};
  ///get the order of the filter
  virtual int GetOrder(void){return 0;}
};

// ----------------------------------------------------------------------------
/**
   Lowpass Filter
*/
class TLowPass:public TFilter{
private:
  // time-domain resolution
  double tdres;
  // cut-off frequency
  double Fc;
  // parameters used in calculation
  double c1LP,c2LP;
  // gain of the filter
  double gain;
  // order of the low-pass filter
  int LPorder;
  double hc[MAX_ORDER],hcl[MAX_ORDER];
public:
  /// Get the order of the filter
  int GetOrder(void){return LPorder;};
  /// Construct the filter with time_resolustion(1/F_s), cutoff frequency, gain, filter order
  TLowPass(double _tdres,double _Fc,double _gain,int _LPorder)
  { init(_tdres,_Fc,_gain,_LPorder); };
  /// Construct the filter with time_resolustion(1/F_s), cutoff frequency, gain=1, filter order
  TLowPass(double _tdres,double _Fc,int _LPorder)
  { init(_tdres,_Fc,1.0,_LPorder); };
  /// Construct the filter with time_resolustion(1/F_s), cutoff frequency, gain, filter order
  void init(double _tdres,double _Fc,double _gain,int _LPorder);
  /// return the cutoff frequency of the filter
  double cutoff(void){ return(Fc); };
  /// filtering
  double run(double x);
};

// ----------------------------------------------------------------------------
/**
  Highpass Filter
*/
class THighPass:public TFilter{
private:
  /*time-domain resolution, cut-off frequency, gain of the filter */
  double tdres,Fc,c1HP,c2HP,gain;
  int HPorder;       //order of the hi-pass
  double hc[MAX_ORDER],hcl[MAX_ORDER];
public:
  /// get the order of the filter
  int GetOrder(void){return HPorder;};
  /// Construct the filter with time_resolution(1/F_s), cut-off freq., gain,order
  THighPass(double _tdres,double _Fc,double _gain,int _HPorder)
  { init(_tdres,_Fc,_gain,_HPorder); };
  /// Construct the filter with time_resolution(1/F_s), cut-off freq., gain=1, order
  THighPass(double _tdres,double _Fc,int _HPorder)
  { init(_tdres,_Fc,1.0,_HPorder); };
  /// Init the filter with time_resolution(1/F_s), cut-off freq., gain,order
  void init(double _tdres,double _Fc,double _gain,int _HPorder);
  /// return the cutoff frequency of the filter
  double cutoff(void){ return(Fc); };
  /// filtering the signal
  double run(double x);
};

// ----------------------------------------------------------------------------
/**
   Gammatone filter
 */
class TGammaTone:public TFilter
{
private:
  double phase;
  /* Cutoff Freq(tau), Shift Freq, ... */
  double tdres,tau,F_shift,delta_phase,gain,c1LP,c2LP;
  COMPLEX gtf[MAX_ORDER],gtfl[MAX_ORDER];
  int gtf_order;
public:
  /// Construct the filter(time_resolution(1/F_s), Center Frequency, tau, gain, order)
  TGammaTone(double _tdres,double _Fshift,double _tau,double _gain,int _order)
  { init(_tdres,_Fshift,_tau,_gain,_order); };
  /// Construct the filter(default gain = 1)
  TGammaTone(double _tdres,double _Fshift,double _tau,int _order)
  { init(_tdres,_Fshift,_tau,1.0,_order); };
  /// Contruct the filter: do not initialize
  TGammaTone(void) {};
  /// get the order of the filter
  int GetOrder(void){ return gtf_order; };
  /// Initialize the filter(time_resolution(1/F_s), Center Frequency, tau, gain, order)
  void init(double _tdres,double _Fshift,double _tau,double _gain,int _order);
  /// Set the tau of the gammatone filter, this is useful for time-varying filter
  void settau(double _tau);
  /// get the tau of the fammatone filter
  double gettau(void){ return tau; };
  /// get the center frequency of the gammatone filter
  double getF(void){ return F_shift; };
  /// set the gain of the gammatone filter, useful for time-varying filter
  void setgain(double _gain){gain = _gain;};
  /// filtering
  double run(double x);
};
#endif
