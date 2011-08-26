#ifndef _HC_HPP
#define _HC_HPP
#include <math.h>
#include <stdlib.h>
#include "filters.hpp"

/* The input of the ihc is the output of the basilar membrane, the
 * gain of the BM is normalized as 1 at CF tone, the output of ihc
 * is the potential
 * This file contains the model of the hair cell and nonlinear
 * funcitons used in the hair cell model
 */

/** Base Class that implements the nonlinear function
 */
class TNonLinear{ /* the Maximum output the nonlinearity is normalized as 1 */
public:
  /// nonlinear output of the function
  virtual double run(double x){
    x = x*1000;
    if(x>0) x = 1.1*x;
    if (x<-0.3330) x = -0.3300;
    else if(x>1) x = 1;
    return(x);
    };
};

/** Logrithm Nonlinear Function\\
    This nonlinearity is used in the IHC to fine tune the AC/DC relationship of the
    IHC output, the output is not saturated
 */
class TLogrithm:public TNonLinear{
private:
  double p_corner,p_slope,p_strength; //nonlinear transform before Boltzman function
  double asym;
  double x0,x1,s0,s1,bias,shift;      //for second order boltzman function
public:
  TLogrithm(double _corner, double _slope, double _strength,double _asym)
  {
    p_corner = _corner;
    p_slope = _slope;
    //p_strength = _strength;
    p_strength = 20.0e6/pow(10,p_corner/20);
    asym = _asym;
  };
  // Could consider just using y = log(1+B*x)
  double run(double x)
  {
    double splx,xx;
    double B = p_strength;
    xx = log(1+B*fabs(x))*p_slope;
    if(x<0)
      {
	//xx = -0.3*xx;
	splx = 20*log10(-x/20e-6);
	//double asym_t = asym-(asym-1)/(1+exp(splx/10.0-2));
	double asym_t = asym-(asym-1)/(1+exp(splx/5.0));
	//double asym_t = asym-(asym-1)/(1+exp((splx-5)/5.0));
	//double asym_t = asym-(asym-1)/(1+exp((splx-10)/8.0));
	//double asym_t = asym;
	xx = -1/asym_t*xx;
      };
    //    xx = x>0?x:0.3*x;
    return(xx);
  };
};

/**  Logrithm Nonliear followed by Second Order Boltzman Nonlinear Function\\
     This nonlinear function is used in the OHC model\\
     The output is normalized to maximum of 1
 */
class TBoltzman:public TNonLinear{
private:
  double p_corner,p_slope,p_strength; //nonlinear transform before Boltzman function
  double x0,x1,s0,s1,bias,shift;      //for second order boltzman function
public:
  /// Construct the nonlinear function
  TBoltzman(double _corner, double _slope, double _strength,
	    double _x0, double _s0, double _x1, double _s1)
  {
    p_corner = _corner;
    p_slope = _slope;
    p_strength = _strength;
    x0 = _x0;
    s0 = _s0;
    x1 = _x1;
    s1 = _s1;
    shift = 1.0/(1.0+exp(x0/s0)*(1+exp(x1/s1)));
  };
  /// Construct the nonlinear function using determined asymetry(by modifying x0)
  TBoltzman(double _corner, double _slope, double _strength,
	    double _x0, double _s0, double _x1, double _s1, double _asym)
    /* asym is the ratio of positive Max to negative Max*/
  {
    p_corner = _corner;
    p_slope = _slope;
    p_strength = _strength;
    //x0 = _x0;
    s0 = _s0;
    x1 = _x1;
    s1 = _s1;
    shift = 1.0/(1.0+_asym);
    x0 = s0*log((1.0/shift-1)/(1+exp(x1/s1)));
  };
  /// get the output of the first nonlinear function
  double displacement(double x)
  {
    double xx,splx;
    if(x<0)
      {
	splx = 20*log10(-x/20e-6);
	xx = -p_slope*log(1+exp((splx-p_corner)*p_strength))/p_strength;
      }
    else if(x>0)
      {
	splx = 20*log10(x/20e-6);
	xx = p_slope*log(1+exp((splx-p_corner)*p_strength))/p_strength;
      }
    else xx = 0;
    return(xx);
  };
  /// get the output of the second nonlinear function(Boltzman Function)
  double boltzman(double x)
  {
    double out;
    out = 1.0/(1.0+exp(-(x-x0)/s0)*(1.0+exp(-(x-x1)/s1)))-shift;
    return(out/(1-shift));
  };
  /// output of the nonlinear function, the output is normalized with maximum value of 1
  double run(double x)
  {
    double x1 = displacement(x);
    double x2 = boltzman(x1);
    //double x2 = (x1>0)?x1:0; //this is just a test
    return(x2);
  };

};

/** A nonlinear function used in the control path after the ohc

    The output controls the tau of the time-varing filter directly, when
    if input is zero the output is TauMax(sharp tuning),if the output is
    saturated around(dc) the output is TauMin(broad tuning)\\
    The nonlinear function : out = TauMax*(minR+(1.0-minR)*exp(-x/s0));\\
    this nonlinearity is determined by the relationship between gain and tauMin/tauMax
 */
class TNL_AfterOHC:public TNonLinear{
private:
  double dc,minR,s0,A;
  double TauMax,TauMin;
public:
  /// Construct the nonlinearity
  TNL_AfterOHC(void){};
  /// Construct the nonlinearity(make the output change from TauMax to TauMin approximately
  TNL_AfterOHC(double _TauMin,double _TauMax)
  { init(_TauMin,_TauMax,0.28,0.1); };
  /** Construct the nonlinearity(make the output change from TauMax to TauMin when the \\
      input changes from 0~dc, and the minimum value of tau is determined by minR*TauMax
   */
  TNL_AfterOHC(double _TauMin,double _TauMax,double _dc,double _minR)
  { init(_TauMin,_TauMax,_dc,_minR); };
  /// called by the constructor to init the nonlinearity
  double init(double _TauMin,double _TauMax,double _dc,double _minR)
  {
    TauMax = _TauMax;
    TauMin = _TauMin;
    double R = _TauMin/_TauMax;
    if(R<_minR) minR = 0.5*R;
    else minR   = _minR;
    A      = minR/(1-minR); /* makes x = 0; output = 1; */
    dc     = _dc;
    R      = R-minR;
    // This is for new nonlinearity
    s0     = -dc/log(R/(1-minR));
  };
  /** output of the nonliearity

      out = TauMax*(minR+(1.0-minR)*exp(-x1/s0));\\
      if the input is zero, the output is TauMax,\\
      if the input is dc, the output is TauMin,\\
      if the input is too large, the output is pushed to TauMax*minR
      if the input is negative, the output is the 2*TauMax-out (keep
      the control signal continuous)
   */
  double run(double x)
  {
    double out;
    double x1 = fabs(x);
    out = TauMax*(minR+(1.0-minR)*exp(-x1/s0));
    return(out);
  };

};
/** Hair Cell CLASS

    the basic strcution of the hc is the nonlinearity followed by low pass filter
    For inner hair cell model, it is logrithm nonlinearity with 7th order low pass filter,\\
    for outer hair cell model, it is two-nonlinear function followed by 5th order low pass filter\\
    for detail, see the cmpa.cpp
 */
class THairCell{
private:
  TLowPass *hclp;
  TNonLinear *hcnl;
public:
  /// construct the hair cell model, determine the low pass filter
  THairCell(double tdres,double Fc,double LPorder)
  {
    hclp = new TLowPass(tdres,Fc,LPorder);
    hcnl = new TNonLinear();
  };
  /// destruct the model
  ~THairCell(void) {delete hclp; delete hcnl;};
public:
  /// set the nonlinearity in the hair cell model is logrithm followed by boltzman (OHC)
  void init_boltzman(double corner,double slope, double strength,
			    double x0,double s0,double x1,double s1)
  {
    delete hcnl;
    hcnl = (TNonLinear*)(new TBoltzman(corner,slope,strength,x0,s0,x1,s1));
  };
  /// set the nonlinearity in the hair cell model is logrithm followed by boltzman (OHC)
  void init_boltzman(double corner,double slope, double strength,
			    double x0,double s0,double x1,double s1,double asym)
  {
    delete hcnl;
    hcnl = (TNonLinear*)(new TBoltzman(corner,slope,strength,x0,s0,x1,s1,asym));
  };
  /// set the nonlineaity in the hair cell model
  void init_logrithm(double corner,double slope,double strength,double asym)
  {
    delete hcnl;
    hcnl = (TNonLinear*)(new TLogrithm(corner,slope,strength,asym));
  };
  /// return the cutoff frequency of the lowpass fiter
  double cutoff(void){return hclp->cutoff();};
  /// get the output of the hair cell model
  double run(double x);
};

#endif

