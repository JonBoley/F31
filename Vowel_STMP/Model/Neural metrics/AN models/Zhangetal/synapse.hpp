#ifndef _SYNAPSE_HPP
#define _SYNAPSE_HPP
#include <math.h>
#include <stdlib.h>
#include <iostream.h>
#include <time.h>
/** \file synapse.hpp
 * \brief model of synapse model, ihc-ppi nonlinear, and spike generator.
 * This file contains routins that converts the ihc potential to spikes trains
 * The class TIhcPpi converts the IHC potential(normalized between[-1/3,1]) to the
 * immediate permibility of the synapse(PPI): it has three derived classes which
 * corespond to the Linear_Rectify(TIhcPpi_linear, Soft_Rectify(TIhcPpi_mgheinz),
 * and Schoohoven's method(TIhcPpi_ian)
 * This class is included by the class Tsynapse to clean up the ppi function in
 * Tsynapse and it should take parameters from TSynapse(PImax?,PIrest?) to initialize

 * The class TSynpase simulate the dynamic process of the synapse, the input is the
 * ihc, and the TIHCPPI converts ihc to ppi(immediate permibility), the output is
 * the possibility of firing rate(before refractory)
 * The derived class TSynapse_WS is a discrete implementation of Westman's model

 * The class TSpikeGenerator get the firing possibility as input and output the
 * spikes times(1). It is a random generator including the refractoriness
 */
class TIHCPPI;
typedef class TIHCPPI TIHCPPI_linear;
class TIHCPPI_mgheinz;
class TIHCPPI_ian;

class TSynapse;
class TSynapse_WS;

class SpikeGenerator;

/** Option 1: Linear equation between ihc (in range [-1/3,1]) and
 * PPI, the half-wave rectify PPI
 * when ihc input = 0, output = pimin
 * when ihc input = 1, output = vsat
 */
class TIHCPPI{
public:
  double cf,spont; /* cf and spont rate of the fiber */
  double kslope,vsat,pimax,pimin;
  virtual double run(double x)
    {
      double out = kslope*x+pimin;
      if(out<0.0) out = 0.0;
      return out;
    };
  virtual void init(double _pimax,double _pimin,double _vsat)
    {
      pimax = _pimax; pimin = _pimin; vsat = _vsat;
      if(kslope<=0) kslope = vsat-pimin;
    };
  virtual void setslope(double _kslope) { kslope = _kslope; };
  TIHCPPI(){kslope = -1;};
};

/** Option 2: NonLinear relation between ihc (in range [-1/3,1])
 * and PPI, such that vsat and Prest are achieved and the PPI
 * goes to 0 smoothly for negative ihc
 */
class TIHCPPI_mgheinz: public TIHCPPI{
public:
  double p_strength,p_corner,p_slope;
  /// Nonlinear function between VIHC and PPI with formulation of y=A*log(1+exp(B*x))
  inline double run(double x)
    {
      double tmp = p_strength*x;         //y = A*log(1+exp(B*x))
      if(tmp<400) tmp = log(1+exp(tmp));
      x = p_slope/p_strength*tmp;
      return(x);
    };
  /// Init the nonlinear function, if the kslope is not set, using vsat to set the parameters
  inline void init(double _pimax,double _pimin,double _vsat)
    {
      if(kslope>=0) _vsat = kslope+_pimin;  //kslope is not set
      p_corner = 0; //dont' use
      double tmp = log(2)*_vsat/_pimin;
      if(tmp<400) p_strength = log(exp(tmp)-1);
      else p_strength = tmp;
      p_slope = _pimin/log(2)*p_strength;
    };
  /// Construct
  TIHCPPI_mgheinz(void){};
};

/** Base class of the synanpse
*/
class TSynapse{
protected:
  double spont;
  double ttdres;
public:
  TSynapse(double tdres){ ttdres = tdres; };
  virtual void init(double _Asp){set_spont(_Asp); return;};
  virtual double run(double x){return x;};
  inline void set_spont(double Asp){spont = Asp;};
  inline void set_tdres(double tdres){ttdres = tdres;};
  inline double get_spont(void) {return spont;};
  inline double get_tdres(void) {return ttdres;};
};
/** Synapse Class that implement the Westman's three pool model with some modification
    The input of synapse(output of the ihc) is rectified using one of the method mentioned above
    (now is the soft-rectifier before treated as ppi) The parameters of the synapse model
    is determined by the characteristic of the PSTH of the fiber
 */
class TSynapse_WS: public TSynapse
{
public:
  //these parameters are pre_set in the synpase model
  double Ass, Asp, TauR, TauST, Ar_Ast, PTS, PI2, Aon;
  //these parameters are derived from the WS's model and equations
  //double Aon,PI1,Ast;
  double Prest,PPI,PImax,PL,PG,CI,CL,CG,VI,VL; //adaptation model parameters
  double CIlast; //used in the differential equations;
  //these are parameters used in IHC-PPI transform
  TIHCPPI *ihcppi;
  char ihcppitype; //l:linear, m:mgheinz, i: ian
  double vsat;
public:
  /// run the synapse model
  double run(double x);
  /// Set the ihc-ppi relationship as linear rectifier
  void setppi_linear(void){
    if(ihcppi!=NULL) delete ihcppi;
    ihcppi = (TIHCPPI*)(new TIHCPPI_linear());
    setppi_linear(-1);
  };
  /// set the ihc-ppi relationship as linear rectifier with slope _ks at zero
  void setppi_linear(double _ks) {
    if(ihcppi!=NULL) delete ihcppi;
    ihcppi = (TIHCPPI*)(new TIHCPPI_linear());
    ihcppi->setslope(_ks);
  };
  /// using mgheinz's soft rectifier
  void setppi_mgheinz(void) {
    if(ihcppi!=NULL) delete ihcppi;
    ihcppi = (TIHCPPI*)(new TIHCPPI_mgheinz());
    ihcppi->setslope(-1);
  };
  /// using mgheinz's soft rectifier with slope _ks at zero
  void setppi_mgheinz(double _ks) {
    if(ihcppi!=NULL) delete ihcppi;
    ihcppi = (TIHCPPI*)(new TIHCPPI_mgheinz());
    ihcppi->setslope(_ks);
  };
  /// init the synapse model, call init_ihcppi
  void init_WS(double _Ass,double _Asp,double _TauR,double _TauST,
	       double _Ar_Ast,double _PI2,double _vsat);
  /// reset the synapse to steady state(at rest)
  void TSynapse_WS::reset(void) {
    CI = Asp/Prest;                             /* CI at rest, from eq.A3,eq.A12 */
    CL = CI*(Prest+PL)/PL;                      /* CL at rest, from eq.1 */
  };
  /// construct the synapse model
  TSynapse_WS(double _tdres):TSynapse(_tdres){
    ihcppi = NULL;
    setppi_linear();
  };
  /// construct the synapse model
  TSynapse_WS(double _tdres,double _Ass,double _Asp,double _TauR,double _TauST,
	      double _Ar_Ast,double _PI2,double _vsat):TSynapse(_tdres)
  {
    ihcppi = NULL;
    setppi_linear();
    init_WS(_Ass,_Asp,_TauR,_TauST,_Ar_Ast,_PI2,_vsat);
  };
};
/** Spike Generator

    This is the spike generator including refractoriness, the spikes are generated by
    nonhemogeneous process. The average rate is modified by the synapse output simutaneously.
    The probability of the spike generating is also affected by the synapse history as absolute
    and relative refractoriness
 */
class TSpikeGenerator{
private:
  double tdres;
  double rtime,rsptime; /* in units of sencond */
  double dead;
  double c0,s0,c1,s1;
public:
  /// Construct the spike generator
  TSpikeGenerator(double _tdres,double _dead,
		 double _c0,double _s0,double _c1,double _s1)
  {
    time_t seed;
    seed = time(NULL);
    //srand48(seed); //this is replaced by srand() for the compatibility
	srand(seed);
    c0 = _c0, s0 = _s0, c1 = _c1, s1 = _s1;
    dead = _dead;
    tdres = _tdres;
    rtime = 0.0; rsptime = 0.0;
  };
  /// reset the spike generator, clear the history of the spike generator and set time as zero
  void init(void) {rtime = rsptime = 0.0; };
  /// reset the spike generator, simulate the situation that the fiber has spontaneous rate
  void init(double spont)
    {
      if(spont>0)
	//rsptime = rtime - drand48()*1/spont;
	rsptime = rtime - rand()/(double)(RAND_MAX)*1/spont;
      else
	{
	  rsptime = rsptime - rtime; rtime = 0;
	};
    };
  /// get the time of spike generator
  double Get_rtime(void){ return rtime; };
  /// generate the spike(1) or none according to the input(synapse output).
  /**
     This function generate the spikes according synapse output(x=sout), also including the\\
     refractoriness
  */
  int run(double x)
    {
      double rint,prob;
      int out = 0;
      rtime += tdres; /* running time */
      /* interval from last spike, including 'wrap around' between trials */
      rint = rtime - rsptime;
      if(rint > dead )
	{
	  prob = x*(1.0-( c0*exp(-(rint - dead)/s0) +
			  c1*exp(-(rint - dead)/s1))) * tdres;
	  //prob = x*pan->tdres;
	  if( prob > drand48() )
	    {
	      rsptime = rtime;
	      out = 1;
	    }
	};
      return(out);
    };
};

#endif
