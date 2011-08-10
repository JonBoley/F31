#ifndef _CMPA_HPP
#define _CMPA_HPP
/** \file cmpa.hpp
 *  \brief This file contain two classes and other most useful routins used by the model
 *
 *   Almost all the parameters of the model are determined in cmpa.cpp
 */
/**
 * \structure structure TBasilarMembrane
 * \structure structure TAuditoryNerve
 */
#include "filters.hpp"
#include "hc.hpp"
#include "synapse.hpp"
#include "complex.hpp"

#define TWOPI 6.2831853
/// maximum order of the filter
#define MAX_ORDER 10

class TBasilarMembrane;
class TAuditoryNerve;

/* Get the parameters of bm */
/// Set the time constant/bandwidth of the tuning filter
double Get_tau(double cf, int order, double &taumax,double &taumin,double &taurange);
/// mapping from CF to basilar memebrane position
double cochlea_f2x(int species,double f);
/// mapping from BM position to CF
double cochlea_x2f(int species,double x);
/// delay of the cat
double delay_cat(double cf);
/// deal with the fatal error
void error(char *fmt);

/** The class the define the basic structure of the time-varing filter

    the class consists of following components: \n
    1. tuning filter(bmfilter), the tau is controlled by the control path\n
    2. wideband pass filter(wbfilter)\n
    3. outer hair cell model(ohc)\n
    4. nonlinear function after the outer haircell(afterohc)
 */
class TBasilarMembrane{ /* class of basilar membrane */
public:
  double tdres; //time domain resolution
  double tau,TauMax,TauMin; //time constant of signal path filter
  double TauWB,TauWBMin;  //time constant of control path filter
  TGammaTone *bmfilter, *wbfilter; //signal and control path filter
  THairCell *ohc; //OHC-like model (nonlinear function + low pass filter)
  TNL_AfterOHC *afterohc; //nonlinear mapping after the ohc-like model
public:
  /// get the tau of the tuning filter
  double GetTau(void){ return tau;};
  /// construct the basilar membrane tuning model
  TBasilarMembrane(double _tdres)
  {
    wbfilter = NULL;
    /* <a name="ihcppirun"></a> */
    ohc = NULL;
    afterohc = NULL;
    tdres = _tdres;
    bmfilter = new TGammaTone;
  };
  /// destruct the basilar membrane tuning model
  ~TBasilarMembrane(void)
  {
    if(wbfilter!=NULL) delete wbfilter;
    if(ohc!=NULL) delete ohc;
    if(afterohc!=NULL) delete afterohc;
    if(bmfilter!=NULL) delete bmfilter;
  };
  /// init the basilar membrane tuning-filter
  void init(double cf,double _TauMax,double _TauMin,int _order)
  {
    TauMax = _TauMax;
    TauMin = _TauMin;
    tau = TauMax;
    bmfilter->init(tdres,cf,tau,1.0,_order); /* 1 is the gain of the filter */
  };
  ///run throught the tuning filter
  virtual double run(double x);
};

/** Class of the auditory nerve fiber, this is a complete model of the fiber

    The class consists of all the parts of the auditory nerve fiber\n
    1. time-varying tuning filter with control path: TBasilarMembrane(bm, 3rd order)\n
    2. the gammatone filter after the 3rd-order nonlinear filter(gffilter)\n
    3. inner hair cell model(ihc)\n
    4. synapse model, from Westman,1986(syn)\n
    5. spike generator, from Carney,1993(sg)
 */
class TAuditoryNerve{
public:
  double tdres;            //time domain resolution
  double cf;               //cf of the auditory nerve
  double spont;            //spontaneous rate of the fiber
  TBasilarMembrane *bm;    //basilar membrane nonlinear fiber

  TGammaTone *gfagain;	   //Linear filter after the nonlinear filter
  THairCell *ihc;          //inner hair cell model, in <hc.hpp>
  TSynapse_WS *syn;        //synapse model, in <synapse.hpp>
  TSpikeGenerator *sg;     //spike generator, in <synapse.hpp>
public:
  /* these are parameters used in the model */
  double ohcasym;
  double tau0,taumax,taurange;
  double inbuf,tau,bmbuf,ihcbuf,sbuf; //output of different stage
public:
  TAuditoryNerve(void)
  { bm = NULL, ihc = NULL, syn = NULL, sg = NULL; };
  double run(double x);
  /// construct the model
  void construct(void);
  /// init the model (set the value of the parameters)
  void init(double _cf,double _spont);
  ///deconstruct the fiber model
  ~TAuditoryNerve(void)
  {
    if(bm!=NULL) delete bm;
    if(ihc!=NULL) delete ihc;
    if(syn!=NULL) delete syn;
    if(sg!=NULL) delete sg;
  };
};

#endif
