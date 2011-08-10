/** \file cmpa.cpp
 * \brief This file contains several useful routines 
 * TAuditoryNerve::construct : construct the fiber model
 * TAuditoryNerve::init      : initialize the parameters of the model (These two must be called)
 * TAuditoryNerve::run       : pass the signal through the Model (not the spikes)

 * synapse() synapse2()      : generate spikes using output of the model
 * Get_tau()                 : How the tau of the time-varing tuning filter is determined
 * TBasilarMembrane::run()   : pass the signal through the tuning filter
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>

#include "cmpa.hpp"
#include "complex.hpp"
#include "filters.hpp"
#include "hc.hpp"
#include "synapse.hpp"


// -------------------------------------------------------------------------
/**
 * run the model\n
 * Input signal->(time-varying filter+1st-order gamma-tone filter)
 * ->inner hair cell->synapse model
 */
double TAuditoryNerve::run(double in)
{ //inbuf,tau,bmbuf,ihcbuf,sbuf are variables in the CLASS
  //bm,gfagain,ihc,syn are variables(CLASS pointer) in the CLASS
  inbuf		= in;
  double bmbuf1	= bm->run(in);
  tau    	= bm->GetTau();
  bmbuf 	= gfagain->run(bmbuf1);
  ihcbuf 	= ihc->run(bmbuf);
  sbuf   	= syn->run(ihcbuf);
  return(sbuf);
};

// --------------------------------------------------------------------------
/** construct the default nerve, alloc all the memory and so on,
   the tdres,mode should be specified, the cf independent parameters should be determined
   in this function, others are specified in init(...) function
*/
void TAuditoryNerve::construct(void)
{ //default
  bm = new TBasilarMembrane(tdres);
  bm->wbfilter = new TGammaTone;
  bm->afterohc = new TNL_AfterOHC;

  /* *** FOLLOWING PARAMETERS ARE NOT CF,SPONT,or SPECIES DEPENDENT *** */
  double absdb = 20; /* The value that the BM starts compression */
  /* ********************* Init OHC model *******************************/
  //bm->ohc = new THairCell(tdres,800.0,5);
  bm->ohc = new THairCell(tdres,800.0,3);

  /* parameter into boltzman is corner,slope,strength,x0,s0,x1,s1,asym */
  /* The corner determines the level that BM have compressive nonlinearity*/
  ohcasym = 7.0; //CLASS variable
  bm->ohc->init_boltzman(absdb-12,0.22,0.08,5,12,5,5,ohcasym);

  /* ******** Init the gammatone filter after the bm tuning filter *******/
  gfagain = new TGammaTone;
  /* ********************* Init Inner Hair Cell *************************/
  ihc = new THairCell(tdres,3800.0,7);
  //Using new nonlinear function p_corner,p_slope, p_str, asym
  ihc->init_logrithm(80,0.1,0.12,3);
  /********************** Init Synapse *********************************/
  syn = new TSynapse_WS(tdres);
  /********************** Init Spike Generator *************************/
  /* Spike history parameters- for two exponents ala Westerman and Smith '89
     these values fit Gaumond & Kim's spontaneous, and also fit noise responses. */
  /* parameter into SG is dead,c0,s0,c1,s1 in units of sec */
  sg = new TSpikeGenerator(tdres,0.00075,0.5,0.001,0.5,0.0125);
};

// --------------------------------------------------------------------------------
void TAuditoryNerve::init(double _cf,double _spont)
{
  int i;
  cf = _cf;
  if(_spont<0) spont = 50;
  else spont = _spont;

  int bmorder = 3;
  /* set up the tau0 and taurange of time-varing gammatone filter */
  Get_tau(cf,bmorder,taumax,tau0,taurange);
  /* *** FOLLOWING PARAMETERS ARE CF,SPONT,or SPECIES DEPENDENT *** */
  bm->init(cf,taumax,tau0,bmorder);
  gfagain->init(tdres,cf,tau0,  1.0,1); //non-time-varying gammatone filter that follows time-varying filter
  /* ******************** Init WideBand Filter ****************************/
  int species = 1; //for CAT
  double x = cochlea_f2x(species,cf);
  double centerfreq = cochlea_x2f(species,x+1.2); //shift the center freq
  double wb_order = 3;
  double tauwb = tau0+0.2*(taumax-tau0);
  bm->TauWB = tauwb;
  bm->TauWBMin = tauwb/taumax*tau0;
  double dtmp = tauwb*TWOPI*(centerfreq-cf);
  double wb_gain = pow((1+dtmp*dtmp), wb_order/2.0);
  bm->wbfilter->init(tdres,centerfreq,tauwb,wb_gain,wb_order);

  /* ********************* Init Nonlinear After OHC *********************/
  /* parameter into nonlinear function is Taumin,taumax,dc,minR        */
  double dc = (ohcasym-1)/(ohcasym+1.0)/2.0;
  bm->afterohc->init(tau0,taumax,dc-0.05,0.05);
  /* ********************* Init Synapse ********************************/
  double PImax = 0.6;
  double cf_factor = 2+3*log10(cf/1000);
  if (cf_factor<1.5) cf_factor = 1.5;
  double kppi = (1+spont)/(5+spont)*cf_factor*20*PImax;
  double vsat = kppi;                           //in fact vsat = kppi+PI1
  syn->setppi_mgheinz(kppi);   // set the IHC-PPI nonlinearity as soft-rectifier
  /* parameter into synapse is Ass,Asp,TauR,TauST,Ar_Ast,PI2,vsat */
  syn->init_WS(350.0,spont,2e-3,60e-3,6,0.6,vsat); //vsat is not used now until kppi is set to negtive
};

// --------------------------------------------------------------------------------
/** Get TauMax, TauMin for the tuning filter. The TauMax is determined by the bandwidth/Q10
    of the tuning filter at low level. The TauMin is determined by the gain change between high
    and low level
    Also the calculation is diffrent for different species
 */
double Get_tau(double cf,int order, double &taumax,double &taumin,double &taurange)
{
  double Q10,bw,gain,ratio;
  double xcf, x1000;
  gain = 20+42*log10(cf/1e3);                // estimate compression gain of the filter
  if(gain>70) gain = 70;
  if(gain<15) gain = 15;
  ratio = pow(10,(-gain/(20.0*order))); // ratio of TauMin/TauMax according to the gain, order

  // Calculate the TauMax according to different species
  /* Universal species from data fitting : From Xuedong Zhang & Ian Bruce */
  /* the value of Q10 determines the taumax(bandwidths at low level) Based on Cat*/
  Q10 = pow(10,0.4708*log10(cf/1e3)+0.4664);
  bw = cf/Q10;
  taumax = 2.0/(TWOPI*bw);
  taumin =  taumax*ratio;
  taurange = taumax-taumin;
  return 0;
};

// --------------------------------------------------------------------------------
/** Calculate the location on Basilar Membrane from best frequency
*/
double cochlea_f2x(int species,double f)
{
  double x;
  switch(species)
    {
    default:
    case 1: //cat
      x = 11.9 * log10(0.80 + f / 456.0);
      break;
    };
  return(x);
};
// --------------------------------------------------------------------------------
/** Calculate the best frequency from the location on basilar membrane
 */
double cochlea_x2f(int species,double x)
{
  double f;
  switch(species)
    {
     default:
    case 1: //cat
      f = 456.0*(pow(10,x/11.9)-0.80);
      break;
    };
  return(f);
};
// --------------------------------------------------------------------------------
/** error: print an error message and die gracefully
    Takes arguments like printf
    Copied from Kernighan and Ritchie, p 174
*/
void error(char *fmt)
{
  fprintf(stderr, "error: ");
  fprintf(stderr, fmt);
  fprintf(stderr, "\n");
  exit(1);  /* closes all open files */
}
// ---------------------------------------------------------------------------
/** Calculate the delay(basilar membrane, synapse for cat*/
double delay_cat(double cf)
{
  /* DELAY THE WAVEFORM (delay buf1, tauf, ihc for display purposes)  */
  /* Note: Latency vs. CF for click responses is available for Cat only (not human) */
  /* Use original fit for Tl (latency vs. CF in msec) from Carney & Yin '88
     and then correct by .75 cycles to go from PEAK delay to ONSET delay */
  double A0 = 8.13; /* from Carney and Yin '88 */
  double A1 = 6.49;
  double x = cochlea_f2x(1,cf); //cat mapping
  double delay = A0 * exp( -x/A1 ) * 1e-3 - 1.0/cf;
return(delay);
};
// ----------------------------------------------------------------------------------------------
/** pass the signal through the tuning filter.
    using different model
    Sharp_Linear | Broad_Linear | Broad_Linear_High | FeedBack_NL | FeedForward_NL
*/
double TBasilarMembrane::run(double x)
{
  double out;
  double x1 = bmfilter->run(x); // pass the signal through the tuning filter

  double wbout;
  wbout = wbfilter->run(x);    //get the output of the wide-band pass as the control signal
  static double A=(TauWB/TauMax-TauWBMin/TauMin)/(TauMax-TauMin);
  static double B=(TauMin*TauMin*TauWB-TauMax*TauMax*TauWBMin)/(TauMax*TauMin*(TauMax-TauMin));
  double taunow = A*tau*tau-B*tau;
  wbfilter->settau(taunow);        //set the tau of the wide-band filter
  // normalize the gain of the wideband pass filter as 0dB at CF
  double dtmp = taunow*TWOPI*(wbfilter->getF()-bmfilter->getF());
  double wb_gain = pow((1+dtmp*dtmp), wbfilter->GetOrder()/2.0);
  wbfilter->setgain(wb_gain);
  double ohcout = ohc->run(wbout); // pass the control signal through OHC model
  tau = afterohc->run(ohcout);     // pass the control signal through nonliearity after OHC
  bmfilter->settau(tau);           // set the tau of the tuning filter
  out = pow((tau/TauMax),bmfilter->GetOrder())*x1;  // Gain Control of the tuning filter
  return(out);
};
