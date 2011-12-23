#include <iostream.h>
#include "synapse.hpp"

/**

   This function initialize the synapse model, derive the paramters needed in the model\\
   from WS 1988's paper
 */
void TSynapse_WS::init_WS(double _Ass,double _Asp,double _TauR,double _TauST,
   double _Ar_Ast,double _PI2,double _vsat)
{ /* The following parameters are derived from the equations in Westerman,1988*/
  /* The parameters are first set to describe the characteristic response
     equation of PSTH (eq.10) */
  Ass = _Ass;              /* Steady State Firing Rate eq.10 */
  Asp = _Asp;              /* Spontaneous Firing Rate  eq.10 */
  TauR = _TauR;            /* Rapid Time Constant eq.10 */
  TauST = _TauST;          /* Short Time Constant eq.10 */
  Ar_Ast = _Ar_Ast;        /* Ratio of Ar/Ast */
  PI2 = _PI2;              /* PI2 : Maximum of the PI(PI at steady state) */
  PTS = 1+9*Asp/(9+Asp);   /* Peak to Steady State Ratio, characteristic of PSTH */

  set_spont(Asp);
  
  //now get the other parameters
  Aon = PTS*Ass;                           /* Onset rate = Ass+Ar+Ast eq.10 */
  double AR = (Aon-Ass)*Ar_Ast/(1+Ar_Ast); /* Rapid component magnitude: eq.10 */
  double AST = Aon-Ass-AR;                 /* Short time component: eq.10 */
  double PI1 = PI2/Aon*Asp;                /* eq.A15 */
  CG = (Asp*(Aon-Asp))/(Aon*PI1*(1-Asp/Ass)); /* eq.A16 */
  double gamma1 = CG/Asp;                  /* eq.A19 */
  double gamma2 = CG/Ass;                  /* eq.A20 */
  double k1 = -1/TauR;                     /* eq.8 & eq.10 */
  double k2 = -1/TauST;                    /* eq.8 & eq.10 */
  /* eq.A21 & eq.A22 */
  double VI0 = (1-PI2/PI1)/(gamma1*(AR*(k1-k2)/CG/PI2+k2/PI1/gamma1-k2/PI2/gamma2));
  double VI1 = (1-PI2/PI1)/(gamma1*(AST*(k2-k1)/CG/PI2+k1/PI1/gamma1-k1/PI2/gamma2));
  VI = (VI0+VI1)/2;
  double alpha = gamma2/k1/k2;             /* eq.A23,eq.A24 or eq.7 */
  double beta = -(k1+k2)*alpha;            /* eq.A23 or eq.7 */
  double theta1 = alpha*PI2/VI;
  double theta2 = VI/PI2;
  double theta3 = gamma2-1/PI2;
  /* from eq.3-5 we get:
     1/PG+1/PL = theta3;                         (eq.5')
     VL = theta1*PG*PL;                          (eq.3')
     beta = theta1*(PL/PI2+1)+theta2*(1/PG+1/PL) (eq.4')
  */
  PL = ((beta-theta2*theta3)/theta1-1)*PI2; /* eq.4' */
  PG = 1/(theta3-1/PL);                     /* eq.5' */
  VL = theta1*PL*PG;                        /* eq.3' */
  CI = Asp/PI1;                             /* CI at rest, from eq.A3,eq.A12 */
  CL = CI*(PI1+PL)/PL;                      /* CL at rest, from eq.1 */
  Prest = PI1;
  PImax = PI2;
  vsat = _vsat;
  //init the ihcppi transform
  ihcppi->init(PImax,Prest,vsat);
  return;
};
/**

   run the synapse model: if the time resolution is not small enough, the concentration of\\
   the immediate pool could be as low as negative, at this time there is an alert message\\
   print out and the concentration is set at saturated level
 */
double TSynapse_WS::run(double x)
{
   double PPI = ihcppi->run(x);
   double CIlast = CI;
   CI = CI + (ttdres/VI)*(-PPI*CI + PL*(CL-CI));
   CL = CL + (ttdres/VL)*(-PL*(CL - CIlast) + PG*(CG - CL));
   if(CI<0){
              double temp = 1/PG+1/PL+1/PPI;
              CI = CG/(PPI*temp);
              cout<<"The immediate concentration in synapse is too low";
              CL = CI*(PPI+PL)/PL;
           };
   return(CI*PPI);
};
