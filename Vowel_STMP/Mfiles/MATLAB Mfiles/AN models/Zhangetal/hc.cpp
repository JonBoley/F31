#include "hc.hpp"
/**

Get the output of the hair cell\\
1. pass the signal through he nonlinear function\\
2. low-pass filter the nonlinear output
*/
double THairCell::run(double x)
{
  double x1 = hcnl->run(x);
  double x2 = hclp->run(x1);
  return(x2);
};

