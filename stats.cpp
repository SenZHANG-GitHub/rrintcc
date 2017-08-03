
#include <cmath>
#include "stats.h"

#define FPMIN 1.0e-30

double normdist(double z)
{
 double sqrt2pi = 2.50662827463;
 double t0, z1, p0 ;
 t0 = 1 / (1 + 0.2316419 * fabs(z));
 z1 = exp(-0.5 * z*z ) / sqrt2pi;
 p0 = z1 * t0
    * (0.31938153 +
    t0 * (-0.356563782 +
    t0 * (1.781477937 +
    t0 * (-1.821255978 +
    1.330274429 * t0))));
 return z >= 0 ? 1 - p0 : p0 ;
}
