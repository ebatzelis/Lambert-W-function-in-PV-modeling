//****************************************************************************************
// This header contains definitions and prototypes of the LambertWs.c
// Version: 1.0
// Date: 5 Mar 2019
// Developer: Stratis Batzelis
//****************************************************************************************

#ifndef LAMBERTWS_H_
#define LAMBERTWS_H_

#include "Math.h"

#define LN2 0.69314718055994530942
#define LN3 1.09861228866810956006
#define LN9 2.19722457733621956422
#define ONE_THIRD 0.33333333333333333333

//////////////////////////////////////////////////////////////////////////////////////////
// Functions prototypes
//////////////////////////////////////////////////////////////////////////////////////////
double lambertWmat(double x, int* ITER);
double lambertWasymp7(double a, double b);
double lambertWasymp4(double a, double b);
double lambertWhybrid(double a, double b);
double lambertWsimple(double a, double b);
double lambertWanalyt(double x);

#endif /* LAMBERTWS_H_ */
