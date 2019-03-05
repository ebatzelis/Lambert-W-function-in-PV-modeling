//****************************************************************************************
// Six implementation alternative of the Lambert W function are given here
// Version: 1.0
// Date: 5 Mar 2019
// Developer: Stratis Batzelis
//****************************************************************************************

#include "LambertWs.h"

//***************************************************************************
//                               Routines
//***************************************************************************

double lambertWmat(double x, int* ITER)
{// Lambert W function W{x}, for any  real x
// Iterative Haley's method (MATLAB function)
// Convergence tolerance set to 1e-6, Max iterations set to 10
// ITER holds the number of iterations needed for monitoring purposes

    double w,newW,EXP,F;
    double i;

    // Init
    w = log(x+1); // This is a simple init strategy. MATLAB uses something much more complicated and cost-intensive

    // Iterate
    for(i=0;i<10;i++) {
        EXP = exp(w);
        F = w*EXP -x;
        // Step
        newW = w -F/( EXP*(w+1) -(w+2)*F/(2*w+2) );
        // Break
        if ( fabs(newW-w) <= 1e-6*fabs(w) ) {
            w = newW;
            break;
        }
        w = newW;
    }
    *ITER = i+1;

    return w;
}

double lambertWasymp7(double a, double b)
{// Lambert W function W{x} = W{a*exp(b)}, where x>=3
// Asymptotic formula - 7 terms

    double L1,L2,one_div_L1;

    L1 = log(a)+b; // log
    if (L1 > LN3) {
        L2 = log(L1); // log
        one_div_L1 = 1/L1; // div

        return L1 -L2 +L2*one_div_L1*( 1 +(-2+L2 +(6-9*L2+2*L2*L2
                +( -12+36*L2-22*L2*L2+3*L2*L2*L2
                +(60-300*L2+350*L2*L2-125*L2*L2*L2+12*L2*L2*L2*L2)*0.2*one_div_L1
                )*0.5*one_div_L1)*ONE_THIRD*one_div_L1 )*0.5*one_div_L1);
    }
    else
        return NAN;
}

double lambertWasymp4(double a, double b)
{// Lambert W function W{x} = W{a*exp(b)}, where x>=3
// Asymptotic formula - 4 terms

    double L1,L2,one_div_L1;

    L1 = log(a)+b; // log
    if (L1 > LN3) {
        L2 = log(L1); // log
        one_div_L1 = 1/L1; // div

        return L1 -L2 +L2*one_div_L1*( 1 + (-2+L2)*0.5*one_div_L1);
    }
    else
        return NAN;
}

double lambertWhybrid(double a, double b)
{// Lambert W function W{x} = W{a*exp(b)}, where x>=0
// Hybrid calculation formula

    double L1,L2,one_div_L1;
    double u,p,r,r2;

    L1 = log(a)+b; // log
    if (L1 >= LN9) { // x >= 9 - Asymptotic expansion
        L2 = log(L1); // log
        one_div_L1 = 1/L1; // div

        return L1 -L2 +L2*one_div_L1*( 1 +(-2+L2 +(6-9*L2+2*L2*L2
                +( -12+36*L2-22*L2*L2+3*L2*L2*L2
                +(60-300*L2+350*L2*L2-125*L2*L2*L2+12*L2*L2*L2*L2)*0.2*one_div_L1
                )*0.5*one_div_L1)*ONE_THIRD*one_div_L1 )*0.5*one_div_L1);
    }
    else { // 0 <= x < 9 - Series expansion
        u = exp(L1-1); // exp
        p = 1-u;
        r = 1/(1+u); // div
        r2 = r*r;

        return u + u*r*p*(1
                +p*0.5*r2*(1
                + p*ONE_THIRD*r2*(-2*u+1
                +p*0.25*r2*(6*u*u-8*u+1
                +p*0.2*r2*(-24*u*u*u+58*u*u-22*u+1) ) ) ) );
    }
}

double lambertWsimple(double a, double b)
{// Lambert W function W{x} = W{a*exp(b)}, where x>=2
// Simple approximation formula

    double L;

    L = log(a)+b; // log
    if (L > LN2) {
        return L -L*log(L)/(L+1); // log+div
    }
    else
        return NAN;
}

double lambertWanalyt(double x)
{// Lambert W function W{x}, where x>=3e-3
// Analytical approximation formula

    const double E = 0.4586887;

    if (x > 3e-3) {
        return (1+E)*log(1.2*x/log(2.4*x/log(1+2.4*x)) ) -E*log(2*x/log(1+2*x)); // 5 log +3divs
    }
    else
        return NAN;
}



