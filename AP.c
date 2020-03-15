/** @file
 *  @brief It is an intermediate layer between the main programme, which works with an abstract electrophysiological model, and specific electro-models.
 *
 *  For each model, we define its header file, some constants needed my the main program, initialization and
 *  main calculation functions.
 */

#include <math.h>
#include "elektro.h"

#ifdef AP

const double volt_init=0;
const double InitialV=0;
const double volt_stim=1;

#ifdef curr_stim
double
    I_stim, t_stim, tipval=0.5;
#endif

double ka;
double a;
double eta;

const int indV=0;
/// eps(u) in the Aliev--Panfilov model (it was ((u)<0.05) ? 1.0 : 0.1 )
#define epsAP(u) (((u)<a) ? 1.0 : eta)

void CommonInit(double dt)
{

}

void InitUV(double uv[Nvar+1])
{
    uv[0]=volt_init;
    uv[1]=InitialV;
}

/// right parts of the system, reaction only
double rightU(double u, double v)
{
    return -ka*u*(u-a)*(u-1)-u*v;
}

double rightV(double u, double v)
{
    return epsAP(u)*(ka*u-v);
}
#endif /// AP

#ifdef APmu

const double volt_init=0; /// steady-state potential
const double InitialV=0;
const double volt_stim=1;

#ifdef curr_stim
double
    I_stim, t_stim, tipval=0.5; /// ???
#endif

double ka; //ka; a are constants of AP model
double a;
double epsAP_v;
double mu1;
double mu2;

#ifdef APmuRusakov
const double b=0.1;
#endif

const int indV=0;

void CommonInit(double dt)
{

}

void InitUV(double uv[Nvar+1])
{
    uv[0]=volt_init;
    uv[1]=InitialV;
}

/// right parts of the system, reaction only
double rightU(double u, double v)
{
    return -ka*u*(u-a)*(u-1)-u*v;
}

double rightV(double u, double v)
{
    return -(epsAP_v+mu1*v/(u+mu2))*(v+ka*u*(u-
    #ifdef APmuRusakov
        b
    #else
        a
    #endif
        -1));
}
#endif ///APmu


#if defined(AP) || defined(APmu)
void CalcNewUV(double uv0[Nvar+1], double uv1[Nvar+1], double DivDGradU, double dt
#ifdef curr_stim
                ,double Istim
#endif
               )
{
    uv1[0] = uv0[0] + dt * (DivDGradU + rightU(uv0[0],uv0[1])
#ifdef curr_stim
                            - Istim
#endif
                            );
    uv1[1] = uv0[1] + dt * rightV(uv0[0],uv0[1]);
}
#endif /// AP, APmu
