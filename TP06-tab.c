#include "elektro.h"

#ifdef TP06

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

const double volt_stim=50; /// stimulation potential (without current)
double
///    I_stim=-80*4, t_stim=0.4;
///    I_stim=-60*4, t_stim=0.5;
///    I_stim=-30*4, t_stim=1;
///    I_stim=-25*2, t_stim=1.5;  <-- it is used!
///    I_stim=-20*4, t_stim=2;
    I_stim, t_stim;

double
    tipval=0,
    c_gNa=1,
    c_gCaL=1,
    c_verapamil=1;
int amiodaron; /// 0, no effect; 1, low concentration 1 micromol; 3, high conc. 3 micromol
const bool debugEl=true;

/*------------------------------------------------------------------------------
FLAG TO CHOOSE BETWEEN EPICARDIAL, ENDOCARDIAL AND MIDMYOCARDIAL CELL TYPES
------------------------------------------------------------------------------*/
#define EPI 0
//#define ENDO 1
//#define MCELL 2

// #define slope07
// #define slope11
///#define slope14
///#define slope18

const int indV=0; /// index of state variable voltage
const int numTa=-1;
/*-----------------------------------------------------------------------------
  ELECTROPHYSIOLOGICAL PARAMETERS:
-----------------------------------------------------------------------------*/

double TAU_F_INACT_FACTOR;
/*double TAU_F_INACT_FACTOR =
#ifdef slope07
0.6;
#elif defined(slope11)
1.0;
#elif defined(slope14)
1.5;
#elif defined(slope18)
2.0;
#else
#endif*/
/// TI moved slope definition to main.c

/// SF commented   double tau_f;

double surface=0.01; //mm²
double Csc; // =0.0185 µF/mm²

///double pKNa = 0.03; //dimensionless
//Constant value to be computed in main
double RTONF; // R*T/F
///double rec_iK1;
double KoSat; //=Ko/(Ko+KmK), dimensionless

//Parameters for Iks
#if defined EPI || defined ENDO
///double Gks_epi = 0.02*2*0.9*4.2*0.4*0.75*sl[slope][GKS]*(1/CAPACITANCE);
double Gks_epi = 0.27;
/*=0.3923027 if SL11 1/ms same in 2006*/
///double Gks_endo=0.02*2*0.9*4.2*0.4*0.75*sl[slope][GKS]*(1/CAPACITANCE);
double Gks_endo=0.27;
/*=0.3923027 if SL11 1/ms same in 2006*/
///double Gks_m=0.02*2*0.9*4.2*0.4*0.75*0.25*(1/CAPACITANCE)*sl[slope][GKS];
double Gks_m=0.27;
/*=0.09807567 if SL11 1/ms same in 2006*/
#elif defined MCELL || defined ENHANCEDMCELL
/// ///double Gks = 0.02*2*0.9*4.2*0.4*0.75*1.6*0.25*(1/CAPACITANCE); /*=0.09807567 1/ms same in 2006*/
#endif

//Parameters for Ito
#if defined EPI || defined MCELL || defined ENHANCEDMCELL
double Gto_epi; /*=0.29372972 1/ms, same in 2006*/
/// SF commented double Gto_m; /*=0.29372972 1/ms, same in 2006*/
/// SF commented double Gto_endo; /*=0.073432432 1/ms, same in 2006*/
#elif defined ENDO
/// ///double Gto; /*=0.073432432 1/ms, same in 2006*/
#endif

double ipmax = 2.7243243;//0.20*0.72*1.75*2*(1/CAPACITANCE);  /*=2.7243243 mV/ms, same pNaK in 2006 */ PK;
//Parameters for INaCa
#if defined EPI || defined ENDO || defined MCELL
double knaca; /*=1000 mV/ms, same knaca in 2006*/
#elif defined ENHANCEDMCELL
/// ///double knaca;
#endif

//To scale the impact of currents on different volumes
double preplength = 0.074;  /*mm*/
double radius = 0.012;  /*mm */
double vc; /*=0.000016403 mm³ */
double vss;  /*=0.000000054678 mm³ */
double vsr; /*=0.0000010935 mm³ */

//parameters for open and close
double k1b=0.15; // 1/(mM²*ms)
double k2b=0.045; // 1/(mM*ms)

/*-------------------ACTION POTENTIAL PROPAGATION----------------------------*/

double S=200; //(1/mm)
double rho=1.75; // (ms/µF)*mm, not the same in 2006

/*-----------------------------------------------------------------------------
  NUMERICAL PARAMETERS GOVERNING ACTION POTENTIAL GENERATION AND PROPAGATION
  SPEED AND PRECISION
  ---------------------------------------------------------------------------*/

/*---------------------PARAMETERS FOR TABULATION-----------------------------*/
double VMIN = -110.;
double VMAX = 80.;
double CAIMIN = 0.;
double CAIMAX = 0.012;

int POINTS=10000;
/*-----------------------------------------------------------------------------
                PARAMETERS FOR INITIAL CONDITIONS
  ---------------------------------------------------------------------------*/
//homogenous initialisation
double Cass_init=0.00007; //mM
//double CaT_init=2.8; //mM PK
double CaT_init=4.272; //mM
double Cai_init=0.00007;
double Nai_init=6; ///1.01338761744302950e+01; /// Tim
double Ki_init=140; ///1.35368636336204105e+02; /// Tim


/*-----------------------------------------------------------------------------
  ELECTROPHYSIOLOGICAL PARAMETERS:
-----------------------------------------------------------------------------*/

//External concentrations
double Ko=5.4; //mM
double Cao=2.0; //mM
double Nao=140.0; //mM

//Intracellular volumes
double Vc=0.000016403;  //mm³
double Vsr=0.0000010935; //mm³
double Vss=0.000000054678; //mm³

//Calcium buffering dynamics
double Bufc=0.2; //mM
double Kbufc=0.001; //mM
double Bufsr=10.; //mM
double Kbufsr=0.3; //mM
double Bufss=0.4; //mM
double Kbufss=0.00025; //mM

//Intracellular calcium flux dynamics
double Vmaxup=0.006375; //mM/ms
double Kup=0.00025; //mM
double Vrel=0.102;//40.8; //1/ms
/// SF commented  double k1_=0.15; //1/(mM²*ms)
/// SF commented  double k2_=0.045; //1/(mM*ms)
double k3=0.060; //1/ms
double k4=0.005;//0.000015; 1/ms
double EC=1.5; //mM
double maxsr=2.5; //dimensionless
double minsr=1.; //dimensionless
double Vleak=0.00036; //1/ms
double Vxfer=0.0038; //1/ms



//Constants
double R=8314.472; //mJ/(K*mol)
double F=96485.3415; //C/mol
double T=310.0; //K

//Cellular capacitance
double Cm=0.000185; //µF
double CAPACITANCE=0.185;

//Parameters for currents
//Parameters for IKr

double Gkr;
/*double Gkr= //1/ms
#ifdef slope07
0.134;
#elif defined(slope11)
0.153;
#elif defined(slope14)
0.172;
#elif defined(slope18)
0.172;
#else
#endif*/
/// TI moved slope definition to main.c


//Parameters for Iks
double pKNa=0.03; //dimensionless

double Gks;
/*double Gks= //1/ms
#ifdef slope07
0.270;
#elif defined(slope11)
0.392;
#elif defined(slope14)
0.441;
#elif defined(slope18)
0.441;
#else
#endif*/
/// TI moved slope definition to main.c


//Parameters for Ik1
double GK1=5.405405; //1/ms
//Parameters for Ito

double Gto=0.29372972; //1/ms

//Parameters for INa
double GNa=14.83783; ///1/ms   coefficient c_gNa added to CommonInit!!!
//Parameters for IbNa
double GbNa=0.00029; //1/ms
//Parameters for INaK
double KmK=1.0; //mM
double KmNa=40.0; //mM
double knak=2.724; //mV/ms
//Parameters for ICaL
double GCaL=0.00003980; //mm³/(ms*µF)
//Parameters for IbCa
double GbCa=0.000592; //1/ms
//Parameters for INaCa
double knaca=1000; //mV/ms
double KmNai=87.5; //mM
double KmCa=1.38; //mM
double ksat=0.1; //mM
double n=0.35; //dimensionless
//Parameters for IpCa

double GpCa;
/*double GpCa= //mV/ms
#ifdef slope07
0.0619;
#elif defined(slope11)
0.1238;
#elif defined(slope14)
0.3714;
#elif defined(slope18)
0.8666;
#else
#endif*/

double KpCa=0.0005; //mM
//Parameters for IpK;

double GpK;
// double GpK= //1/ms
// #ifdef slope07
// 0.0730;
// #elif defined(slope11)
// 0.0146;
// #elif defined(slope14)
// 0.0073;
// #elif defined(slope18)
// 0.00219;
// #else
// #endif


/*-----------------------------------------------------------------------------
                PARAMETERS FOR INITIAL CONDITIONS
------------------------------------------------------------------------------*/
//Initial values of state variables
const double V_init =-8.53798558432631864e+01; /// Tim
const double volt_init=-8.53798558432631864e+01;

//these values are the same at all timesteps and for all grid points
//so their value is static and they can be shared in parallel applications

double Kupsquare;

double VTabConst1;
double VTabConst2;
double CaiTabConst1;
double CaiTabConst2;

typedef struct
{
    //  parameters for tabulation
    int POINTS;
    double VMIN, VMAX, StepV, CAIMIN, CAIMAX, StepCAI;
    //  arrays with the results
    //asymptotic values for time dependent gates: f(V)
    double *M_INF, *H_INF, *J_INF, *Xr1_INF, *Xr2_INF, *Xs_INF, *R_INF,
    *S_INF, *D_INF, *F_INF, *F2_INF;
    //exponents of timeconstants of gates f(V):
    double *EXPTAU_M, *EXPTAU_H, *EXPTAU_J, *EXPTAU_Xr1, *EXPTAU_Xr2,
    *EXPTAU_Xs, /** SF commented  *EXPTAU_Sepi,*EXPTAU_Sendo,*EXPTAU_Sm,  SF added */ *EXPTAU_S,
    *EXPTAU_R, *EXPTAU_D, *EXPTAU_F, *EXPTAU_F2;
    //other voltage dependent terms:
    double /** SF commented *rec_iK1, */ *ghk_iCaL1_abs, *ghk_iCaL2_abs, *exp_iNaCa1_abs, *exp_iNaCa2_abs, /// SF added _abs
    /** SF changed *rec_ip */ *rec_ip_abs; /// , *rec_ipK;   SF commented
    //other calcium dependent terms:
    /// double *ECA;  SF commented
} TTabulation;

TTabulation MyTabulation;

#define OURINA
//#define LRINA

#define OURICAL
//#define LRICAL
//#define PBICAL

TTabulation CalcTabulation(int aPOINTS, double aVMIN, double aVMAX, double aCAIMIN, double aCAIMAX, double HT)
{
    TTabulation res;
    double vv;
    int i;
    /// SF commented   double Ek = RTONF*(log((Ko/Ki_init)));
    /// SF commented   double halfRTONF = RTONF/2;

    res.POINTS = aPOINTS;
    res.VMIN = aVMIN;
    res.VMAX = aVMAX;
    res.StepV = (res.VMAX-res.VMIN)/res.POINTS;
    res.CAIMIN = aCAIMIN;
    res.CAIMAX = aCAIMAX;
    res.StepCAI = (res.CAIMAX-res.CAIMIN)/res.POINTS;


    //  Memory allocations for tables
//    try {
    res.M_INF = calloc(res.POINTS, sizeof(double));
    res.H_INF  = calloc(res.POINTS, sizeof(double));
    res.J_INF  = calloc(res.POINTS, sizeof(double));
    res.Xr1_INF  = calloc(res.POINTS, sizeof(double));
    res.Xr2_INF  = calloc(res.POINTS, sizeof(double));
    res.Xs_INF  = calloc(res.POINTS, sizeof(double));
    res.R_INF  = calloc(res.POINTS, sizeof(double));
    res.S_INF = calloc(res.POINTS, sizeof(double)); /// SF added
    /// SF commented
/// ///    res.S_INFepi = calloc(res.POINTS, sizeof(double));
/// ///    res.S_INFendo = calloc(res.POINTS, sizeof(double));
/// ///    res.S_INFm = calloc(res.POINTS, sizeof(double));
    res.D_INF  = calloc(res.POINTS, sizeof(double));
    res.F_INF  = calloc(res.POINTS, sizeof(double));
    res.F2_INF  = calloc(res.POINTS, sizeof(double));
    res.EXPTAU_M  = calloc(res.POINTS, sizeof(double));
    res.EXPTAU_H  = calloc(res.POINTS, sizeof(double));
    res.EXPTAU_J  = calloc(res.POINTS, sizeof(double));
    res.EXPTAU_Xr1  = calloc(res.POINTS, sizeof(double));
    res.EXPTAU_Xr2  = calloc(res.POINTS, sizeof(double));
    res.EXPTAU_Xs  = calloc(res.POINTS, sizeof(double));
    /// SF commented
///    res.EXPTAU_Sepi  = calloc(res.POINTS, sizeof(double));
///    res.EXPTAU_Sendo  = calloc(res.POINTS, sizeof(double));
///    res.EXPTAU_Sm  = calloc(res.POINTS, sizeof(double));
    res.EXPTAU_S  = calloc(res.POINTS, sizeof(double));
    res.EXPTAU_R  = calloc(res.POINTS, sizeof(double));
    res.EXPTAU_D  = calloc(res.POINTS, sizeof(double));
    res.EXPTAU_F  = calloc(res.POINTS, sizeof(double));
    res.EXPTAU_F2 = calloc(res.POINTS, sizeof(double));
    ///res.rec_iK1 = calloc(res.POINTS, sizeof(double));  SF commented
    res.ghk_iCaL1_abs = calloc(res.POINTS, sizeof(double));
    res.ghk_iCaL2_abs = calloc(res.POINTS, sizeof(double));
    res.exp_iNaCa1_abs  = calloc(res.POINTS, sizeof(double)); /// SF added _abs
    res.exp_iNaCa2_abs  = calloc(res.POINTS, sizeof(double)); /// SF added _abs
    /** SF changed res.rec_ip = */ res.rec_ip_abs = calloc(res.POINTS, sizeof(double));
    /// res.rec_ipK = calloc(res.POINTS, sizeof(double)); SF commented
    /// res.ECA = calloc(res.POINTS, sizeof(double));  SF commented

//    } catch (bad_alloc) {
//        if(!quiet) cerr << "error allocating memory" << endl;
//        exit(1);
//    }

    //Fill the voltage tabulation tables:
    for (i=0; i<res.POINTS; i++)
    {
        vv = res.VMIN + i*res.StepV;

#if defined OURINA
        double AM = 1./(1.+exp((-60.-vv)/5.));
        double BM = 0.1/(1.+exp((vv+35.)/5.)) + 0.10/(1.+exp((vv-50.)/200.));
        res.EXPTAU_M[i] = AM*BM;
//        res.EXPTAU_M[i] = res.EXPTAU_M[i]; // /1000 removed
        res.M_INF[i] =1./sqr(1.+exp((-56.86-vv)/9.03));
        double AH_1 = 0.;
        double BH_1 = 0.77/(0.13*(1.+exp(-(vv+10.66)/11.1)));
        double AH_2 = 0.057*exp(-(vv+80.)/6.8);
        double BH_2 = 2.7*exp(0.079*vv) + 3.1e5*exp(0.3485*vv);
        if (vv >= -40.)
            res.EXPTAU_H[i] = 1.0/(AH_1+BH_1);
        else
            res.EXPTAU_H[i] = 1.0/(AH_2+BH_2);
//        res.EXPTAU_H[i] = res.EXPTAU_H[i];  // /1000 removed
    #if defined EPI
            res.H_INF[i] = 1./sqr(1.+exp((vv+71.55)/7.43));
    /// SF commented
/// ///    #elif defined ENDO
/// ///            res.H_INF[i] = 1./((1.+exp((vv+71.55)/7.43))*(1.+exp((vv+71.55)/7.43)));
/// ///    #elif defined MCELL
/// ///            res.H_INF[i] = 1./((1.+exp((vv+71.55)/7.43))*(1.+exp((vv+71.55)/7.43)));
/// ///    #elif defined ENHANCEDMCELL
/// ///            res.H_INF[i] = 0.985/((1.+exp((vv+71.55)/7.43))*(1.+exp((vv+71.55)/7.43)))
/// ///                           + 0.015;
    #endif
        double AJ_1 = 0.;
        double BJ_1 = (0.6*exp((0.057)*vv)/(1.+exp(-0.1*(vv+32.))));
        double AJ_2 = (-2.5428e4*exp(0.2444*vv) - 6.948e-6*exp(-0.04391*vv)) *
                      (vv+37.78)/(1.+exp(0.311*(vv+79.23)));
        double BJ_2 = 0.02424*exp(-0.01052*vv) / (1.+exp(-0.1378*(vv+40.14)));
        if (vv >= -40.)
            res.EXPTAU_J[i] = 1.0/(AJ_1+BJ_1);
        else
            res.EXPTAU_J[i] = 1.0/(AJ_2+BJ_2);
        /* TODO: Check what's going on here */
//        res.EXPTAU_J[i]=res.EXPTAU_J[i];  // /1000 removed
        //EXPTAU_J[i]=EXPTAU_J[i]*0.1;
        //EXPTAU_J[i]=EXPTAU_J[i]*0.5;

        res.J_INF[i]=res.H_INF[i];
        /// SF commented
/// ///#elif defined  LRINA

/// ///#define AH_1(V) 0.
/// ///#define AJ_1(V) 0.
/// ///#define BH_1(V) (1./(0.13*(1.+exp(-(V+10.66)/11.1))))
/// ///#define BJ_1(V) (0.3*exp((-2.535e-7)*V)/(1.+exp(-0.1*(V+32.))))
/// ///#define AH_2(V) (0.135*exp(-(V+80.)/6.8))
/// ///#define AJ_2(V) (((-1.2714e5)*exp(0.2444*V)-(3.474e-5)*exp(-0.04391*V))*(V+37.78)/(1.+exp(0.311*(V+79.23))))
/// ///#define BH_2(V) (3.56*exp(0.079*V)+(3.1e5)*exp(0.35*V))
/// ///#define BJ_2(V) (0.1212*exp(-0.01052*V)/(1.+exp(-0.1378*(V+40.14))))
/// ///#define AM(V) (0.32*(V+47.13)/(1.-exp(-0.1*(V+47.13))))
/// ///#define BM(V) (0.08*exp(-V/11.))

/// ///        if (vv<-47.129 && vv>-47.131)    //vv 47.13 gives division by zero in AM
/// ///            EXPTAU_M[i] = 1.0/(3.2 + BM(vv));//limit value of AM equals 3.2
/// ///        else
/// ///            EXPTAU_M[i] = 1.0/(AM(vv)+BM(vv));
/// ///        EXPTAU_M[i] = EXPTAU_M[i]/1000.;
/// ///        //asymptotic value:
/// ///        if (vv<-47.129 && vv>-47.131)   //vv 47.13 gives division by zero in AM
/// ///            M_INF[i] = 3.2/(3.2 + BM(vv));//limit value of AM equals 3.2
/// ///        else
/// ///            M_INF[i] = AM(vv)/(AM(vv )+BM(vv ));
/// ///        //H gate:
/// ///        //time constant:
/// ///        if (vv>=-40.)
/// ///            EXPTAU_H[i] = 1.0/(AH_1(vv )+BH_1(vv ));
/// ///        else
/// ///            EXPTAU_H[i] = 1.0/(AH_2(vv )+BH_2(vv ));
/// ///        EXPTAU_H[i] =EXPTAU_H[i]/1000.;
/// ///        //asymptotic value:
/// ///        if (vv>=-40.)
/// ///            H_INF[i] = AH_1(vv )/(AH_1(vv )+BH_1(vv ));
/// ///        else
/// ///            H_INF[i] = AH_2(vv)/(AH_2(vv )+BH_2(vv ));
/// ///        //J gate
/// ///        //time constant:
/// ///        if (vv >= -40.)
/// ///            EXPTAU_J[i] = 1.0/(AJ_1(vv )+BJ_1(vv ));
/// ///        else
/// ///            EXPTAU_J[i] = 1.0/(AJ_2(vv )+BJ_2(vv ));
/// ///        EXPTAU_J[i] =EXPTAU_J[i]/1000.;
/// ///        //asymptotic value
/// ///        if (vv >= -40.)
/// ///            J_INF[i] = AJ_1(vv )/(AJ_1(vv )+BJ_1(vv ));
/// ///        else
/// ///            J_INF[i] = AJ_2(vv )/(AJ_2(vv )+BJ_2(vv ));
#endif

        //IKr
        res.Xr1_INF[i]=1./(1.+exp((-26.-vv)/7.));
        double axr1=450./(1.+exp((-45.-vv)/10.));
        double bxr1=6./(1.+exp((vv-(-30.))/11.5));
        res.EXPTAU_Xr1[i]=axr1*bxr1;//1122*exp(-(vv+35)*(vv+35)/400)+100;
//        res.EXPTAU_Xr1[i]=res.EXPTAU_Xr1[i];  // /1000 removed
        res.Xr2_INF[i]=1./(1.+exp((vv-(-88.))/24.));
        double axr2=3./(1.+exp((-60.-vv)/20.));
        double bxr2=1.12/(1.+exp((vv-60.)/20.));
        res.EXPTAU_Xr2[i]=axr2*bxr2;
//        res.EXPTAU_Xr2[i]=res.EXPTAU_Xr2[i];  // /1000 removed


        //IKs
        //Xs_INF[i]=1./(1.+exp((-45.-vv)/14.));
        res.Xs_INF[i]=1./(1.+exp((-5.-vv)/14.));
        //Normal fit to data:
        //TODO: Ask Sasha about Axs and Bxs
        //double Axs=1100./(sqrt(1.+exp((-10.-vv)/6)));
        //double Bxs=1./(1.+exp((vv-60.)/20.));
        //EXPTAU_Xs[i]=Axs*Bxs;
        res.EXPTAU_Xs[i] = 1400./sqrt(1.+exp((5.-vv)/6)) /
                           (1.+exp((vv-35.)/15.)) + 80;
//        res.EXPTAU_Xs[i] = res.EXPTAU_Xs[i]; // /1000 removed
        //if(vv>0)
        //EXPTAU_Xs[i]=EXPTAU_Xs[i]*0.7;
        //if(vv<-40)
        //EXPTAU_Xs[i]=EXPTAU_Xs[i]*0.5;

        //IK1
        //voltage dependent exponent
        /// SF commented
        ///double Ak1 = 0.1/(1.+exp(0.06*(vv-Ek-200)));
        ///double Bk1 = (3.*exp(0.0002*(vv-Ek+100)) + exp(0.1*(vv-Ek-10))) / (1.+exp(-0.5*(vv-Ek)));
        ///res.rec_iK1[i] = Ak1/(Ak1+Bk1);

        //Ip
///        res.rec_ip[i] = 1./(1.+0.1245*exp(-0.1*vv/RTONF)+0.0353*exp(-vv/RTONF));   SF commented
        res.rec_ip_abs[i] = (ipmax*(Ko/(Ko+KmK))*(Nai_init/(Nai_init+KmNa)))/(1.+0.1245*exp(-0.1*vv/RTONF)+0.0353*exp(-vv/RTONF));  /// SF added

        //IpK
        ///res.rec_ipK[i] = GpK/(1.+exp((25-vv)/5.98)); SF commented

        //Ito
#ifdef EPI
        res.R_INF[i] = 1./(1.+exp(((30-10)-vv)/6.));
        res.S_INF[i] = 1./(1.+exp((vv-(-10-10))/5.)); /// SF added
        /// SF commented
/// ///        res.S_INFepi[i] = 1./(1.+exp((vv-(-10-10))/5.));
/// ///        res.S_INFendo[i]=1./(1.+exp((vv-(-18-10))/5.));
/// ///        res.S_INFm[i]=1./(1.+exp((vv-(-10-10))/5.));
        res.EXPTAU_R[i] = 9.5*exp(-(vv+40.)*(vv+40.)/1800.)+0.8;
//        res.EXPTAU_R[i] = EXPTAU_R[i];  // /1000 removed
///        double EXPTAU_S_ = 85.*exp(-(vv+45.)*(vv+45.)/320.)+ 5./(1.+exp((vv-20.)/5.))+3.; SF changed
        res.EXPTAU_S[i] = 85.*exp(-(vv+45.)*(vv+45.)/320.)+ 5./(1.+exp((vv-20.)/5.))+3.;
        /// SF commented
///        res.EXPTAU_Sepi[i] = EXPTAU_S_; /// SF changed 85.*exp(-(vv+45.)*(vv+45.)/320.)+ 5./(1.+exp((vv-20.)/5.))+3.;
///        res.EXPTAU_Sendo[i] = EXPTAU_S_; /// SF changed 1000.*exp(-(vv+67)*(vv+67)/1000.)+8.;
///        res.EXPTAU_Sm[i] = EXPTAU_S_; /// SF changed 85.*exp(-(vv+45.)*(vv+45.)/320.) + 5./(1.+exp((vv-20.)/5.))+3.;
        //   EXPTAU_S[i] = EXPTAU_S[i];   // /1000 removed
/// ///#elif defined ENDO
/// ///        res.R_INF[i] = 1./(1.+exp(((30-10)-vv)/6.));
/// ///        res.S_INF[i] = 1./(1.+exp((vv-(-18-10))/5.));
/// ///        res.EXPTAU_R[i] = 9.5*exp(-(vv+40.)*(vv+40.)/1800.)+0.8;
/// /////        res.EXPTAU_R[i] = res.EXPTAU_R[i];  // /1000 removed
/// ///        res.EXPTAU_S[i] = 1000.*exp(-(vv+67)*(vv+67)/1000.)+8.;
/// /////        res.EXPTAU_S[i] = res.EXPTAU_S[i];  // /1000 removed
/// ///#elif defined MCELL
/// ///        R_INF[i] = 1./(1.+exp(((30-10)-vv)/6.));
/// ///        S_INF[i] = 1./(1.+exp((vv-(-10-10))/5.));
/// ///        EXPTAU_R[i] = 9.5*exp(-(vv+40.)*(vv+40.)/1800.)+0.8;
/// /////        EXPTAU_R[i] = EXPTAU_R[i];  // /1000 removed
/// ///        EXPTAU_S[i] = 85.*exp(-(vv+45.)*(vv+45.)/320.) +
/// ///                      5./(1.+exp((vv-20.)/5.))+3.;
/// /////        EXPTAU_S[i]=EXPTAU_S[i];  // /1000 removed
/// ///#elif defined ENHANCEDMCELL
/// ///        R_INF[i] = 1./(1.+exp(((30-10)-vv)/6.));
/// ///        S_INF[i] = 1./(1.+exp((vv-(-10-10))/5.));
/// ///        EXPTAU_R[i] = 9.5*exp(-(vv+40.)*(vv+40.)/1800.) + 0.8;
/// /////        EXPTAU_R[i] = EXPTAU_R[i];  // /1000 removed
/// ///        EXPTAU_S[i] = 85.*exp(-(vv+45.)*(vv+45.)/320.) +
/// ///                      5./(1.+exp((vv-20.)/5.))+3.;
/// /////        EXPTAU_S[i] = EXPTAU_S[i]; // /1000 removed
#endif

#ifdef OURICAL
        //ICa
        res.D_INF[i] = 1./(1.+exp((-8-vv)/7.5));
        double Ad = 1.4/(1.+exp((-35-vv)/13))+0.25;
        double Bd = 1.4/(1.+exp((vv+5)/5));
        double Cd = 1./(1.+exp((50-vv)/20));
        res.EXPTAU_D[i] = Ad*Bd+Cd;
        /// ///res.EXPTAU_D[i] = res.EXPTAU_D[i]; // /1000 removed
        res.F_INF[i] = 1./(1.+exp((vv+20)/7));
        double Af = 2.45*450*exp(-(vv+27)*(vv+27)/225);
        double Bf = 200/(1.+exp((13-vv)/10));
        double Cf = 180/(1.+exp((vv+30)/10));
        res.EXPTAU_F[i] = Af+Bf+Cf+20;
        /// ///res.EXPTAU_F[i] = res.EXPTAU_F[i];  // /1000 removed
        /** if (vv>0) SF commented */ res.EXPTAU_F[i] = res.EXPTAU_F[i]*TAU_F_INACT_FACTOR;//*0.6;//*3;
        res.F2_INF[i] = 0.67/(1.+exp((vv+35)/7))+0.33;
        //F2_INF[i] = 0.77/(1.+exp((vv+35)/7))+0.23;
        double Af2 = 562*exp(-(vv+27)*(vv+27)/240);//PK
        double Bf2 = 31/(1.+exp((25-vv)/10));
        double Cf2 = 80/(1.+exp((vv+30)/10));//PK
        res.EXPTAU_F2[i] = Af2+Bf2+Cf2;
        if (vv<-40) res.EXPTAU_F2[i] = 0.2*res.EXPTAU_F2[i];//PK
        /// ///res.EXPTAU_F2[i] = res.EXPTAU_F2[i];  // /1000 removed
        // default:   if (vv<-40)
        //   EXPTAU_F2[i] = EXPTAU_F2[i]*0.2;//*0.3;
///        res.ghk_iCaL1[i] = 4*(vv-15.)*(F/RTONF)*(0.25*exp(2*(vv-15.)/RTONF)) / (exp(2*(vv-15.)/RTONF)-1.); SF commented
///        res.ghk_iCaL2[i] = 4*(vv-15.)*(F/RTONF)*(1*Cao)/(exp(2*(vv-15.)/RTONF)-1.);                        SF commented
        res.ghk_iCaL1_abs[i] = GCaL*4*(vv-15.)*(F/RTONF)*(0.25*exp(2*(vv-15.)/RTONF)) / (exp(2*(vv-15.)/RTONF)-1.);
        res.ghk_iCaL2_abs[i] = GCaL*4*(vv-15.)*(F/RTONF)*(1*Cao)/(exp(2*(vv-15.)/RTONF)-1.);
/// ///#elif defined LRICAL
/// ///        double ad=(0.095*exp(-0.01*(vv-5.))/(1.+exp(-0.072*(vv-5.))));
/// ///        double bd=(0.07*exp(-0.017*(vv+44.))/(1.+exp(0.05*(vv+44.))));
/// ///        double af=(0.012*exp(-0.008*(vv+28.))/(1.+exp(0.15*(vv+28.))));
/// ///        double bf=(0.0065*exp(-0.02*(vv+30.))/(1.+exp(-0.2*(vv+30.))));
/// ///        D_INF[i]=ad/(ad+bd);
/// ///        EXPTAU_D[i]=1./(ad+bd);
/// ///        EXPTAU_D[i]=EXPTAU_D[i]/1000.;
/// ///        F_INF[i]=af/(af+bf);
/// ///        EXPTAU_F[i]=1./(af+bf);
/// ///        EXPTAU_F[i]=EXPTAU_F[i]/1000.;
/// ///#elif defined PBICAL
/// ///        double ad = 14.9859/(16.6813*sqrt(2.0*M_PI)) *
/// ///                    exp((-0.5*(vv-22.36)*(vv-22.36))/(16.6813*16.6813));
/// ///        double bd = 0.1471 - 5.3/(14.93*sqrt(2.0*M_PI)) *
/// ///                    exp((-0.5*(vv-6.2744)*(vv-6.2744))/(14.93*14.93));
/// ///        double af = 0.006872/(1.0+exp((vv-6.1546)/6.1223));
/// ///        double bf = 0.0005474 + (0.0687*exp(-0.1081*(vv+9.8255)) + 0.0112) /
/// ///                    (1.0+exp(-0.2779*(vv+9.8255)));
/// ///        D_INF[i] = ad/(ad+bd);
/// ///        EXPTAU_D[i] = 1./(ad+bd);
/// ///        EXPTAU_D[i] = EXPTAU_D[i]/1000.;
/// ///        F_INF[i] = af/(af+bf);
/// ///        EXPTAU_F[i] = 1./(af+bf);
/// ///        EXPTAU_F[i] = EXPTAU_F[i]/1000.;
#endif

        //INaCa
        //voltage dependent exponents
        res.exp_iNaCa1_abs[i] = knaca /(KmNai*KmNai*KmNai+Nao*Nao*Nao)
                            /(KmCa+Cao) /(1+ksat*exp((n-1)*vv/RTONF)) * exp(n*vv/RTONF) * Cao * kub(Nai_init);
        res.exp_iNaCa2_abs[i] = knaca /(KmNai*KmNai*KmNai+Nao*Nao*Nao) /(KmCa+Cao)
                            /(1+ksat*exp((n-1)*vv/RTONF)) * exp((n-1)*vv/RTONF) *
                            Nao*Nao*Nao * 2.5;

        //Instead of the timeconstant the exp(-HT/timeconstants) is tabulated
        //to even further limit computations during simulation.
        res.EXPTAU_M[i] = exp(-HT/res.EXPTAU_M[i]);
        res.EXPTAU_H[i] = exp(-HT/res.EXPTAU_H[i]);
        res.EXPTAU_J[i] = exp(-HT/res.EXPTAU_J[i]);
        res.EXPTAU_Xr1[i] = exp(-HT/res.EXPTAU_Xr1[i]);
        res.EXPTAU_Xr2[i] = exp(-HT/res.EXPTAU_Xr2[i]);
        res.EXPTAU_Xs[i] = exp(-HT/res.EXPTAU_Xs[i]);
        /// SF commented
///        res.EXPTAU_Sepi[i] = exp(-HT/res.EXPTAU_Sepi[i]);
///        res.EXPTAU_Sendo[i] = exp(-HT/res.EXPTAU_Sendo[i]);
///        res.EXPTAU_Sm[i] = exp(-HT/res.EXPTAU_Sm[i]);
        /// SF added
        res.EXPTAU_S[i] = exp(-HT/res.EXPTAU_S[i]);
        res.EXPTAU_R[i] = exp(-HT/res.EXPTAU_R[i]);
        res.EXPTAU_D[i] = exp(-HT/res.EXPTAU_D[i]);
        res.EXPTAU_F[i] = exp(-HT/res.EXPTAU_F[i]);
        res.EXPTAU_F2[i] = exp(-HT/res.EXPTAU_F2[i]);


        ///inclusion of maximum conductance values in activation steady
        ///state: so maximum of steady state activation is not 1 but Gx
/// SF commented
/// ///        res.M_INF[i] = res.M_INF[i];
/// ///        res.Xr1_INF[i] = res.Xr1_INF[i];
/// ///        res.Xs_INF[i] = res.Xs_INF[i];
/// ///        res.R_INF[i] = res.R_INF[i];
/// ///        res.D_INF[i] = res.D_INF[i];
/// ///        res.rec_iK1[i] = res.rec_iK1[i];
        //   S_INF[i] = S_INF[i];
    }
/**  SF commented
    double Cai;
    for (i=0; i<res.POINTS; i++)
    {
        Cai = res.CAIMIN + i*res.StepCAI;
        res.ECA[i] = halfRTONF*(log(Cao/Cai));
    }
*/
    return res;
}

void CommonInit(double HT)
{
    RTONF = R*T/F;

    //these values are the same at all timesteps and for all grid points
    //so their value is static and they can be shared in parallel applications

    Kupsquare = Kup*Kup;

    // from Parameters.h
    /// SF commented   tau_f = TAU_F_INACT_FACTOR;
    Csc=(Cm/surface); // =0.0185 µF/mm²
    vc=0.49*3.14159265*radius*radius*preplength; /*=0.000016403 mm³ */
    vss=vc/300;  /*=0.000000054678 mm³ */
    vsr=vc/15; /*=0.0000010935 mm³ */

    //Parameters for Ito
#if defined EPI || defined MCELL || defined ENHANCEDMCELL
    Gto_epi= 0.05434*(1/CAPACITANCE); /*=0.29372972 1/ms, same in 2006*/
    ///SF commented  Gto_m=0.05434*(1/CAPACITANCE); /*=0.29372972 1/ms, same in 2006*/
    ///SF commented  Gto_endo=0.0247*1.1*2*0.25*(1/CAPACITANCE); /*=0.073432432 1/ms, same in 2006*/
/// ///#elif defined ENDO
/// ///    Gto = 0.0247*1.1*2*0.25*(1/CAPACITANCE); /*=0.073432432 1/ms, same in 2006*/
#endif

    //Parameters for INaCa
#if defined EPI || defined ENDO || defined MCELL
    knaca = 1000*0.185*(1/CAPACITANCE); /*=1000 mV/ms, same knaca in 2006*/
/// ///#elif defined ENHANCEDMCELL
/// ///    knaca = 1000*0.185*0.6;
#endif

    if (amiodaron==1)
    {
        Gkr*=0.7128;
        GNa*=0.7633;
        knak*=0.9398;
        GCaL*=0.8529;
        knaca*=0.7674;
        Gks*=0.6942;
    }
    if (amiodaron==3)
    {
        Gkr*=0.4841;
        GNa*=0.5878;
        knak*=0.8387;
        GCaL*=0.6591;
        knaca*=0.5238;
        Gks*=0.5373;
    }

    GNa*=c_gNa; ///1/ms   coefficient c_gNa added!!!
    GCaL*=c_gCaL; ///1/ms   coefficient c_gCaL added!!!
    GCaL*=c_verapamil; Gkr*=c_verapamil;

    MyTabulation = CalcTabulation(POINTS,VMIN,VMAX,CAIMIN,CAIMAX,HT);
/**         FOR DEBUG
    int i;
    printf ("print tabulation data...");
    FILE *ftab=fopen("tab.dat","w");
    //Выводим заголовки столбцов
    fprintf(ftab, "Volt M_INF H_INF J_INF Xr1_INF Xr2_INF Xs_INF R_INF S_INF D_INF F_INF F2_INF ");
    fprintf(ftab, "EXPTAU_M EXPTAU_H EXPTAU_J EXPTAU_Xr1 EXPTAU_Xr2 EXPTAU_Xs EXPTAU_S EXPTAU_R EXPTAU_D EXPTAU_F EXPTAU_F2 ");
    fprintf(ftab, "ghk_iCaL1_abs ghk_iCaL2_abs exp_iNaCa1_abs exp_iNaCa2_abs rec_ip_abs\n"); ///" rec_ipK\n");
    /// fprintf(ftab, "Ca ECA\n");

    #define t MyTabulation

    for (i=0; i<POINTS; i++)
    {
        fprintf(ftab, "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f ",
                t.VMIN + i*t.StepV,
                t.M_INF[i], t.H_INF[i], t.J_INF[i], t.Xr1_INF[i], t.Xr2_INF[i], t.Xs_INF[i],
                t.R_INF[i], t.S_INF[i], t.D_INF[i], t.F_INF[i], t.F2_INF[i]);
        fprintf(ftab, "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f ",
                t.EXPTAU_M[i], t.EXPTAU_H[i], t.EXPTAU_J[i], t.EXPTAU_Xr1[i], t.EXPTAU_Xr2[i],
                t.EXPTAU_Xs[i], t.EXPTAU_S[i], t.EXPTAU_R[i],
                t.EXPTAU_D[i], t.EXPTAU_F[i], t.EXPTAU_F2[i]);
        fprintf(ftab, "%.8f %.8f %.8f %.8f %.8f\n",
                /// t.rec_iK1[i],
                t.ghk_iCaL1_abs[i], t.ghk_iCaL2_abs[i], t.exp_iNaCa1_abs[i], t.exp_iNaCa2_abs[i],
                t.rec_ip_abs[i]); ///, t.rec_ipK[i]);
        ///  fprintf(ftab, "%.8f %.8f\n", t.CAIMIN + i*t.StepCAI, t.ECA[i]);
    }
    fclose (ftab);
    #undef t
    printf ("ok!\n");
*/
    VTabConst1 = -(MyTabulation.VMIN / MyTabulation.StepV) + 0.5;
    VTabConst2 = 1. / MyTabulation.StepV;
    CaiTabConst1 = -(MyTabulation.CAIMIN / MyTabulation.StepCAI) + 0.5;
    CaiTabConst2 = 1. / MyTabulation.StepCAI;
}

double total_of_free(double free, double buf, double k)
{
    return free * (1. + buf / (free + k));
}

double ca_i_total_of_free(double ca)
{
    return total_of_free(ca, Bufc, Kbufc);
}

double ca_sr_total_of_free(double ca)
{
    return total_of_free(ca, Bufsr, Kbufsr);
}

double ca_ss_total_of_free(double ca)
{
    return total_of_free(ca, Bufss, Kbufss);
}
#if 0
void InitializeNewVariables(double* V, TTabulation *Ta)
{
    V/*->pVolt*/[0] = V_init;
    V/*->pCai_total*/[1] = ca_i_total_of_free(Cai_init); /// SF changed Cai_init;
    V/*->pCaSR_total*/[2] = ca_sr_total_of_free(CaT_init); /// SF changed CaT_init;
    V/*->pCaSS_total*/[3] = ca_ss_total_of_free(Cass_init); /// SF changed Cass_init;
    V/*->pM*/[4] = 1.66171456717857527e-03; /// Tim
    V/*->pH*/[5] = 7.48991583078272738e-01; /// Tim
    V/*->pJ*/[6] = 7.48599753512242727e-01; /// Tim
    V/*->pXr1*/[7] = 2.07650883539678194e-04; /// Tim
    V/*->pXr2*/[8] = 4.72733193318403078e-01; /// Tim
    V/*->pXs*/[9] = 3.23090984071628057e-03; /// Tim
    V/*->pR*/[10] = 2.35709886030767176e-08; /// Tim
    V/*->pS*/[11] = 9.99997904693742057e-01; /// Tim
    V/*->pD*/[12] = 3.30548499869733964e-05; /// Tim
    V/*->pF*/[13] = 9.77158843631106722e-01; /// Tim
    V/*->pF2*/[14] = 9.99406290178190937e-01; /// Tim
    V/*->pFCASS*/[15] = 9.99972178434477055e-01; /// Tim
    V/*->pRbar*/[16] = 9.89066126100902498e-01; /// Tim
}
#endif // 0
#if 1
void InitializeNewVariables(double* V, TTabulation *Ta)
{
    V/*->pVolt*/[0] = V_init;
    V/*->pCai_total*/[1] = ca_i_total_of_free(Cai_init); /// SF changed Cai_init;
    V/*->pCaSR_total*/[2] = ca_sr_total_of_free(CaT_init); /// SF changed CaT_init;
    V/*->pCaSS_total*/[3] = ca_ss_total_of_free(Cass_init); /// SF changed Cass_init;
    V/*->pM*/[4] = 0; /// SF changed
    V/*->pH*/[5] = 0.77; /// SF changed
    V/*->pJ*/[6] = 0.77; /// SF changed
    V/*->pXr1*/[7] = 0.1; /// SF changed
    V/*->pXr2*/[8] = 1; /// SF changed
    V/*->pXs*/[9] = 0; /// SF changed
    V/*->pR*/[10] = 0; /// SF changed
    V/*->pS*/[11] = 1; /// SF changed
    V/*->pD*/[12] = 0; /// SF changed
    V/*->pF*/[13] = 0.5;
    V/*->pF2*/[14] = 1; /// SF changed
    V/*->pFCASS*/[15] = 1.;
    V/*->pRbar*/[16] = 1.;
}
#endif

void InitUV(double* uv)
{
    InitializeNewVariables (uv, &MyTabulation);
}

double qsolve(double b, double minus_c)
{
    double sqrt_d;
    sqrt_d = sqrt(sqr(b) + 4 * minus_c);
    if (b > 0)
        return 2. * minus_c / (sqrt_d + b);
    else
        return 0.5 * (sqrt_d - b);
}

double free_of_total(double total, double buf, double k)
{
    double b, minus_c;
    b = buf - total + k;
    minus_c = k * total;
    return qsolve(b, minus_c);
}

double rev_potential_ca(double ca_in)
{
    return RTONF * 0.5 * log(Cao / ca_in);
}

double RecIpK(double v) /// k_plateau_abs
{
    return GpK / (1. + exp((25. - v) / 5.98));
}

void Step(double *V, double *Vn, TTabulation *Ta, double DivDGradU, double HT, double Istim)
{

    double IKr, IKs, IK, IK1, Ito, INa, IbNa, ICaL, IbCa, INaCa, IpCa, IpK, Ip,
    Irel, Ileak, Jxfer;

    int Vtab; /// SF commented  , Caitab;

    double /** dNai, dKi, SF commented */ dCai, dCaT, dCaSS, SERCA, VminEK, VminENa,
    /// CaCSQN, bjsr, cjsr, CaSSBuf, bcss, ccss, CaBuf, bc, cc,     SF commented
    Ak1, Bk1, rec_iK1,
    taufcass, fcassinf, Ek, Ena, Eks; ///, axialcurrent, tauf,finf;     SF commented

    // define all variables
#define       svolt        (V[0] )
#define       sCai_total   (V[1] )
#define       sCaSR_total  (V[2] )
#define       sCaSS_total  (V[3] )
#define       sm           (V[4] )
#define       sh           (V[5] )
#define       sj           (V[6] )
#define       sxr1         (V[7] )
#define       sxr2         (V[8] )
#define       sxs          (V[9] )
#define       ss           (V[11] )
#define       sr           (V[10] )
#define       sd           (V[12] )
#define       sf           (V[13] )
#define       sf2          (V[14] )
#define       sfcass       (V[15] )
#define       sRbar        (V[16] )

#define       svolt_n      (Vn[0] )
#define       sCai_total_n (Vn[1] )
#define       sCaSR_total_n (Vn[2] )
#define       sCaSS_total_n (Vn[3] )
#define       sm_n         (Vn[4] )
#define       sh_n         (Vn[5] )
#define       sj_n         (Vn[6] )
#define       sxr1_n       (Vn[7] )
#define       sxr2_n       (Vn[8] )
#define       sxs_n        (Vn[9] )
#define       ss_n         (Vn[11] )
#define       sr_n         (Vn[10] )
#define       sd_n         (Vn[12] )
#define       sf_n         (Vn[13] )
#define       sf2_n        (Vn[14] )
#define       sfcass_n     (Vn[15])
#define       sRbar_n      (Vn[16] )

double sItot_n;

///#define       RecIK1      Ta->rec_iK1[Vtab]    SF commented
///#define       GHKICaL1    Ta->ghk_iCaL1[Vtab]  SF commented
///#define       GHKICaL2    Ta->ghk_iCaL2[Vtab]  SF commented
#define       GHKICaL1_abs    Ta->ghk_iCaL1_abs[Vtab]
#define       GHKICaL2_abs    Ta->ghk_iCaL2_abs[Vtab]
/// SF changed: twice _abs added
#define       ExpINaCa1_abs   Ta->exp_iNaCa1_abs[Vtab]
#define       ExpINaCa2_abs   Ta->exp_iNaCa2_abs[Vtab]
///#define       RecIp       Ta->rec_ip[Vtab]     SF commented
#define       RecIp_abs       Ta->rec_ip_abs[Vtab]
///#define       RecIpK      Ta->rec_ipK[Vtab]    SF commented

    double dRR;
    double kCaSR;
/** SF commented
    Vtab = (int)(VTabConst2*svolt+VTabConst1);
    if (Vtab > Ta->POINTS-1||Vtab<0)
    {
        if (debugEl)
        {
            printf("ERROR: Vtab outside range of tabulation\n");
            printf("Current value: %f\n",svolt);
            printf("Range: %f to %f\n",(double)Ta->VMIN,(double)Ta->VMAX);
//                    printf("Step is: %i\n",step);
        }
        exit(EXIT_FAILURE);
    }
*/
/** SF commented
    Caitab=(int)(CaiTabConst2*sCai_total+CaiTabConst1);
    if (Caitab > Ta->POINTS-1||Caitab<0)
    {
        if (debugEl)
        {
            printf("ERROR: Caitab outside range of tabulation\n");
            printf("Current value: %f\n",sCai_total);
            printf("Range: %f to %f\n",
                   (double)Ta->CAIMIN, (double)Ta->CAIMAX);
//                    printf("Step is: %i\n",step);
///            printf("i = %d j = %d\n", i, j);
        }
        exit(EXIT_FAILURE);
    }
*/

    if ((svolt<VMIN) || (svolt>VMAX))
    {
        if (debugEl)
        {
            printf("ERROR: svolt outside range of tabulation\n");
            printf("Current value: %f\n",svolt);
            printf("Range: %f to %f\n",(double)VMIN,(double)VMAX);
        }
        exit (EXIT_FAILURE);
    }

    #define TAB_INDEX(x, tmin, tmax) (((x) - tmin) / (tmax - tmin))

    double index=TAB_INDEX(svolt,VMIN,VMAX);
    Vtab=(int)(index*POINTS);
    double ca_total=sCai_total;
    double ca=free_of_total(ca_total,Bufc,Kbufc);

    if ((ca<CAIMIN) || (ca>CAIMAX))
    {
        if (debugEl)
        {
            printf("ERROR: ca outside Cai range\n");
            printf("Current value: %f\n",ca);
            printf("Range: %f to %f\n", (double)CAIMIN, (double)CAIMAX);
        }
        exit (EXIT_FAILURE);
    }
    /// double Gks_a=1.6;
    /// GbNa = 0.000176*0.305*(1/CAPACITANCE);                     /// we use a constant value defined above
    /// Gkr = 0.00945*5*0.5*0.75*1.6*sqrt(Ko/5.4)*(1/CAPACITANCE); /// we use a constant value defined above

    /// Gks_epi =0.02*2*0.9*4.2*0.4*0.75*Gks_a*(1/CAPACITANCE);    /// we use a constant value defined above


    Ek  = RTONF*(log((Ko/Ki_init))); /// SF changed Ek  = RTONF*(log((Ko/Ki)));
    Ena = RTONF*(log((Nao/Nai_init))); /// SF changed Ena = RTONF*(log((Nao/Nai)));
    Eks = RTONF*(log((Ko+pKNa*Nao)/(Ki_init+pKNa*Nai_init))); /// SF changed Eks = RTONF*(log((Ko+pKNa*Nao)/(Ki+pKNa*Nai)));

    VminEK  = svolt-Ek;
    VminENa = svolt-Ena;

    IKr = Gkr*sxr1*sxr2*VminEK;
    IKs=Gks_epi*sxs*sxs*(svolt-Eks);
    IK  = IKr+IKs;

    Ak1 = 0.1/(1.+exp(0.06*(VminEK-200)));
    Bk1 = (3.*exp(0.0002*(VminEK+100)) +
           exp(0.1*(VminEK-10)))/(1.+exp(-0.5*(VminEK)));
    rec_iK1 = Ak1/(Ak1+Bk1);
    IK1 = GK1*rec_iK1*VminEK;
    Ito=Gto_epi*sr*ss*VminEK; /// SF changed Ito=Gto_epi*sr*ssep*VminEK;

    INa  = GNa*sm*sm*sm*sh*sj*VminENa;
    IbNa = GbNa*VminENa;

    double ca_ss_total=sCaSS_total;
    double ca_ss=free_of_total(ca_ss_total,Bufss,Kbufss);

    ICaL  = sd*sf*sf2*sfcass*(GHKICaL1_abs*ca_ss-GHKICaL2_abs); /// SF changed   ICaL = GCaL*sd*sf*sf2*sfcass*(GHKICaL1*CaSS-GHKICaL2);
    IbCa  = GbCa*(svolt - rev_potential_ca(ca)); /// SF changed   IbCa  = GbCa*(svolt - (Ta->ECA[Caitab]));
    INaCa = ExpINaCa1_abs-ExpINaCa2_abs*ca; /// SF changed   INaCa = (ExpINaCa1*Nai*Nai*Nai-ExpINaCa2*Cai);
//    Ip    = ipmax*KoSat*(Nai/(Nai+KmNa))*RecIp;  /*=Inak from 2006 */ PK
    Ip    = RecIp_abs; /// SF changed /*=Inak from 2006 */    Ip = ipmax*(Ko/(Ko+KmK))*(Nai/(Nai+KmNa))*RecIp;
    IpCa  = GpCa*ca/(KpCa+ca); /// SF changed   IpCa  = GpCa*Cai/(KpCa+Cai);
    IpK   = RecIpK(svolt)*VminEK; /// SF changed  RecIpK*VminEK;


//            ISTIM=(i>=i_low && i<=i_high && j>=j_low && j<=j_high)?Istim:0;


    //Determine total current
    sItot_n = IK + IK1 + Ito + INa + IbNa + ICaL + IbCa + Ip + INaCa +
              IpCa + IpK + Istim;

    double ca_sr_total=sCaSR_total;
    double ca_sr=free_of_total(ca_sr_total,Bufsr,Kbufsr);

    kCaSR   = maxsr - (maxsr-minsr)/(1+sqr(EC)/sqr(ca_sr)); /// SF changed   kCaSR = maxsr - (maxsr-minsr)/(1+(EC/CaT)*(EC/CaT));
    dRR     = k4*(1-sRbar)-k2b*kCaSR*ca_ss*sRbar; /// SF changed   dRR = k4*(1-sRR)-k2b*kCaSR*CaSS*sRR;
    sRbar_n = sRbar + HT*dRR; /// SF changed   sRR_n = sRR + HT*dRR;

    double temp=(k1b/kCaSR)*sqr(ca_ss); /// SF added
    double o_gate = temp*sRbar_n/(k3+temp); /// SF changed   sO_n = ((k1b/kCaSR)*CaSS*CaSS*sRR_n)/(k3+(k1b/kCaSR)*CaSS*CaSS);
    Irel   = Vrel*o_gate*(ca_sr-ca_ss); /// SF changed   Irel = Vrel*sO_n*(CaT-CaSS);
    Ileak  = Vleak*(ca_sr-ca);  /// SF changed   Ileak = Vleak*(CaT-Cai);
    SERCA  = Vmaxup/(1.+(Kupsquare/sqr(ca))); /// SF changed   SERCA  = Vmaxup/(1.+(Kupsquare/(Cai*Cai)));

///    CaCSQN = Bufsr*CaT/(CaT+Kbufsr);  /*this is Ca_srbufsr*/ SF commented
    dCaT   = HT*(SERCA-Irel-Ileak);

    ca_sr_total += dCaT;
    sCaSR_total_n = ca_sr_total;

///    bjsr   = Bufsr /** -CaCSQN */  -dCaT-CaT+Kbufsr;  /// b = CA_BUF_SR - ca_sr_total + K_BUF_SR; /// SF commented
///    cjsr   = Kbufsr*(/**CaCSQN+ */dCaT+CaT); ///PK /// SF commented
///    CaT_n  = (sqrt(bjsr*bjsr+4*cjsr)-bjsr)*0.5; /// SF commented

    Jxfer   = Vxfer*(ca_ss-ca); /// SF changed   Jxfer   = Vxfer*(CaSS-Cai);
///    CaSSBuf = Bufss*CaSS/(CaSS+Kbufss); SF commented
    dCaSS   = HT*(-Jxfer*(vc/vss)+Irel*(vsr/vss)+(-ICaL*(1/(2*vss*F))*Cm));

    ca_ss_total += dCaSS;
    sCaSS_total_n = ca_ss_total;

///    bcss   = Bufss /** -CaSSBuf */ -dCaSS-CaSS+Kbufss;   ///  b = CA_BUF_SS - ca_ss_total + K_BUF_SS; /// SF commented
///    ccss   = Kbufss*(/**CaSSBuf+*/dCaSS+CaSS); ///PK /// SF commented
///    CaSS_n = (sqrt(bcss*bcss+4*ccss)-bcss)*0.5; /// SF commented

///    CaBuf = Bufc*Cai/(Cai+Kbufc); /// SF commented
    dCai  = HT*((-(IbCa+IpCa-2*INaCa)*(1/(2*vc*F))*Cm)+(Ileak-SERCA)*	(vsr/vc)+Jxfer);

    ca_total += dCai;
    sCai_total_n = ca_total;

///    bc = Bufc-CaBuf-dCai-Cai+Kbufc;
///    cc = Kbufc*(CaBuf+dCai+Cai);
///    Cai_n = (sqrt(bc*bc+4*cc)-bc)*0.5;

///    dNai = -(INa+IbNa+3*Ip+3*INaCa)*(1/(vc*F))*Cm;
///    Nai_n = Nai + HT*dNai;

///    axialcurrent=0;//diffusion*( Drhx*(volt((i+1),j)+ volt((i-1),j)-2.0*volt(i,j))+
//				     (0.25)*Drhy*(volt(i,(j+1))+volt(i,(j-1))-2.0*volt(i,j)));

///    dKi = -(///Istim+
///            IK1+Ito+IK-2*Ip+IpK-axialcurrent)*(1/(vc*F))*Cm;
///    Ki_n = Ki + HT*dKi;

    //Update gates
    sm_n = Ta->M_INF[Vtab]-(Ta->M_INF[Vtab]-sm)*Ta->EXPTAU_M[Vtab];
    sh_n = Ta->H_INF[Vtab]-(Ta->H_INF[Vtab]-sh)*Ta->EXPTAU_H[Vtab];
    sj_n = Ta->J_INF[Vtab]-(Ta->J_INF[Vtab]-sj)*Ta->EXPTAU_J[Vtab];
    sxr1_n = Ta->Xr1_INF[Vtab]-(Ta->Xr1_INF[Vtab]-sxr1)*Ta->EXPTAU_Xr1[Vtab];
    sxr2_n = Ta->Xr2_INF[Vtab]-(Ta->Xr2_INF[Vtab]-sxr2)*Ta->EXPTAU_Xr2[Vtab];
    sxs_n = Ta->Xs_INF[Vtab]-(Ta->Xs_INF[Vtab]-sxs)*Ta->EXPTAU_Xs[Vtab];
    ss_n = Ta->S_INF[Vtab]-(Ta->S_INF[Vtab]-ss)*Ta->EXPTAU_S[Vtab]; /// SF changed    ss_n = Ta->S_INF[Vtab]-(Ta->S_INF[Vtab]-ssep)*Ta->EXPTAU_Sepi[Vtab]; /// SF added
///    ssep_n= ss_n; /// SF changed Ta->S_INFepi[Vtab]-(Ta->S_INFepi[Vtab]-ssep)*Ta->EXPTAU_Sepi[Vtab];  /// SF commented
///    ssen_n= ss_n; /// SF changed Ta->S_INFendo[Vtab]-(Ta->S_INFendo[Vtab]-ssen)*Ta->EXPTAU_Sendo[Vtab]; /// SF commented
///    ssm_n= ss_n; /// SF changed Ta->S_INFm[Vtab]-(Ta->S_INFm[Vtab]-ssm)*Ta->EXPTAU_Sm[Vtab]; /// SF commented
    sr_n= Ta->R_INF[Vtab]-(Ta->R_INF[Vtab]-sr)*Ta->EXPTAU_R[Vtab];
    sd_n = Ta->D_INF[Vtab]-(Ta->D_INF[Vtab]-sd)*Ta->EXPTAU_D[Vtab];
    sf_n =Ta->F_INF[Vtab]-(Ta->F_INF[Vtab]-sf)* Ta->EXPTAU_F[Vtab]; /// SF decommented
    sf2_n =Ta->F2_INF[Vtab]-(Ta->F2_INF[Vtab]-sf2)* Ta->EXPTAU_F2[Vtab];

    temp=1./(1+sqr(ca_ss*20)); /// SF added
    taufcass=80*temp+2; /// SF changed    taufcass=80/(1+(CaSS_n*20)*(CaSS_n*20))+2; // /1000 removed
    fcassinf=0.6*temp+0.4; /// SF changed    fcassinf=0.6/(1+pow((CaSS_n*20),2))+0.4;
    sfcass_n = fcassinf-(fcassinf-sfcass)*exp(-HT/taufcass);

///    double Af = 2.45*450*exp(-(svolt+27)*(svolt+27)/225);
///    double Bf = 200/(1.+exp((13-svolt)*0.1));
///    double Cf = 180/(1.+exp((svolt+30)*0.1));
///    tauf=Af+Bf+Cf+20;
    //if(svolt>0) tauf=tauf*tau_f;
///    tauf=tauf*tau_f;//PK
///    finf=1./(1.+exp((svolt+20)/7));
///    sf_n=finf-(finf-sf)*exp(-HT/tauf);

    svolt_n = svolt + HT * (DivDGradU - sItot_n); /// SF changed   svolt_n = svolt + HT*(axialcurrent-sItot_n);
    /*in this equation capacitance is now removed*/
}

void CalcNewUV(double uv0[Nvar+1], double uv1[Nvar+1], double DivDGradU, double dt, double Istim)
{
    Step (uv0, uv1, &MyTabulation, DivDGradU, dt, Istim);
}

#endif /// TP06
