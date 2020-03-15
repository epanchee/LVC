/** @file
 *  @brief It is an intermediate layer between the main programme, which works with an abstract electrophysiological model, and specific electro-models.
 *
 *  For each model, we define its header file, some constants needed by the main program, initialization and
 *  main calculation functions.
 */

#ifndef ELECTRO_H_INCLUDED
#define ELECTRO_H_INCLUDED

// #define AP
// #define APmu
	// #define APmuRusakov
// #define TP06
	// #define slope07
	// #define slope11
#define LR

/// if curr_stim is defined, we use stimulation current; instantaneous stimulation otherwise
#define curr_stim
#ifdef curr_stim
extern double
    I_stim, t_stim, tipval;
#endif

#define sqr(x) ((x)*(x))
#define kub(x) ((x)*(x)*(x))

#ifdef TP06

/// number of phase variables except potential u
#define Nvar 16
extern double c_gNa, c_gCaL, c_verapamil;
extern int amiodaron; /// 0, no effect; 1, low concentration 1 micromol; 3, high conc. 3 micromol
extern double 
	TAU_F_INACT_FACTOR, Gkr, Gks, GpCa, GpK;
#endif /// TP06


#ifdef AP

extern const double InitialV;

extern double ka; //ka; a are constants of AP model
extern double a;
extern double eta;

/// number of phase variables except potential u
#define Nvar 1

#endif ///AP


#ifdef APmu

extern const double InitialV;

extern double ka; //ka; a are constants of AP model
extern double a;
extern double epsAP_v;
extern double mu1;
extern double mu2;

/// number of phase variables except potential u
#define Nvar 1

#endif ///APmu

#ifdef LR

/// number of phase variables except potential u
#define Nvar 7
extern double c_gNa; /// parameters for diminishing conductivities for INa, Isi, IK
extern double c_gsi;
extern double c_gK;

#endif /// LR

extern const int indV;  /// index of state variable voltage
extern const double volt_init; /// steady-state potential
extern const double volt_stim;

void CommonInit(double dt);
void InitUV(double uv[Nvar+1]);

void CalcNewUV(double uv0[Nvar+1],
               double uv1[Nvar+1],
               double DivDGradU, double dt
#ifdef curr_stim
                ,double Istim
#endif
               );

#endif // ELECTRO_H_INCLUDED
