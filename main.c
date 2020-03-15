#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "elektro.h"
#include "tipHans.h"

#define delta 0.001

int
    nx,ny; /// we have nodes (0...nx, 0...ny)

/// only one of the pacing modes has to be defined
///#define pacecondit
#define paceperiod
// #define S2time
#define S2cond
// #define S1S2_by_stimulation
#define obstacle

int electrode,
    chirality = 0; // 0 - clockwise SW, 1 - counterclockwise SW
double 
    electrode_width,
    S1_width, /// in percents of the field
    S2_width;
bool (*where)(int, int); 

bool Snow(int i, int j, bool (*where)(int, int)){
    return where(i,j);
}

bool S1(int i, int j)
{
    return chirality==0? i<S1_width*nx : j<S1_width*ny;
}

bool S2(int i, int j)
{
    return chirality==0? j<S2_width*ny : i<S2_width*nx;
}

bool S_last(int i, int j){
    // return i<S1_width*nx && j > 0.5*ny;
    return (i<10) && (j>ny-10);
}

#ifdef S2cond
bool forefront = false;
double h_limit_S2;
double l_limit_S2;
int S2_x = 200, S2_y = 200;
bool S2condition(double potential){
    if(potential > h_limit_S2 && !forefront)
        forefront = true;
    return (potential < l_limit_S2) && forefront;
}
#endif

/// this stimulus will be applied indefinitely many times with some frequency
bool Sn(int i, int j)
{
    if (electrode==1) return (i<10) && (j<10);
    if (electrode==2) return (fabs(i-0.5*nx)<5) && (fabs(j-0.5*ny)<5);
    if (electrode==3) return (i<10) && ((j<10) || (ny-j<10));
    if (electrode==4) return fabs(i-0.5*nx)<5 && ((j<10) || (ny-j<10));
    if (electrode==5) return (i<=10);
    if (electrode==6) return (i<=20);
    if (electrode==100) return (i<=(int)electrode_width);
    if (electrode==101) return (i>=nx-10);
    return false;
}

double pace_start, pace_stop, /// pacing starts after ... and is stopped at ...
       S_last_t=-1;

#ifdef paceperiod
double pace_period;
///int pace_phase;
#endif // paceperiod
int
    space_mult_x,space_mult_y, ///print every ... point
    print_every, ///print data and seek tips every ... timesteps
    scar_stencil_i, scar_stencil_j; // left bottom corner of scar stencil
long int
    pos; // position in input file
double
    dt, /// time step     !!! dt*print_every should be integer !!!
    t1, /// simulation time
    t_S2,
    drx, dry,
    Diffuz1, /// Diffusion coefficient along fibres
    Diffuz2; /// across fibres
// #ifdef smartstim
double
    u_low = -50,
    u_high; /// low and high limits for smart stimulation
bool
    good; /// boolean variable for smart stimulation
int
    smartstim;
FILE *stimout;
// #endif
#ifdef TP06
    double slope = .7;
#endif // TP06

#ifdef AP

#endif // AP

#ifdef APmu

#endif // APmu

#ifdef TP06

#endif // TP06

#ifdef LR

#endif /// LR || TP06

void *mycalloc(size_t number, size_t size)
{
    void *res=calloc(number,size);
    if (res==NULL)
    {
        printf ("error! cannot allocate memory %zu x %zu\n",number,size);
        exit (1);
    }
    return res;
}

///Allocate memory for a (0..n1 x 0..n2) real array
void CreateMatrix2(double*** M, int n1, int n2)
{
    *M=(double**)mycalloc((size_t) (n1 + 1), sizeof(double*));
    int i;
    for (i=0; i<=n1; i++)
        (*M)[i]=(double*)mycalloc((size_t) (n2 + 1), sizeof(double));
}

///Allocate memory for a (0..n1 x 0..n2) real array
void CreateIntMatrix2(int*** M, int n1, int n2)
{
    *M=(int**)mycalloc((int) (n1 + 1), sizeof(int*));
    int i;
    for (i=0; i<=n1; i++)
        (*M)[i]=(int*)mycalloc((int) (n2 + 1), sizeof(int));
}

///Deallocate memory of an array (0..n1 * 0..n2), it doesn't depend on n2
void DestroyMatrix2(void** M, int n1)
{
    int i;
    for (i=0; i<=n1; i++)
        free (M[i]);
    free (M);
}
///Allocate memory for a (0..n1 x 0..n2 x 0..n3) real array
void CreateMatrix3(double**** M, int n1, int n2, int n3)
{
    *M=(double***)mycalloc((size_t) (n1 + 1), sizeof(double**));
    int i,j;
    for (i=0; i<=n1; i++)
    {
        (*M)[i]=(double**)mycalloc((size_t) (n2 + 1), sizeof(double*));
        for (j=0; j<=n2; j++)
            (*M)[i][j]=(double*)mycalloc((size_t) (n3 + 1), sizeof(double));
    }
}
///Deallocate memory of an array (0..n1 * 0..n2 * 0..n3), it doesn't depend on n3
void DestroyMatrix3(void*** M, int n1, int n2)
{
    int i,j;
    for (i=0; i<=n1; i++)
    {
        for (j=0; j<=n2; j++)
            free (M[i][j]);
        free (M[i]);
    }
    free (M);
}

void init_boundary_mask(int** mask, int nx, int ny){
    for (int i = 0; i <= nx; ++i)
    {
        for (int j = 0; j <= ny; ++j)
        {
            mask[i][j] = 0;
        }
    }

    for (int i = 0; i <= nx; ++i)
    {
        mask[i][0] = 1;
        mask[i][ny] = 1;
    }

    for (int j = 0; j <= ny; ++j)
    {
        mask[0][j] = 1;
        mask[nx][j] = 1;
    }
}

int main(void)
{
    FILE *fin=fopen("input.txt","r");
    if (fin==NULL)
    {
        printf ("error! cannot open input.txt\n");
        exit (1);
    }
    #define readint(var) if(fscanf(fin,#var "=%d ",&var)!=1) {printf("error during reading " #var "!\n");exit(2);} else printf(#var "=%d\n",var);
    #define readdbl(var) if(fscanf(fin,#var "=%lf ",&var)!=1) {printf("error during reading " #var "!\n");exit(2);} else printf(#var "=%f\n",var);
    #define readint_with_default(var,defval) pos=ftell(fin);  if(fscanf(fin,#var "=%d ",&var)!=1) {printf("cannot read " #var ", set default value %i\n",defval); var=defval; fseek(fin, pos, SEEK_SET);} else printf(#var "=%d\n",var);
    #define readdbl_with_default(var,defval) pos=ftell(fin); if(fscanf(fin,#var "=%lf ",&var)!=1) {printf("cannot read " #var ", set default value %f\n",defval); var=defval; fseek(fin, pos, SEEK_SET);} else printf(#var "=%f\n",var);

    readdbl(pace_start)
#ifdef paceperiod
    readdbl(pace_period)
///    readint(pace_phase)
#endif /// paceperiod
    readint(space_mult_x)
    readint(space_mult_y)
    readint(print_every)
    readdbl(dt)
    readdbl(t1)
    readdbl(drx)
    readdbl(dry)
    readint(nx)
    readint(ny)
    nx+=2;
    ny+=2;
#ifdef S2time
    readdbl(t_S2)
#endif
#ifdef S2cond
    readint(S2_x)
    double S2_y_default = (ny - 1) / 2;
    readint_with_default(S2_y, S2_y_default)
    readdbl_with_default(h_limit_S2, -13.0)
    readdbl_with_default(l_limit_S2, -50.0)
#endif
    readint_with_default(chirality, 0)
    readdbl(Diffuz1)
    readdbl(Diffuz2)
    readdbl_with_default(S1_width, 0.1)
    readdbl_with_default(S2_width, 0.5)
#ifdef curr_stim
    readdbl(I_stim)
    readdbl(t_stim)
#endif /// curr_stim

/*#ifdef obstacle
    readint(scar_stencil_i)
    readint(scar_stencil_j)
#endif*/

// #ifdef smartstim
readint_with_default(smartstim, 0)
readdbl_with_default(u_low, -50)
    // readdbl(u_high)
// #endif /// smartstim

#ifdef AP
    readdbl(ka)
    readdbl(a)
    readdbl(eta)

#endif // AP

#ifdef APmu
    readdbl(ka)
    readdbl(a)
    readdbl(epsAP_v)
    readdbl(mu1)
    readdbl(mu2)

#endif // APmu

#ifdef LR
    readdbl_with_default(c_gNa,1.)
    readdbl_with_default(c_gsi,1.)
    readdbl_with_default(c_gK,1.)
#endif // LR

#ifdef TP06
    readdbl_with_default(slope,.7)
    readdbl_with_default(c_gNa,1.)
    readdbl_with_default(c_gCaL,1.)
    readint_with_default(amiodaron,0)
    readdbl_with_default(c_verapamil,1.)
#endif // TP06

    readint(electrode)
    if (electrode==100)
        readdbl(electrode_width)
    readdbl_with_default(pace_stop,t1)
    readdbl_with_default(S_last_t,-1)

    int assim;
    readint_with_default(assim,0)
     /// we give Nbegin periods equal Tbegin first, then we decrease period to Ttarget
     /// we decrease by deltaT at a time
     /// we give Nrepeat stimuli for each series

     /// example: Nrepeat=3, deltaT=-1, Tbegin=100, Ttarget=96, Nbegin=7
     /// 100 100 100 100 100 100 100 99 99 99 98 98 98 97 97 97 96 96 96
    int Nbegin,Nrepeat;
    double Tbegin=pace_period,Ttarget,deltaT;
    readint(Nbegin)
    readdbl(deltaT)
    readdbl(Ttarget)
    readint(Nrepeat)

    int read_from_file,save_to_file;
    readint_with_default(read_from_file,0)
    readint_with_default(save_to_file,0)

    #undef readdbl
    #undef readint
    fclose(fin);

    #ifdef S2cond
        printf("Using S2 by condition\n");
    #endif

    #ifdef TP06
        if (slope == 0.7){
            TAU_F_INACT_FACTOR = .6;
            Gkr = 0.134;
            Gks = 0.270;
            GpCa = 0.0619;
            GpK = 0.0730;
        } else if (slope == 1.1){
            TAU_F_INACT_FACTOR = 1.;
            Gkr = 0.153;
            Gks = 0.392;
            GpCa = 0.1238;
            GpK = 0.0146;
        } else if (slope == 1.4){
            TAU_F_INACT_FACTOR = 1.5;
            Gkr = 0.172;
            Gks = 0.441;
            GpCa = 0.3714;
            GpK = 0.0073;
        } else if (slope == 1.8){
            TAU_F_INACT_FACTOR = 2.;
            Gkr = 0.172;
            Gks = 0.441;
            GpCa = 0.8666;
            GpK = 0.00219;
        } else {
            printf ("There are no parameters for this slope. Exiting...\n");
            return 0;            
        }

    #endif

    int step_num=0; // t=step_num*dt
    char *fnu="u.dat";

    if ((drx<0) || (dry<0))
    {
        printf ("Bad input: drx or dry<0\n");
        return 0;
    }
    if ((nx<0) || (ny<0))
    {
        printf ("Bad input: nx or ny<0\n");
        return 0;
    }
    if (t1<0)
    {
        printf ("Bad input: t1<0\n");
        return 0;
    }
    if (dt<0)
    {
        printf ("Bad input: dt<0\n");
        return 0;
    }

    FILE *fparam=fopen("param.dat","w");
    fprintf (fparam,"%f\n%f\n%i\n%f\n%f\n%i \n%f\n%i\n%i\n",
             dt*print_every,
             drx,
             nx-2,
             Diffuz1,
             dt,
             space_mult_x,
             dry, ny-2, space_mult_y);

    FILE *fout=fopen(fnu,"w");
    if (fout==NULL)
    {
        printf ("error! cannot open file %s for output\n",fnu);
        return 0;
    }

    FILE *f_stimtable=fopen("stim_table.info","w");
    if (f_stimtable==NULL)
    {
        printf ("error! cannot open file stim_table.info for stimulus field output\n");
        return 0;
    }

    if(smartstim == 1) {
        stimout=fopen("stimulus.info","w");
        if (stimout==NULL)
        {
            printf ("error! cannot open file stimulus.info for stimulus field output\n");
            return 0;
        }
    }

    double start_time, end_time;
    int i,j;
    double ***uv0, ***uv1, ***temp;
    double **st;
    int **mask;

    double **u_curr, **u_very_old;
    CreateMatrix2 (&u_curr,nx,ny);
    CreateMatrix2 (&u_very_old,nx,ny);

    CreateMatrix3 (&uv0,nx,ny,Nvar);
    CreateMatrix3 (&uv1,nx,ny,Nvar);

    CreateMatrix2 (&st,nx,ny);
    CreateIntMatrix2 (&mask,nx,ny);
    init_boundary_mask(mask, nx, ny);

    CommonInit(dt);

    #ifdef obstacle

    FILE *f_scar=fopen("scar.dat","r");
    if (!f_scar)
    {
        printf ("Error while opening file scar.dat! Make sure it exists\n");
        exit(1);
    }
    int scar_stencil_size_i, scar_stencil_size_j;
    fscanf (f_scar, "%d %d %d %d", &scar_stencil_j, &scar_stencil_i, &scar_stencil_size_j, &scar_stencil_size_i);

    if (scar_stencil_i+scar_stencil_size_i >= nx || scar_stencil_j+scar_stencil_size_j >= ny){
        printf ("Scar size exceeds field size! Exiting ...\n");
        exit(1);
    }
    for (i=scar_stencil_i; i<scar_stencil_i+scar_stencil_size_i; i++){
        for (j=scar_stencil_j; j<scar_stencil_j+scar_stencil_size_j; j++){
            if (fscanf (f_scar, "%d ", &mask[i][j]) != 1){
                printf("Error during reading scar mask\n");
                exit(2);
            }
        }
    }

    fclose(f_scar);    

    #endif

    if (read_from_file == 1)
    {
        FILE *fin=fopen("uv0.bin","rb");
        if (!fin)
        {
            printf ("Error while opening input binary file uv0.bin! Make sure it exists\n");
            exit(1);
        }
        for (i=0; i<=nx; i++)
            for (j=0; j<=ny; j++)
                if (fread(uv0[i][j],sizeof(double),Nvar+1,fin)!=Nvar+1)
                {
                    printf ("Error while reading input binary file! i=%i j=%i\n",i,j);
                    exit(1);
                }
        fclose(fin);
    }
    else
        for (i=0; i<=nx; i++)
            for (j=0; j<=ny; j++)
            {
                InitUV (uv0[i][j]);
                #if !defined(S1S2_by_stimulation)
                /// S1
                if (Snow(i,j, S1))
                    uv0[i][j][indV]=volt_stim;
                #endif
            }

    FILE *ftip;
    if ((ftip=fopen("tipsHans.info","w"))==NULL)
    {
        printf ("error! cannot open file tipsHans.txt for writing!\n");
        exit (1);
    }

    /// First- and second-order derivatives of U (div(D*grad u) == axby*Ux + bxcy*Uy + A Uxx + 2B Uxy + C Uyy
    double u_x,u_y,u_xx,u_xy,u_yy, DivDGradU;

    double _4drxdry = 1./(4.*drx*dry);
    double _sqrdrx=1./(drx*drx);
    double _sqrdry=1./(dry*dry);
    double _2drx=1./(2.*drx);
    double _2dry=1./(2.*dry);

    double tipdata[maxntips+1][2];
    bool _S2=false;
    bool fst_stim_tstep=false; // for stimulus field output
#ifdef curr_stim
    bool pacingNow=false;
#endif
    start_time = omp_get_wtime();

    /// we are sure that Nstim stimulations are enough for all our simulations
    #define Nstim 10000
    double periods[Nstim], moment_stim[Nstim+1], stim_end=0;
    moment_stim[0]=pace_start;

    if (assim == 1)
    {
        int Ntransition=(int)((Ttarget-Tbegin)/deltaT)*Nrepeat;
        for (i=0; i<Nstim; i++)
            if (i<Nbegin)
                periods[i] = Tbegin;
            else
            {
                if (i<Nbegin+Ntransition)
                {
                    if ((i-Nbegin)%Nrepeat==0)
                        periods[i] = periods[i-1]+deltaT;
                    else
                        periods[i] = periods[i-1];
                }
                else
                    periods[i] = Ttarget;
            }
    }
    else
    {
        for (i=0; i<Nstim; i++)
            periods[i] = pace_period;
    }

    for (i=1; i<=Nstim && moment_stim[i-1] + periods[i-1] < pace_stop; i++)
        moment_stim[i] = moment_stim[i-1] + periods[i-1];

    int stim_number = i;

    for (i=0; (i<stim_number) && (moment_stim[i]<t1) && (moment_stim[i]<pace_stop); i++){
        printf ("%i period=%.1f moment_stim=%.1f\n",i,periods[i],moment_stim[i]);
        fprintf(f_stimtable, "%i %.1f %.1f\n",i,periods[i],moment_stim[i]);
    }
    fflush(f_stimtable);
    fclose (f_stimtable);

    int stim_num=0;

    while (step_num*dt<t1)
    {
        // S1
        #ifdef S1S2_by_stimulation
        if (!pacingNow && step_num == 0 && read_from_file == 0 && S1_width > 0)
        {
            pacingNow = true;
            where = &S1;
            stim_end = t_stim;
        }
        #endif

        /// S2
        if (!_S2 && read_from_file == 0 && electrode != 100)
            #ifdef S2time
            if (step_num*dt > t_S2)
            #endif
            #ifdef S2cond
            if (S2condition(uv0[S2_x][S2_y][indV]))
            #endif
            {
                printf("Applying an S2 stimulus\n");

                #ifdef S1S2_by_stimulation
                // S2
                if (!pacingNow && S2_width > 0)
                {
                    pacingNow = true;
                    where = &S2;
                    stim_end = step_num*dt+t_stim;
                }
                #else
                for (i=0; i<=nx; i++)
                    for (j=0; j<=ny; j++)
                        if (Snow(i,j,S2))
                            uv0[i][j][indV]=volt_stim;
                #endif

                _S2=true;
            }


        /// pacing
        if (step_num*dt > moment_stim[stim_num] && step_num*dt < pace_stop && stim_num < stim_number)
        {
#ifdef curr_stim
            if (!pacingNow)
            {
                pacingNow=true;
                where = &Sn;
                fst_stim_tstep = true;
                stim_end = moment_stim[stim_num] + t_stim;
                printf ("+ time=%f, period=%f, stim_end=%f, stim_num=%i\n",step_num*dt,periods[stim_num],stim_end,stim_num);
                stim_num++;
            }
#else
            for (i=0; i<=n; i++)
                for (j=0; j<=n; j++)
                    if (Snow(i,j, Sn))
                        uv0[i][j][indV]=volt_stim;
            stim_num++;
#endif
        }

        /// last stim
        if (step_num*dt > S_last_t && step_num*dt < S_last_t + t_stim){
            if (!pacingNow)
            {
                pacingNow=true;
                where = &S_last;
                stim_end = S_last_t + t_stim;
            }
        }

        /// File output
        if (step_num%print_every==0)
        {
            for (i=0; i<=nx; i++)
                for (j=0; j<=ny; j++)
                    u_curr[i][j] = uv0[i][j][indV];

            if (step_num>0)
            {
                double tipvals[2];
                tipvals[0] = tipval; /// for old voltage
                tipvals[1] = tipval; /// for new voltage

                int tipsfound=0;

                track_tipline (nx, ny, u_very_old, u_curr, tipvals, tipdata, &tipsfound);

                printf ("Found %i tip(s)\n",tipsfound); fflush(stdout);
                if (tipsfound>0)
                {
                    for (i=0; i<tipsfound; i++)
                    {
                        fprintf (ftip,"%.0f %i ",step_num*dt,i);

                        fprintf (ftip,"%5.2f %5.2f\n",tipdata[i][0]*drx,tipdata[i][1]*dry);
                    }
                    fflush (ftip);
                }
            }
            for (i=0; i<=nx; i++)
                memcpy (u_very_old[i], u_curr[i], sizeof(double)*(ny+1));

            for (i=1; i<nx; i+=space_mult_x)
            {
                for (j=1; j<ny; j+=space_mult_y)
                    fprintf (fout,"%.0f ",uv0[i][j][indV]
#if defined(AP) || defined(APmu)
                             *100
#endif
                             );
                fprintf (fout,"\n");
            }
            fflush(fout);

            printf("%6.2f%%\n",100.0*step_num*dt/t1); fflush (stdout);
        }

        double ux,uy;

        bool u_xx_done, u_yy_done;

        #pragma omp parallel for private(i, j, u_xx, u_yy, u_xx_done, u_yy_done, DivDGradU)
        for (i=1; i<nx; i++)
        {
            for (j=1; j<ny; j++)
            {

                /// we do not calculate potential on border or inside obstacle
                if (mask[i][j] != 0){
                    uv1[i][j][indV] = volt_init;
                    continue;
                }

                u_xx_done = u_yy_done = false;

                if (mask[i-1][j] == 1){
                    u_xx = (2*uv0[i+1][j][indV] - 2*uv0[i][j][indV])*_sqrdrx;
                    u_xx_done = true;
                }
                if (mask[i+1][j] == 1){
                    u_xx = (2*uv0[i-1][j][indV] - 2*uv0[i][j][indV])*_sqrdrx;
                    u_xx_done = true;
                }
                if (mask[i][j+1] == 1){
                    u_yy = (2*uv0[i][j-1][indV] - 2*uv0[i][j][indV])*_sqrdry;
                    u_yy_done = true;
                }
                if (mask[i][j-1] == 1){
                    u_yy = (2*uv0[i][j+1][indV] - 2*uv0[i][j][indV])*_sqrdry;
                    u_yy_done = true;
                }

                if(!u_xx_done)
                    u_xx=(uv0[i-1][j][indV] - 2*uv0[i][j][indV] + uv0[i+1][j][indV])*_sqrdrx;
                if(!u_yy_done)
                    u_yy=(uv0[i][j-1][indV] - 2*uv0[i][j][indV] + uv0[i][j+1][indV])*_sqrdry;
                DivDGradU = Diffuz1*u_xx + Diffuz2*u_yy;

                if(smartstim == 1){
                    CalcNewUV(uv0[i][j], uv1[i][j], DivDGradU, dt, 0);
                    if(pacingNow && Snow(i,j,where)){
                        // good = (uv0[i][j][indV] < u_low || uv0[i][j][indV] > u_high || uv1[i][j][indV] - uv0[i][j][indV] > 0);
                        if(fst_stim_tstep){
                            good = (uv0[i][j][indV] < u_low || uv1[i][j][indV] - uv0[i][j][indV] > 0);
                            st[i][j] = good? 1 : 0;
                        }

                        if(st[i][j] == 1)
                            #ifdef LR
                            uv1[i][j][indV] = I_stim*dt;
                            #else
                            uv1[i][j][indV] -= I_stim*dt;
                            #endif
                    }
                } else
                CalcNewUV(uv0[i][j], uv1[i][j], DivDGradU, dt
#ifdef curr_stim
                                                                ,pacingNow ? (Snow(i,j,where) ? I_stim : 0) : 0
#endif
                          );             

                if (fpclassify(uv1[i][j][indV])==FP_NAN)
                {
                    printf ("\n NAN! \n");
                    exit (1);
                }
            }
        }

        temp=uv1;
        uv1=uv0;
        uv0=temp;

if (smartstim == 1)
{
    if (pacingNow && fst_stim_tstep){
        // stimulus field output
        fst_stim_tstep = false;
        for (i=1; i<nx; i++)
        {
            for (j=1; j<ny; j++)
                fprintf (stimout,"%.0f ",st[i][j]);
            fprintf (stimout,"\n");
        }
        fflush(stimout);
    }
}

#ifdef curr_stim
        if (pacingNow && (step_num*dt > stim_end))
        {
            pacingNow=false;
            printf ("- stim end=%f, time=%f\n",stim_end,step_num*dt);
        }
#endif /// curr_stim

        step_num++;
    }

    end_time=omp_get_wtime();
    printf ("Time of calculation = %f ms\n", end_time - start_time);
    fprintf(fparam, "%i\n", 1); // calculation is succesfully done
    fprintf(fparam, "%i\n", (int) end_time - start_time); // time of calculation

    fclose (fout);
    fclose (fparam);
    if(smartstim == 1) fclose(stimout);

    if (save_to_file == 1)
    {
        fout=fopen("uv1.bin","wb");
        if (!fout)
        {
            printf ("Error while opening output binary file uv1.bin!\n");
            exit(1);
        }
        for (i=0; i<=nx; i++)
            for (j=0; j<=ny; j++)
                if (fwrite(uv0[i][j],sizeof(double),Nvar+1,fout)!=Nvar+1)
                {
                    printf ("Error while writing input binary file! i=%i j=%i\n",i,j);
                    exit(1);
                }
        fclose (fout);
    }

    DestroyMatrix3 ((void***)uv0,nx,ny);
    DestroyMatrix3 ((void***)uv1,nx,ny);
#ifndef staticArr
    DestroyMatrix2 ((void**)u_curr,nx);
    DestroyMatrix2 ((void**)u_very_old,nx);
#endif
    DestroyMatrix2 ((void**)st,nx);
    DestroyMatrix2 ((void**)mask,nx);

    return 0;
}
