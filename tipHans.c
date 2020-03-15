// class structure for medium geometry
// by Hans Dierckx 4 Jan 2011

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include "tipHans.h"

// calculate tip position at subgrid level (using bilinear interpolation); based on calculate.cc in tip.cc
int calc_tippos(double vij, double vi1j, double vi1j1, double vij1, double vnewij, double vnewi1j,
                double vnewi1j1, double vnewij1, double V_iso1, double V_iso2, double *xy)
{
    double AC, BC, GC, DC;
    double AD, BD, GD, DD;
    double Disc;
    double Q, R, S, T, U, V;
    double xn, xp, yn, yp;
    double QOnR, SOnR;
    double T2;
    double sqrtDisc;

    //For old voltage values of point and neighbours
    AC=((vij)-(vij1)+(vi1j1)-(vi1j));//upleft-downleft + downright-upright
    GC=((vij1)-(vij));                 //downleft-upleft
    BC=((vi1j)-(vij));                //upright-upleft
    DC=(vij-V_iso1); // adapted here  //upleft-iso

    //For current voltage values of point and neighbours
    AD=((vnewij)-(vnewij1)+(vnewi1j1)-(vnewi1j));
    GD=((vnewij1)-(vnewij));
    BD=((vnewi1j)-(vnewij));
    DD=(vnewij-V_iso2); // adapted here

    Q=(BC*AD-BD*AC);
    R=(GC*AD-GD*AC);
    S=(DC*AD-DD*AC);

    QOnR=Q/R;
    SOnR=S/R;

    T=AC*QOnR;
    U=(AC*SOnR-BC+GC*QOnR);
    V=(GC*SOnR)-DC;

    //Compute discriminant of the abc formula
    //with a=T, b=U and c=V
    Disc=U*U-4.*T*V;
    if(Disc<0)
    {
        //If the discriminant is smaller than
        //zero there is no solution and a
        //failure flag should be returned
        return(0);
    }
    else
    {


// Otherwise two solutions for xvalues
        T2=2. * T;
        sqrtDisc=sqrt(Disc);

        xn=(-U-sqrtDisc)/T2;
        xp=(-U+sqrtDisc)/T2;
        //Leading to two solutions for yvalues
        yn=-QOnR*xn-SOnR;
        yp=-QOnR*xp-SOnR;

// demand that fractions lie in interval [0,1]
        if(xn>=0 && xn<=1 && yn>=0 && yn<=1)
        {
            //If the first point fulfills these
            //conditions take that point
            xy[0]=xn;
            xy[1]=yn;
            return(1);
        }
        else if(xp>=0 && xp<=1 && yp>=0 && yp<=1)
        {
            //If the second point fulfills these
            //conditions take that point
            xy[0]=xp;
            xy[1]=yp;
            return(1);
        }
        else
        {
            //If neither point fulfills these
            //conditions return a failure flag.
            return(0);
        }
    }
}

void track_tipline(int nx, int ny, double **var1, double **var2,
                   double* tipvals, double tipdata_[maxntips+1][2], int *tipsfound_)
{
    int tipsfound;

#define tipdata(i,j) tipdata_[i][j]
#define u(volt, ind1, ind2) (volt[(ind1)][(ind2)])

    double iso1 = tipvals[0];
    double iso2 = tipvals[1];


// each row contains data for a point of the tipline
// stored in columns: type, posx posy posz, dxu dyu dzu, dxv dyv dzv, tx ty tx
// possible types: 1 = YZ, 11 = YZ at lower medium boundary, 21 = YZ at upper medium boundary, 31 at lower domain boundary, 41 at upper domain boundary
//                 2 = XZ, 12 = XZ at lower medium boundary, 22 = XZ at upper medium boundary
//                 3 = XY, 13 = XY at lower medium boundary, 23 = XY at upper medium boundary
// sign of type equals < T , e_i > (i.e. entering or leaving the voxel)

    tipsfound = 0;
    int counter,xpos,ypos, delta=10;

    int size[2]= {nx,ny};

    // check XY-planes
    for (xpos = delta; xpos < size[0]-delta; xpos++)
    {
        for (ypos = delta; ypos < size[1]-delta; ypos++)
        {
            if (tipsfound >= maxntips)
                goto gotovo;
            counter=0;
            if (      (u(var1, xpos, ypos) >= iso1) &&
                    ( (u(var1, xpos+1, ypos) < iso1) ||
                      (u(var1, xpos, ypos+1) < iso1) ||
                      (u(var1, xpos+1, ypos+1) < iso1))     )
                counter=1;
            else
                if  ( (u(var1, xpos, ypos) < iso1) &&
                      ((u(var1, xpos+1, ypos) >= iso1) ||
                       (u(var1, xpos, ypos+1) >= iso1) ||
                       (u(var1, xpos+1, ypos+1) >= iso1))     )
                    counter=1;

            if(counter==1)
            {
                if ( (u(var2, xpos, ypos) >= iso2) &&
                        ((u(var2, xpos+1, ypos) < iso2) ||
                         (u(var2, xpos, ypos+1) < iso2) ||
                         (u(var2, xpos+1, ypos+1) < iso2))     )
                    counter=2;
                else
                    if  ( (u(var2, xpos, ypos) < iso2) &&
                          ((u(var2, xpos+1, ypos) >= iso2) ||
                           (u(var2, xpos, ypos+1) >= iso2) ||
                           (u(var2, xpos+1, ypos+1) >= iso2))     )
                        counter=2;

                if (counter==2)
                {
                    double tipfco [2];
                    if(calc_tippos(u(var1, xpos, ypos), u(var1, xpos+1, ypos), u(var1, xpos+1, ypos+1), u(var1, xpos, ypos+1),
                                   u(var2, xpos, ypos), u(var2, xpos+1, ypos), u(var2, xpos+1, ypos+1), u(var2, xpos, ypos+1), iso1, iso2, tipfco)==1)
                    {
                        tipdata(tipsfound, 0) = xpos+tipfco[0];
                        tipdata(tipsfound, 1) = ypos+tipfco[1];

                        tipsfound++;

                    } // if calc_tippos
                } // if counter 2
            } // if counter 1
            // end loop ZXY
        }
    }

gotovo: ;
    if(tipsfound>maxntips-2)
    {
        printf ("Warning from tracktipline.m: maximum number of tips tracked reached:%i\n",maxntips);
    }
    *tipsfound_ = tipsfound;
}
