#ifndef _DERIVATIVES_H

/* Second order accurate first-derivative operators */
#define D2CEN(f,i,di,h) ((-f[-di + i] + f[di + i])*oo2##h)


/* Fourth order accurate first-derivative operators */
#define D4CEN(f,i,di,h) ((f[-2*di + i] - 8*f[-di + i] + 8*f[di + i] -\
     f[2*di + i])*oo12##h)

/* Sixth order accurate first-derivative operators */
#define D6CEN(f,i,di,h) ((-f[-3*di + i] + 9*f[-2*di + i] -\
     9*(5*f[-di + i] - 5*f[di + i] + f[2*di + i]) +\
     f[3*di + i])*oo60##h)


/* Eighth order accurate first-derivative operators */
#define D8CEN(f,i,di,h) ((3*f[-4*di + i] - 32*f[-3*di + i] +\
     168*(f[-2*di + i] - 4*f[-di + i] + 4*f[di + i] -\
     f[2*di + i]) + 32*f[3*di + i] - 3*f[4*di + i])*oo840##h)


/* Second order accurate second-derivative operators */
#define D2xy(f,i,di,dj,dx,dy) ((f[-di - dj + i] - f[di - dj + i] -\
     f[-di + dj + i] + f[di + dj + i])*oo4##dx##dy)
#define D2xx(f,i,di,dx) ((-2*f[i] + f[-di + i] +\
     f[di + i])*oo1##dx##2)


/* Fourth order accurate second-derivative operators */
#define D4xy(f,i,di,dj,dx,dy) ((-8*f[-di - 2*dj + i] +\
     8*f[di - 2*dj + i] - f[2*di - 2*dj + i] -\
     8*f[-2*di - dj + i] + 64*f[-di - dj + i] -\
     64*f[di - dj + i] + 8*f[2*di - dj + i] +\
     8*f[-2*di + dj + i] - 64*f[-di + dj + i] +\
     64*f[di + dj + i] - 8*f[2*di + dj + i] -\
     f[-2*di + 2*dj + i] + 8*f[-di + 2*dj + i] -\
     8*f[di + 2*dj + i] + f[-2*(di + dj) + i] +\
     f[2*(di + dj) + i])*oo144##dx##dy)
#define D4xx(f,i,di,dx) ((-30*f[i] - f[-2*di + i] + 16*(f[-di + i] +\
     f[di + i]) - f[2*di + i])*oo12##dx##2)

/* Sixth order accurate second-derivative operators */
#define D6xy(f,i,di,dj,dx,dy) ((-9*f[-2*di - 3*dj + i] +\
     45*f[-di - 3*dj + i] - 45*f[di - 3*dj + i] +\
     9*f[2*di - 3*dj + i] - f[3*di - 3*dj + i] -\
     9*f[-3*di - 2*dj + i] - 405*f[-di - 2*dj + i] +\
     405*f[di - 2*dj + i] - 81*f[2*di - 2*dj + i] +\
     9*f[3*di - 2*dj + i] + 45*f[-3*di - dj + i] -\
     405*f[-2*di - dj + i] + 2025*f[-di - dj + i] -\
     2025*f[di - dj + i] + 405*f[2*di - dj + i] -\
     45*f[3*di - dj + i] - 45*f[-3*di + dj + i] +\
     405*f[-2*di + dj + i] - 2025*f[-di + dj + i] +\
     2025*f[di + dj + i] - 405*f[2*di + dj + i] +\
     45*f[3*di + dj + i] + 9*f[-3*di + 2*dj + i] -\
     81*f[-2*di + 2*dj + i] + 405*f[-di + 2*dj + i] -\
     405*f[di + 2*dj + i] - 9*f[3*di + 2*dj + i] -\
     f[-3*di + 3*dj + i] + 9*f[-2*di + 3*dj + i] -\
     45*f[-di + 3*dj + i] + 45*f[di + 3*dj + i] -\
     9*f[2*di + 3*dj + i] + f[-3*(di + dj) + i] +\
     81*(f[-2*(di + dj) + i] + f[2*(di + dj) + i]) +\
     f[3*(di + dj) + i])*oo3600##dx##dy)
#define D6xx(f,i,di,dx) ((-490*f[i] + 2*f[-3*di + i] -\
     27*(f[-2*di + i] - 10*(f[-di + i] + f[di + i]) +\
     f[2*di + i]) + 2*f[3*di + i])*oo180##dx##2)



/* Eighth order accurate second-derivative operators */
#define D8xy(f,i,di,dj,dx,dy) ((-96*f[-3*di - 4*dj + i] +\
     504*f[-2*di - 4*dj + i] - 2016*f[-di - 4*dj + i] +\
     2016*f[di - 4*dj + i] - 504*f[2*di - 4*dj + i] +\
     96*f[3*di - 4*dj + i] - 9*f[4*di - 4*dj + i] -\
     96*f[-4*di - 3*dj + i] - 5376*f[-2*di - 3*dj + i] +\
     21504*f[-di - 3*dj + i] - 21504*f[di - 3*dj + i] +\
     5376*f[2*di - 3*dj + i] - 1024*f[3*di - 3*dj + i] +\
     96*f[4*di - 3*dj + i] + 504*f[-4*di - 2*dj + i] -\
     5376*f[-3*di - 2*dj + i] - 112896*f[-di - 2*dj + i] +\
     112896*f[di - 2*dj + i] - 28224*f[2*di - 2*dj + i] +\
     5376*f[3*di - 2*dj + i] - 504*f[4*di - 2*dj + i] -\
     2016*f[-4*di - dj + i] + 21504*f[-3*di - dj + i] -\
     112896*f[-2*di - dj + i] + 451584*f[-di - dj + i] -\
     451584*f[di - dj + i] + 112896*f[2*di - dj + i] -\
     21504*f[3*di - dj + i] + 2016*f[4*di - dj + i] +\
     2016*f[-4*di + dj + i] - 21504*f[-3*di + dj + i] +\
     112896*f[-2*di + dj + i] - 451584*f[-di + dj + i] +\
     451584*f[di + dj + i] - 112896*f[2*di + dj + i] +\
     21504*f[3*di + dj + i] - 2016*f[4*di + dj + i] -\
     504*f[-4*di + 2*dj + i] + 5376*f[-3*di + 2*dj + i] -\
     28224*f[-2*di + 2*dj + i] + 112896*f[-di + 2*dj + i] -\
     112896*f[di + 2*dj + i] - 5376*f[3*di + 2*dj + i] +\
     504*f[4*di + 2*dj + i] + 96*f[-4*di + 3*dj + i] -\
     1024*f[-3*di + 3*dj + i] + 5376*f[-2*di + 3*dj + i] -\
     21504*f[-di + 3*dj + i] + 21504*f[di + 3*dj + i] -\
     5376*f[2*di + 3*dj + i] - 96*f[4*di + 3*dj + i] -\
     9*f[-4*di + 4*dj + i] + 96*f[-3*di + 4*dj + i] -\
     504*f[-2*di + 4*dj + i] + 2016*f[-di + 4*dj + i] -\
     2016*f[di + 4*dj + i] + 504*f[2*di + 4*dj + i] -\
     96*f[3*di + 4*dj + i] + 9*f[-4*(di + dj) + i] +\
     64*(16*f[-3*(di + dj) + i] + 441*(f[-2*(di + dj) + i] +\
     f[2*(di + dj) + i]) + 16*f[3*(di + dj) + i]) +\
     9*f[4*(di + dj) + i])*oo705600##dx##dy)
#define D8xx(f,i,di,dx) ((-14350*f[i] - 9*f[-4*di + i] +\
     16*(8*f[-3*di + i] - 63*(f[-2*di + i] - 8*(f[-di + i] +\
     f[di + i]) + f[2*di + i]) + 8*f[3*di + i]) -\
     9*f[4*di + i])*oo5040##dx##2)



/***********************************************************
** Kreiss-Oliger DISIPATION OPERATORS                     **
************************************************************/


/* First order dissipation */
#define DISSIP_any_1(f,i,di) (1.0/16.0 * (-2*f[i] + f[-di + i] +\
     f[di + i]))
#define DISS1(f) (DISSIP_any_1(f,ijk,di) +\
     DISSIP_any_1(f,ijk,dj) + DISSIP_any_1(f,ijk,dk) )

/* Third order dissipation */
#define DISSIP_any_3(f,i,di) (1.0/16.0 * (6*f[i] + f[-2*di + i] -\
     4*f[-di + i] - 4*f[di + i] + f[2*di + i]))
#define DISS3(f) (DISSIP_any_3(f,ijk,di) +\
     DISSIP_any_3(f,ijk,dj) + DISSIP_any_3(f,ijk,dk) )

/* Fifth order dissipation */
#define DISSIP_any_5(f,i,di) (1.0/16.0 * (-20*f[i] +\
     f[-3*di + i] - 6*f[-2*di + i] + 15*f[-di + i] +\
     15*f[di + i] - 6*f[2*di + i] + f[3*di + i]))
#define DISS5(f) (DISSIP_any_5(f,ijk,di) +\
     DISSIP_any_5(f,ijk,dj) + DISSIP_any_5(f,ijk,dk) )

/* Seventh order dissipation */
#define DISSIP_any_7(f,i,di) (1.0/16.0 * (70*f[i] + f[-4*di + i] -\
     8*f[-3*di + i] + 28*f[-2*di + i] - 56*f[-di + i] -\
     56*f[di + i] + 28*f[2*di + i] - 8*f[3*di + i] + f[4*di + i]))
#define DISS7(f) (DISSIP_any_7(f,ijk,di) +\
     DISSIP_any_7(f,ijk,dj) + DISSIP_any_7(f,ijk,dk) )

/* Ninth order dissipation */
#define DISSIP_any_9(f,i,di) (1.0/16.0 * (-252*f[i] +\
     f[-5*di + i] - 10*f[-4*di + i] + 45*f[-3*di + i] -\
     120*f[-2*di + i] + 210*f[-di + i] + 210*f[di + i] -\
     120*f[2*di + i] + 45*f[3*di + i] - 10*f[4*di + i] +\
     f[5*di + i]))
#define DISS9(f) (DISSIP_any_9(f,ijk,di) +\
     DISSIP_any_9(f,ijk,dj) + DISSIP_any_9(f,ijk,dk) )

#define _DERIVATIVES_H
#endif

