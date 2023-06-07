#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>  // a useful debugging tool
#include <sys/time.h>  // needed for timing

#include "derivatives.h"
#include "looptest.h"
double get_time(void);

// These are used by the derivative macros.
#define oo12h2 1.0
#define oo5040h2 1.0
#define oo705600hh 1.0

#define di (1)
#define dj (nxyz)
#define dk (nxyz*nxyz)

/* A note about the derivative macro:
 * D4xx(fin, ijk, stride, h)
 * can be used to calculate x, y, or z second derivatives
 * The 'stride' determines the actual derivative. If stride=1, then
 * we get the xx derivative, if it's nxyz, then we get the yy derivative, 
 * and if it's nxyz*nxyz, then  we get the zz derivative. For
 * convenience, we define di, dj, dk. Note for these tests, we don't
 * actually care about dividing by the 12 * h * h, hence we set oo12h2
 * to 1
*/

/* std_loop_4th_order_deriv_xx
 * perform calculation associated with second "x" derivative
 * The arrays are constucted such that x_{i+1} is located right next
 * to x_i. This uses the standard looping order. The 3D data is packed
 * into a 1D array.
*/
double std_loop_4th_order_deriv_xx(const int gs, const int nxyz, double * restrict fin, double * restrict fout)
{
  const double tstart = get_time();
  for (int k = gs; k < nxyz - gs; ++k)
  {
    for (int j = gs; j < nxyz - gs; ++j)
    {
      for (int i = gs; i < nxyz - gs; ++i)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz; // the index of point (i,j,k)
        fout[ijk] =  D4xx(fin, ijk, di, h);
      }
    }
  }
  const double tend = get_time();
  return tend - tstart;
}

/* rev_loop_4th_order_deriv_xx
 * perform calculation associated with second "x" derivative
 * The arrays are constucted such that x_{i+1} is located right next
 * to x_i. This uses the reverse looping order
*/
double rev_loop_4th_order_deriv_xx(const int gs, const int nxyz, double * restrict fin, double * restrict fout)
{
  const double tstart = get_time();
  for (int i = gs; i < nxyz - gs; ++i)
  {
    for (int j = gs; j < nxyz - gs; ++j)
    {
      for (int k = gs; k < nxyz - gs; ++k)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz;
        fout[ijk] =  D4xx(fin, ijk, di, h);
      }
    }
  }
  const double tend = get_time();
  return tend - tstart;
}

/* std_loop_4th_order_deriv_yy
 * perform calculation associated with second "y" derivative
 * The arrays are constucted such that x_{i+1} is located right next
 * to x_i and y_{i+1} is located "nxyz" elements away from y_i.
 * This uses the standard looping order
*/
double std_loop_4th_order_deriv_yy(const int gs, const int nxyz, double * restrict fin, double * restrict fout)
{
  const double tstart = get_time();
  for (int k = gs; k < nxyz - gs; ++k)
  {
    for (int j = gs; j < nxyz - gs; ++j)
    {
      for (int i = gs; i < nxyz - gs; ++i)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz;
        fout[ijk] =  D4xx(fin, ijk, dj, h);
      }
    }
  }
  const double tend = get_time();
  return tend - tstart;
}

/* rev_loop_4th_order_deriv_yy
 * perform calculation associated with second "y" derivative
 * The arrays are constucted such that x_{i+1} is located right next
 * to x_i and y_{i+1} is located "nxyz" elements away from y_i.
 * This uses the reverse looping order
*/
double rev_loop_4th_order_deriv_yy(const int gs, const int nxyz, double * restrict fin, double * restrict fout)
{
  const double tstart = get_time();
  for (int i = gs; i < nxyz - gs; ++i)
  {
    for (int j = gs; j < nxyz - gs; ++j)
    {
      for (int k = gs; k < nxyz - gs; ++k)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz;
        fout[ijk] =  D4xx(fin, ijk, dj, h);
      }
    }
  }
  const double tend = get_time();
  return tend - tstart;
}

/* std_loop_4th_order_deriv_zz
 * perform calculation associated with second "z" derivative
 * The arrays are constucted such that x_{i+1} is located right next
 * to x_i and z_{i+1} is located "nxyz * nxyz" elements away from z_i.
 * This uses the standard looping order
*/
double std_loop_4th_order_deriv_zz(const int gs, const int nxyz, double * restrict fin, double * restrict fout)
{
  const double tstart = get_time();
  for (int k = gs; k < nxyz - gs; ++k)
  {
    for (int j = gs; j < nxyz - gs; ++j)
    {
      for (int i = gs; i < nxyz - gs; ++i)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz;
        fout[ijk] =  D4xx(fin, ijk, dk, h);
      }
    }
  }
  const double tend = get_time();
  return tend - tstart;
}

/* std_loop_4th_order_deriv_zz
 * perform calculation associated with second "z" derivative
 * The arrays are constucted such that x_{i+1} is located right next
 * to x_i and z_{i+1} is located "nxyz * nxyz" elements away from z_i.
 * This uses the reverse looping order
*/
double rev_loop_4th_order_deriv_zz(const int gs, const int nxyz, double * restrict fin, double * restrict fout)
{
  const double tstart = get_time();
  for (int i = gs; i < nxyz - gs; ++i)
  {
    for (int j = gs; j < nxyz - gs; ++j)
    {
      for (int k = gs; k < nxyz - gs; ++k)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz;
        fout[ijk] =  D4xx(fin, ijk, dk, h);
      }
    }
  }
  const double tend = get_time();
  return tend - tstart;
}

double std_loop_8th_order_deriv_xx(const int gs, const int nxyz, double * restrict fin, double * restrict fout)
{
  const double tstart = get_time();
  for (int k = gs; k < nxyz - gs; ++k)
  {
    for (int j = gs; j < nxyz - gs; ++j)
    {
      for (int i = gs; i < nxyz - gs; ++i)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz; // the index of point (i,j,k)
        fout[ijk] =  D8xx(fin, ijk, di, h);
      }
    }
  }
  const double tend = get_time();
  return tend - tstart;
}

double rev_loop_8th_order_deriv_xx(const int gs, const int nxyz, double * restrict fin, double * restrict fout)
{
  const double tstart = get_time();
  for (int i = gs; i < nxyz - gs; ++i)
  {
    for (int j = gs; j < nxyz - gs; ++j)
    {
      for (int k = gs; k < nxyz - gs; ++k)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz;
        fout[ijk] =  D8xx(fin, ijk, di, h);
      }
    }
  }
  const double tend = get_time();
  return tend - tstart;
}

double std_loop_8th_order_deriv_yy(const int gs, const int nxyz, double * restrict fin, double * restrict fout)
{
  const double tstart = get_time();
  for (int k = gs; k < nxyz - gs; ++k)
  {
    for (int j = gs; j < nxyz - gs; ++j)
    {
      for (int i = gs; i < nxyz - gs; ++i)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz;
        fout[ijk] =  D8xx(fin, ijk, dj, h);
      }
    }
  }
  const double tend = get_time();
  return tend - tstart;
}

double rev_loop_8th_order_deriv_yy(const int gs, const int nxyz, double * restrict fin, double * restrict fout)
{
  const double tstart = get_time();
  for (int i = gs; i < nxyz - gs; ++i)
  {
    for (int j = gs; j < nxyz - gs; ++j)
    {
      for (int k = gs; k < nxyz - gs; ++k)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz;
        fout[ijk] =  D8xx(fin, ijk, dj, h);
      }
    }
  }
  const double tend = get_time();
  return tend - tstart;
}

double std_loop_8th_order_deriv_zz(const int gs, const int nxyz, double * restrict fin, double * restrict fout)
{
  const double tstart = get_time();
  for (int k = gs; k < nxyz - gs; ++k)
  {
    for (int j = gs; j < nxyz - gs; ++j)
    {
      for (int i = gs; i < nxyz - gs; ++i)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz;
        fout[ijk] =  D8xx(fin, ijk, dk, h);
      }
    }
  }
  const double tend = get_time();
  return tend - tstart;
}

double rev_loop_8th_order_deriv_zz(const int gs, const int nxyz, double * restrict fin, double * restrict fout)
{
  const double tstart = get_time();
  for (int i = gs; i < nxyz - gs; ++i)
  {
    for (int j = gs; j < nxyz - gs; ++j)
    {
      for (int k = gs; k < nxyz - gs; ++k)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz;
        fout[ijk] =  D8xx(fin, ijk, dk, h);
      }
    }
  }
  const double tend = get_time();
  return tend - tstart;
}

double std_loop_8th_order_deriv_xy(const int gs, const int nxyz, double * restrict fin, double * restrict fout)
{
  const double tstart = get_time();
  for (int k = gs; k < nxyz - gs; ++k)
  {
    for (int j = gs; j < nxyz - gs; ++j)
    {
      for (int i = gs; i < nxyz - gs; ++i)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz;
        fout[ijk] = D8xy(fin, ijk, di, dj, h, h);
      }
    }
  }
  const double tend = get_time();
  return tend - tstart;
}

double rev_loop_8th_order_deriv_xy(const int gs, const int nxyz, double * restrict fin, double * restrict fout)
{
  const double tstart = get_time();
  for (int i = gs; i < nxyz - gs; ++i)
  {
    for (int j = gs; j < nxyz - gs; ++j)
    {
      for (int k = gs; k < nxyz - gs; ++k)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz;
        fout[ijk] =  D8xy(fin, ijk, di, dj, h, h);
      }
    }
  }
  const double tend = get_time();
  return tend - tstart;
}

double std_loop_8th_order_deriv_xz(const int gs, const int nxyz, double * restrict fin, double * restrict fout)
{
  const double tstart = get_time();
  for (int k = gs; k < nxyz - gs; ++k)
  {
    for (int j = gs; j < nxyz - gs; ++j)
    {
      for (int i = gs; i < nxyz - gs; ++i)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz;
        fout[ijk] = D8xy(fin, ijk, di, dk, h, h);
      }
    }
  }
  const double tend = get_time();
  return tend - tstart;
}

double rev_loop_8th_order_deriv_xz(const int gs, const int nxyz, double * restrict fin, double * restrict fout)
{
  const double tstart = get_time();
  for (int i = gs; i < nxyz - gs; ++i)
  {
    for (int j = gs; j < nxyz - gs; ++j)
    {
      for (int k = gs; k < nxyz - gs; ++k)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz;
        fout[ijk] =  D8xy(fin, ijk, di, dk, h, h);
      }
    }
  }
  const double tend = get_time();
  return tend - tstart;
}

double std_loop_8th_order_deriv_yz(const int gs, const int nxyz, double * restrict fin, double * restrict fout)
{
  const double tstart = get_time();
  for (int k = gs; k < nxyz - gs; ++k)
  {
    for (int j = gs; j < nxyz - gs; ++j)
    {
      for (int i = gs; i < nxyz - gs; ++i)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz;
        fout[ijk] = D8xy(fin, ijk, dj, dk, h, h);
      }
    }
  }
  const double tend = get_time();
  return tend - tstart;
}

double rev_loop_8th_order_deriv_yz(const int gs, const int nxyz, double * restrict fin, double * restrict fout)
{
  const double tstart = get_time();
  for (int i = gs; i < nxyz - gs; ++i)
  {
    for (int j = gs; j < nxyz - gs; ++j)
    {
      for (int k = gs; k < nxyz - gs; ++k)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz;
        fout[ijk] =  D8xy(fin, ijk, dj, dk, h, h);
      }
    }
  }
  const double tend = get_time();
  return tend - tstart;
}

void zero_out_arrays(const int nxyz, double * restrict fin, double * restrict fout)
{
  for (int k = 0; k < nxyz; ++k)
  {
    for (int j = 0; j < nxyz; ++j)
    {
      for (int i = 0; i < nxyz; ++i)
      {
        const size_t ijk = i + (j + k * nxyz) * nxyz;
        fout[ijk] =  0;
        fin[ijk] =  0;
      }
    }
  }
}

// get_time will return a double containing the current time in
// seconds.
double get_time(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);

  return (tv.tv_sec) + 1.0e-6 * tv.tv_usec;
}
