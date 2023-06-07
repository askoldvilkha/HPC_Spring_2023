#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "looptest.h"

#define gs 4 // keep this at 4
int main(int argc, char **argv)
{
  if (argc != 2)
  {
    fprintf(stderr, "Must supply exacly one argument, the number of points per dimension in the computational domain\n");
    exit(-1);
  }

  const int nxyz = atoi(argv[1]);
  const size_t one_d_s = nxyz;

  const size_t three_d_size = one_d_s * one_d_s * one_d_s;  // needs to be size_t, could overflow otherwise

  double *fin = calloc(three_d_size, sizeof(double));
  double *fout = calloc(three_d_size, sizeof(double));

  assert(fin);
  assert(fout);


  /************************************************************************
   * Choose the type of loop
   * Uncomment only one of the following
   */

  //const double time = std_loop_4th_order_deriv_xx(gs, nxyz , fin, fout);
  //const double time = rev_loop_4th_order_deriv_xx(gs, nxyz , fin, fout);

  //const double time = std_loop_4th_order_deriv_yy(gs, nxyz , fin, fout);
  //const double time = rev_loop_4th_order_deriv_yy(gs, nxyz , fin, fout);

  //const double time = std_loop_4th_order_deriv_zz(gs, nxyz , fin, fout);
  //const double time = rev_loop_4th_order_deriv_zz(gs, nxyz , fin, fout);

  //const double time = std_loop_8th_order_deriv_xx(gs, nxyz , fin, fout);
  //const double time = rev_loop_8th_order_deriv_xx(gs, nxyz , fin, fout);

  //const double time = std_loop_8th_order_deriv_yy(gs, nxyz , fin, fout);
  //const double time = rev_loop_8th_order_deriv_yy(gs, nxyz , fin, fout);

  //const double time = std_loop_8th_order_deriv_zz(gs, nxyz , fin, fout);
  //const double time = rev_loop_8th_order_deriv_zz(gs, nxyz , fin, fout);

  const double time = std_loop_8th_order_deriv_xy(gs, nxyz , fin, fout);
  //const double time = rev_loop_8th_order_deriv_xy(gs, nxyz , fin, fout);

  //const double time = std_loop_8th_order_deriv_xz(gs, nxyz , fin, fout);
  //const double time = rev_loop_8th_order_deriv_xz(gs, nxyz , fin, fout);

  //const double time = std_loop_8th_order_deriv_yz(gs, nxyz , fin, fout);
  //const double time = rev_loop_8th_order_deriv_yz(gs, nxyz , fin, fout);

  for (int i=0; i<10; i++)
  {
    printf ("%d %e\n", nxyz, time);  
  }
  

  return 0;
}
