//SixTrackLib
//
//Authors: R. De Maria, G. Iadarola, D. Pellegrini, H. Jasim
//
//Copyright 2017 CERN. This software is distributed under the terms of the GNU
//Lesser General Public License version 2.1, copied verbatim in the file
//`COPYING''.
//
//In applying this licence, CERN does not waive the privileges and immunities
//granted to it by virtue of its status as an Intergovernmental Organization or
//submit itself to any jurisdiction.


#ifndef _PARTICLE_
#define _PARTICLE_

#include <stdint.h>


typedef struct Particles {
  // reference quantities
  double m0; // eV
  double q0; // C
  double q; // C
  double beta0;
  double gamma0;
  double p0c; //eV  
  
  //coordinate arrays
  uint64_t *partid;
  uint64_t *elemid;
  uint64_t *turn;
  uint64_t *state; //negativeparticle lost
  double *s;
  double *x;
  double *px; // Px/P0
  double *y;
  double *py; // Px/P0
  double *sigma; // s-beta0*c*t //t is the time since the beginning of the simulation
  double *psigma; // (E-E0)/ (beta0 P0c)
  double *chi; // q/q0 * m/m0
  double *delta; // P/P0-1 = 1/rpp-1
  double *rpp; // ratio P0/P
  double *rvv; // ratio beta / beta0
} Particles;
 

#endif
