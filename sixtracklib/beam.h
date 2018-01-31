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



#ifndef _BEAM_
#define _BEAM_
#include "value.h"

#include "particle.h"

#ifdef VECTOR

typedef struct Beam {
  uint64_t npart;
  Particles* particles;
} Beam;

#else

typedef struct Beam {
  uint64_t npart;
  Particle* particles;
} Beam;

#endif

#ifndef _GPUCODE

Beam *Beam_new(uint64_t npart);

#endif

#endif
