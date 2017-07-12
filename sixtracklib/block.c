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


#include <stdio.h>
#include "block.h"

#ifdef VECTOR
#include "track_vec.c"
#else
#include "track.c"
#endif


// Data management

static inline type_t get_type(CLGLOBAL value_t *data, uint64_t elemid ) {
  return (type_t) data[elemid].i64;
}

//Block

static inline uint64_t Block_get_nelen(CLGLOBAL value_t *data, size_t elemid ) {
  return data[elemid + 1].i64;
}

static inline CLGLOBAL uint64_t *Block_get_elemids(CLGLOBAL value_t *data, size_t elemid ) {
  return &data[elemid + 2].u64 ;
}

// Tracking single

#ifdef VECTOR
int track_single(CLGLOBAL value_t *data,
                 CLGLOBAL Particles *particles,
                 uint64_t elemid,
                 uint64_t npart,
                 uint64_t elembyelemoff, uint64_t turnbyturnoff) {

   enum type_t typeid = get_type(data, elemid);
   CLGLOBAL value_t *elem = data+elemid+1; //Data starts after typeid
   // diagnostics, will be always 0 if diagnostics disabled in Block_track
   // TODO port to new SOA data structure
   /*if (turnbyturnoff>0) {
     for (uint64_t i_part=0; i_part < npart; i_part++) {
     Particle *p = &particles[i_part];
       uint64_t dataoff=turnbyturnoff+sizeof(Particle)/8 * i_part;
       for (int i_attr=0;i_attr<sizeof(Particle)/8;i_attr++) {
         data[dataoff + i_attr] = ((value_t *) p)[i_attr];
       }
     }
   }*/
   // for each particle ...
//       _DP("Block_track: elemid=%zu typedid=%u\n",elemid,typeid);
     switch (typeid) {
       // 50%
         case DriftID:
              for (uint64_t i_part=0; i_part < npart; i_part++) {
                if (particles->state[i_part] > 0 ) {
                  Drift_track(particles, i_part, (CLGLOBAL Drift*) elem);
                }
              }
         break;
         // 30%
         case MultipoleID:
              for (uint64_t i_part=0; i_part < npart; i_part++) {
                if (particles->state[i_part] > 0 ) {
                  Multipole_track(particles, i_part, (CLGLOBAL Multipole*) elem);
                }
              }
         break;
         // 0.1%
         case CavityID:
              for (uint64_t i_part=0; i_part < npart; i_part++) {
                if (particles->state[i_part] > 0 ) {
                  Cavity_track(particles, i_part, (CLGLOBAL Cavity*) elem);
                }
              }
         break;
         // 0 - 50%
         case AlignID:
              for (uint64_t i_part=0; i_part < npart; i_part++) {
                if (particles->state[i_part] > 0 ) {
                  Align_track(particles, i_part, (CLGLOBAL Align*) elem);
                }
              }
         break;
         case IntegerID: break;
         case DoubleID: break;
         case BlockID: break;
         case DriftExactID:
              for (uint64_t i_part=0; i_part < npart; i_part++) {
                if (particles->state[i_part] > 0 ) {
                  DriftExact_track(particles, i_part, (CLGLOBAL DriftExact*) elem);
                }
              }
         break;
     }
     // diagnostics, will be always 0 if diagnostics disabled in Block_track
     // TODO port to new SOA data structure
     /*if (elembyelemoff>0){
       for (uint64_t i_part=0; i_part < npart; i_part++) {
         Particle *p = &particles[i_part];
         uint64_t dataoff=elembyelemoff+sizeof(Particle)/8 * i_part;
         for (int i_attr=0;i_attr<sizeof(Particle)/8;i_attr++) {
            data[dataoff + i_attr] = ((value_t *) p)[i_attr];
         }
       }
    }*/
   return 1;
}
#else

int track_single(CLGLOBAL value_t *data,
                 CLGLOBAL Particle *particles,
                 uint64_t elemid,
                 uint64_t npart,
                 uint64_t elembyelemoff, uint64_t turnbyturnoff) {

   enum type_t typeid = get_type(data, elemid);
   CLGLOBAL value_t *elem = data+elemid+1; //Data starts after typeid
   // diagnostics, will be always 0 if diagnostics disabled in Block_track
   if (turnbyturnoff>0) {
     for (uint64_t i_part=0; i_part < npart; i_part++) {
     Particle *p = &particles[i_part];
       uint64_t dataoff=turnbyturnoff+sizeof(Particle)/8 * i_part;
       for (int i_attr=0;i_attr<sizeof(Particle)/8;i_attr++) {
         data[dataoff + i_attr] = ((value_t *) p)[i_attr];
       }
     }
   }
   // for each particle ...
//       _DP("Block_track: elemid=%zu typedid=%u\n",elemid,typeid);
     switch (typeid) {
       // 50%
         case DriftID:
              for (uint64_t i_part=0; i_part < npart; i_part++) {
                Particle *p = &particles[i_part];
                if (p->state > 0 ) {
                  Drift_track(p, (CLGLOBAL Drift*) elem);
                }
              }
         break;
         // 30%
         case MultipoleID:
              for (uint64_t i_part=0; i_part < npart; i_part++) {
                Particle *p = &particles[i_part];
                if (p->state > 0 ) {
                  Multipole_track(p, (CLGLOBAL Multipole*) elem);
                }
              }
         break;
         // 0.1%
         case CavityID:
              for (uint64_t i_part=0; i_part < npart; i_part++) {
                Particle *p = &particles[i_part];
                if (p->state > 0 ) {
                  Cavity_track(p, (CLGLOBAL Cavity*) elem);
                }
              }
         break;
         // 0 - 50%
         case AlignID:
              for (uint64_t i_part=0; i_part < npart; i_part++) {
                Particle *p = &particles[i_part];
                if (p->state > 0 ) {
                  Align_track(p, (CLGLOBAL Align*) elem);
                }
              }
         break;
         case IntegerID: break;
         case DoubleID: break;
         case BlockID: break;
         case DriftExactID:
              for (uint64_t i_part=0; i_part < npart; i_part++) {
                Particle *p = &particles[i_part];
                if (p->state > 0 ) {
                  DriftExact_track(p, (CLGLOBAL DriftExact*) elem);
                }
              }
         break;
     }
     // diagnostics, will be always 0 if diagnostics disabled in Block_track
     if (elembyelemoff>0){
       for (uint64_t i_part=0; i_part < npart; i_part++) {
         Particle *p = &particles[i_part];
         uint64_t dataoff=elembyelemoff+sizeof(Particle)/8 * i_part;
         for (int i_attr=0;i_attr<sizeof(Particle)/8;i_attr++) {
            data[dataoff + i_attr] = ((value_t *) p)[i_attr];
         }
       }
    }
   return 1;
}

#endif

// Tracking loop

#ifdef _GPUCODE

CLKERNEL void Block_track(
                 CLGLOBAL value_t *data, CLGLOBAL Particle *particles,
                 uint64_t blockid, uint64_t nturn, uint64_t npart,
                 uint64_t elembyelemid, uint64_t turnbyturnid){
   uint64_t nelem    = Block_get_nelen(data, blockid);
   CLGLOBAL uint64_t *elemids = Block_get_elemids(data, blockid);
   uint64_t i_part = get_global_id(0);
   uint64_t elembyelemoff=0;
   uint64_t turnbyturnoff=0;
   Particle pp=particles[i_part];
   for (int i_turn=0; i_turn< nturn; i_turn++){
     for (int i_elem=0; i_elem< nelem; i_elem++) {
       if (elembyelemid>0){
         elembyelemoff=elembyelemid +
                      sizeof(Particle)/8 * npart * i_turn +
                      sizeof(Particle)/8 * npart * nturn  * i_elem ;
//            printf("%lu \n",elembyelemoff);
       }
       if (turnbyturnid>0){
         turnbyturnoff=turnbyturnid +
                        sizeof(Particle)/8 * npart * i_turn;
       }
      track_single(data, particles, elemids,
                   &pp, i_elem, i_part,elembyelemoff, turnbyturnoff);
    }
    if (particles[i_part].state>=0) {
      particles[i_part].turn++;
    }
  }
  particles[i_part]=pp;
}

#else


#ifdef VECTOR

int Block_track(value_t *data, Beam* restrict beam,
                uint64_t blockid, uint64_t nturn,
                uint64_t elembyelemid, uint64_t turnbyturnid){
   uint64_t nelem    = Block_get_nelen(data, blockid);
   uint64_t *elemids = Block_get_elemids(data, blockid);
   uint64_t npart=beam->npart;
   //uint64_t elembyelemoff=0; XXX
   uint64_t turnbyturnoff=0;
   // for each revolution around accelerator ...
   for (int i_turn=0; i_turn< nturn; i_turn++) {
     // for each accelerator element ...
     for (int i_elem=0; i_elem< nelem; i_elem++) {
       uint64_t elemid=elemids[i_elem];
          // diagnostics
          // TODO port diagnostics code
          /*
          if (elembyelemid>0){
            elembyelemoff=elembyelemid +
                         sizeof(Particle)/8 * npart * i_turn +
                         sizeof(Particle)/8 * npart * nturn  * i_elem ;
//            printf("cpu %lu \n",elembyelemoff);
          }
          // diagnostics
          if (turnbyturnid>0){
            turnbyturnoff=turnbyturnid +
                         sizeof(Particle)/8 * npart * i_turn;
//            printf("%lu \n",turnbyturnoff);
          }*/
          track_single(data, beam->particles, elemid, npart,
              elembyelemid, turnbyturnoff);
       }
     for (uint64_t i_part=0; i_part < npart; i_part++){
       if (beam->particles[i_part].state >= 0) beam->particles[i_part].turn++;
     }
   }
   return 1;
}

#else

int Block_track(value_t *data, Beam *restrict beam,
                uint64_t blockid, uint64_t nturn,
                uint64_t elembyelemid, uint64_t turnbyturnid){
   uint64_t nelem    = Block_get_nelen(data, blockid);
   uint64_t *elemids = Block_get_elemids(data, blockid);
   uint64_t npart=beam->npart;
   uint64_t elembyelemoff=0;
   uint64_t turnbyturnoff=0;
   // for each revolution around accelerator ...
   for (int i_turn=0; i_turn< nturn; i_turn++) {
     // for each accelerator element ...
     for (int i_elem=0; i_elem< nelem; i_elem++) {
       uint64_t elemid=elemids[i_elem];
          // diagnostics
          if (elembyelemid>0){
            elembyelemoff=elembyelemid +
                         sizeof(Particle)/8 * npart * i_turn +
                         sizeof(Particle)/8 * npart * nturn  * i_elem ;
//            printf("cpu %lu \n",elembyelemoff);
          }
          // diagnostics
          if (turnbyturnid>0){
            turnbyturnoff=turnbyturnid +
                         sizeof(Particle)/8 * npart * i_turn;
//            printf("%lu \n",turnbyturnoff);
          }
          track_single(data, beam->particles, elemid, npart,
              elembyelemid, turnbyturnoff);
       }
     for (uint64_t i_part=0; i_part < npart; i_part++){
       if (beam->particles[i_part].state >= 0) beam->particles[i_part].turn++;
     }
   }
   return 1;
}
#endif // !VECTOR

#endif // !GPUCODE
