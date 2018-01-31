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


#ifndef _GPUCODE

#include <math.h>
#include <stdio.h>

#ifndef __APPLE__
#define M_PI 3.14159265358979323846
#endif

#endif


#define CLIGHT 299792458

#include "particle.h"

#include "track.h"

static inline int Drift_track(Particles* p, uint64_t ip, CLGLOBAL Drift *el){
  double xp, yp;
  double length=el->length;
  xp = p->px[ip] * p->rpp[ip];
  yp = p->py[ip] * p->rpp[ip];
  p->x[ip] += xp * length;
  p->y[ip] += yp * length;
  p->sigma[ip] += length * (1 - p->rvv[ip]*( 1 + (xp*xp+yp*yp)/2 ) );
  p->s[ip]+=length;
//  _DP("Drift_track: length=%g\n",length);
  return 1;
};


static inline int DriftExact_track(Particles* p, uint64_t ip, CLGLOBAL DriftExact *el){
  double lpzi, lbzi, px, py, opd;
  double length = el->length;
  opd=1+p->delta[ip];
  px=p->px[ip]; py=p->py[ip];
  lpzi= length/sqrt(opd*opd-px*px-py*py);
  lbzi=(p->beta0[ip]*p->beta0[ip]*p->psigma[ip]+1)*lpzi;
  p->x[ip] += px*lpzi ;
  p->y[ip] += py*lpzi ;
  p->sigma[ip] += length - lbzi;
  p->s[ip] += length ;
  return 1;
}

static inline int Multipole_track(Particles* p, uint64_t ip, CLGLOBAL Multipole *el){
  double x,y,chi,dpx,dpy,zre,zim,b1l,a1l,hxx,hyy;
  long int order=el->order;
  double hxl=el->hxl;
  double hyl=el->hyl;
  double l=el->l;
  CLGLOBAL double *bal = el->bal;
  dpx=bal[order*2];
  dpy=bal[order*2+1];
  x=p->x[ip]; y=p->y[ip]; chi=p->chi[ip];
//  _DP("Multipole_track: dpx,y=%g %G\n",dpx,dpy);
  for (int ii=order-1;ii>=0;ii--){
    zre=(dpx*x-dpy*y);
    zim=(dpx*y+dpy*x);
//    _DP("Multipole_track: y,x=%g %G\n",x,y);
    dpx=bal[ii*2]+zre;
    dpy=bal[ii*2+1]+zim;
//    _DP("Multipole_track: dpx,y=%g %G\n",dpx,dpy);
  }
  dpx=-chi*dpx ;
  dpy=chi*dpy ;
//  _DP("Multipole_track: dpx,y=%g %G\n",dpx,dpy);
  if (l>0){
     b1l=chi*bal[0]; a1l=chi*bal[1];
     hxx=hxl/l*x; hyy=hyl/l*y;
     dpx+=hxl + hxl*p->delta[ip] - b1l*hxx;
     dpy-=hyl + hyl*p->delta[ip] - a1l*hyy;
     p->sigma[ip]-=chi*(hxx-hyy)*l*p->rvv[ip];
  }
  p->px[ip]+=dpx ;  p->py[ip]+=dpy ;
  return 1 ;
}

static inline int Cavity_track(Particles* p, uint64_t ip, CLGLOBAL Cavity *el){
  double volt = el->volt;
  double freq = el->freq;
  double lag = el->lag;
  double phase, pt, opd;
  phase=lag-2*M_PI/CLIGHT*freq*p->sigma[ip]/p->beta0[ip];
  //printf("ggg00 %e %e\n",p->psigma,p->psigma+p->chi*volt/(p->p0c));
  p->psigma[ip]+=p->chi[ip]*volt*sin(phase)/(p->p0c[ip]*p->beta0[ip]);
  pt=p->psigma[ip] * p->beta0[ip];
  opd=sqrt( pt*pt+ 2*p->psigma[ip] + 1 );
  p->delta[ip]=opd - 1;
  p->beta[ip]=opd/(1/p->beta0[ip]+pt);
  //p->gamma=1/sqrt(1-p->beta*p->beta);
  p->gamma[ip]=(pt*p->beta0[ip]+1)*p->gamma0[ip];
  p->rpp[ip]=1/opd;
  p->rvv[ip]=p->beta0[ip]/p->beta[ip];
  //printf("ggg2 %e %e %e\n",pt,opd,p->delta);
  return 1;
}

static inline int Align_track(Particles* p, uint64_t ip, CLGLOBAL Align *el){
  double xn,yn;
  double cz = el->cz;
  double sz = el->sz;
  double dx = el->dx;
  double dy = el->dy;
  xn= cz*p->x[ip]-sz*p->y[ip] - dx;
  yn= sz*p->x[ip]+cz*p->y[ip] - dy;
  p->x[ip]=xn;
  p->y[ip]=yn;
  xn= cz*p->px[ip]+sz*p->py[ip];
  yn=-sz*p->px[ip]+cz*p->py[ip];
  p->px[ip]=xn;
  p->py[ip]=yn;
  return 1;
};
