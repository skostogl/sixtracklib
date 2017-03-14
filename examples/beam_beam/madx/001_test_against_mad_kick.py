import numpy as np
from scipy.constants import e as qe
from scipy.constants import m_p, c

intensity_pbun = 1e11
energy_GeV = 7000.
nemittx = 2.5e-6
nemitty = 1.e-6

tune_x = 0.
tune_y = 0.
beta_s = 0.40

include_beambeam = True
offsetx_s = 5e-5
offsety_s = -1e-5

#compute sigmas
mp_GeV = m_p*c**2/qe/1e9
gamma = energy_GeV/mp_GeV
sigmax_s = np.sqrt(beta_s*nemittx/gamma)
sigmay_s = np.sqrt(beta_s*nemitty/gamma)

sigma_avg = 0.5*(sigmax_s + sigmay_s)

theta_obs = np.pi/10
#~ theta_obs = np.pi/2.
#~ theta_obs = np.pi/3.


# set initial conditions
n_points = 100
r0_particles = np.array([0.] + list(np.linspace(-15*sigma_avg, 15*sigma_avg, n_points)))
x0_particles = r0_particles*np.cos(theta_obs)
px0_particles = 0*x0_particles
y0_particles = r0_particles*np.sin(theta_obs)
py0_particles = 0*x0_particles


    

import sys, os
BIN = os.path.expanduser("../../../")
sys.path.append(BIN)
import sixtracklib


### Build beam
beam=sixtracklib.cBeam(npart=len(x0_particles))
for ii in xrange(len(beam.particles)):
  beam.particles[ii]['partid'] = ii
  beam.particles[ii]['elemid'] = 0
  beam.particles[ii]['turn'] = 0
  beam.particles[ii]['state'] = 0
  beam.particles[ii]['s'] = 0
  beam.particles[ii]['x'] = x0_particles[ii]
  beam.particles[ii]['px'] = px0_particles[ii]
  beam.particles[ii]['y'] = y0_particles[ii]
  beam.particles[ii]['py'] = py0_particles[ii]
  beam.particles[ii]['sigma'] = 0.
  beam.particles[ii]['psigma'] = 0.
  beam.particles[ii]['chi'] = 1.
  beam.particles[ii]['delta'] = 0.
  beam.particles[ii]['rpp'] = 1.
  beam.particles[ii]['rvv'] = 1.
  beam.particles[ii]['beta'] = 1.
  beam.particles[ii]['gamma'] = gamma
  beam.particles[ii]['mass0'] = m_p*c**2/qe
  beam.particles[ii]['charge0'] = qe
  beam.particles[ii]['beta0'] = 1.
  beam.particles[ii]['gamma0'] = gamma
  beam.particles[ii]['p0c'] = energy_GeV*1e9
  
###  Build the ring
block=sixtracklib.cBlock(size=50)
#block.Multipole(bal=np.array([1.,2.,3.,4.,5.,6.]),l=0,hx=0,hy=0)
block.LinMap(alpha_x_s0=0., beta_x_s0=beta_s, D_x_s0=0., 
             alpha_x_s1=0., beta_x_s1=beta_s, D_x_s1=0.,
             alpha_y_s0=0., beta_y_s0=beta_s, D_y_s0=0.,
             alpha_y_s1=0., beta_y_s1=beta_s, D_y_s1=0.,
             dQ_x=tune_x, dQ_y=tune_y)
             
if include_beambeam:
  if np.abs((sigmax_s-sigmay_s)/((sigmax_s+sigmay_s)/2.))<1e-3:
    print "round beam"
    block.BB4D(N_s = intensity_pbun, beta_s = beta_s, q_s = qe, 
              transv_field_data = {'type':'gauss_round',  'sigma': (sigmax_s+sigmay_s)/2., 'Delta_x': offsetx_s, 'Delta_y': offsety_s})             
  else:
    print "elliptic beam"
    block.BB4D(N_s = intensity_pbun, beta_s = beta_s, q_s = qe, transv_field_data = {'type':'gauss_ellip',  'sigma_x': sigmax_s, 'sigma_y': sigmay_s, 'Delta_x': offsetx_s, 'Delta_y': offsety_s})             

block.Block() 


### Tracking stlb
track_fun =  block.track
# test OPENCL:
# track_fun =  block.track_cl
track_fun(beam)

#Remove Dipole kick
kick_x = beam.px[1:]-beam.px[0]
kick_y = beam.py[1:]-beam.py[0]



  
### Tracking MAD
import track_mad as tm
_, _, px_particles_mad, py_particles_mad = tm.track_mad_linmap_and_beambeam(intensity_pbun, energy_GeV, nemittx, nemitty, 
                    tune_x, tune_y, beta_s, include_beambeam, offsetx_s, offsety_s, sigmax_s, sigmay_s, 
                    x0_particles, px0_particles, y0_particles, py0_particles, nturns=1)
kick_x_mad = px_particles_mad[1:, 1]
kick_y_mad = py_particles_mad[1:, 1]
import pylab as pl
pl.close('all')
pl.figure(1)
ax1 = pl.subplot(2,1,1)
pl.plot(r0_particles[1:], kick_x_mad, 'b', label='mad')
pl.plot(r0_particles[1:], kick_x, '.r', label='sixtracklib')
pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='x') 
pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='y') 
pl.legend(loc='best')

pl.grid('on')
pl.subplot(2,1,2, sharex = ax1)
pl.plot(r0_particles[1:], kick_y_mad, 'b')
pl.plot(r0_particles[1:], kick_y, '.r')
pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='x') 
pl.gca().ticklabel_format(style='sci', scilimits=(0,0),axis='y') 

pl.suptitle('theta = %.1f deg'%(theta_obs*90/(np.pi/2)))
pl.grid('on')
pl.show()

wurst
