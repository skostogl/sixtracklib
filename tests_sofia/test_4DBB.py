import numpy as np
import pandas as pd
from scipy import special
import matplotlib.pyplot as plt
from scipy.constants import e as qe

import sys; sys.path.append('../')

import sixtracklib
from utils import *


def test_track():

    fRF      = 400e6
    V_RF     = 16e6
    lag_deg  = 180
    p0c_eV   = 6.5e12  
    qx       = 0.31
    qy       = 0.32
    betax    = 122.21
    betay    = 210.1048
    h        = 35640.
    alfax    = 0.0
    alfay    = 0.
    gamma_tr = 55.68
    pmass_eV = 938.272046e6
    exn      = 2.5e-6  
    eyn      = 2.5e-6  
    npart    = 100**2
    nturns   = 2000
    
    #ap = 1./(gamma_tr**2)
    ap = 0.0
    beta0, gamma0 = beta_gamma(p0c_eV, pmass_eV)
    sigma_x, sigma_px, sigma_y, sigma_py = sigmas(exn, eyn, p0c_eV, pmass_eV, betax, betay)

    machine=sixtracklib.CBlock()
    machine.add_LinearMap(qx=qx,qy=qy,betax=betax,betay=betay,alfax=alfax,alfay=alfay,ap=ap,h=h, fRF=fRF)

    q_part = qe
    N_part = 2.e11
    beta_s = 0.40
    min_sigma_diff = 1e-16

    machine.add_BeamBeam4D(name='bb4d',
    q_part = q_part,
    N_part = N_part,
    sigma_x = sigma_x,
    sigma_y = sigma_y,
    beta_s = beta_s,
    min_sigma_diff = min_sigma_diff)
    
   
    machine.add_Multipole(knl=[0.0,0.0,0.0,500.])
    machine.add_Multipole(knl=[0.0,0.0,0.2])

    #machine.add_Cavity(voltage=V_RF,frequencx=fRF,lag=lag_deg)

    bunch=sixtracklib.CParticles(npart=npart,
                  p0c=p0c_eV,
                  beta0 = beta0,
                  gamma0 = gamma0)
    #X,Y =rect_grid(sigma_x,sigma_y, npart)
    X,Y =polar_grid(0.01*sigma_x, 6.*sigma_x, 0.01*sigma_y, 6.*sigma_y, sigma_x,sigma_y, npart)
    bunch.x = X
    bunch.y = Y

    x0 = bunch.x/sigma_x
    y0 = bunch.y/sigma_y
    bunch.set_delta(27e-5)

    '''
    weightx = 1.0-special.erf(x0/(np.sqrt(2)))
    weighty = 1.0-special.erf(y0/(np.sqrt(2)))
    weight = weightx*weighty
    weight = weight/max(weight)
    plot_init_hist(x0,y0,weight, weights_flag=True)
    '''
    
    particles,ebe,tbt=machine.track_cl(bunch,nturns=nturns,
                                    elembyelem=None,turnbyturn=True)
    x,px,y,py = tbt.x, tbt.px,tbt.y,tbt.py
    df_fma = plot_fma(x,px,y,py,x0,y0,half=1000)

    result_folder = '../postprocessing/data/test_data'
    save_data(df_fma,result_folder)

    return machine,particles,ebe,tbt


if __name__=='__main__':
    machine,particles,ebe,tbt=test_track()

