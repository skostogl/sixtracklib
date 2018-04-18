import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

sys.path.append('../postprocessing/NAFF_cpp')
sys.path.append('../NAFF_cpp')
from NAFF import *

sys.path.append('../postprocessing/NAFF_cpp/modules')
sys.path.append('../NAFF_cpp/modules')
from tune import *
from footprint import *


params = {'legend.fontsize': 12,
         'axes.labelsize':  12,
         'axes.titlesize':  12,
         'xtick.labelsize': 12,
         'ytick.labelsize': 12,
         'image.cmap':'jet',
         'lines.linewidth':1.,
         'lines.markersize': 5 }
plt.rcParams.update(params)



def beta_gamma(p0c_eV, pmass_eV):
  return p0c_eV/np.sqrt(p0c_eV**2+pmass_eV**2),np.sqrt(p0c_eV**2+pmass_eV**2)/pmass_eV

def sigmas(exn, eyn, p0c_eV, pmass_eV, betax, betay):
  beta0, gamma0 = beta_gamma(p0c_eV, pmass_eV)
  egeomx   = exn/(beta0*gamma0)
  egeomy   = eyn/(beta0*gamma0)
  sigma_x  = np.sqrt(egeomx*betax)
  sigma_px = np.sqrt(egeomx/betax)
  sigma_y  = np.sqrt(egeomy*betay)
  sigma_py = np.sqrt(egeomy/betay) 
  return sigma_x, sigma_px, sigma_y, sigma_py

#### Create an initial rectangular grid
def rect_grid(sigma_x,sigma_y, npart):
  x  = np.linspace(0.1*sigma_x, 6.1*sigma_x, np.sqrt(npart))
  y  = np.linspace(0.1*sigma_y, 6.1*sigma_y, np.sqrt(npart))
  xx, yy  = np.meshgrid(x,y)
  return xx.flatten(), yy.flatten()

#### Create an initial polar grid
def polar_grid(sigma_x_init, sigma_x_final, sigma_y_init, sigma_y_final,sigma_x,sigma_y, npart):
  r,t   = np.meshgrid(np.linspace(sigma_x_init, sigma_x_final, np.sqrt(npart)),np.linspace(0,  np.pi/2, np.sqrt(npart)))
  r2,t2   = np.meshgrid(np.linspace(sigma_y_init, sigma_y_final, np.sqrt(npart)),np.linspace(0,  np.pi/2, np.sqrt(npart)))
  X = r * np.cos(t)
  Y = r2 * np.sin(t2)
  return X.flatten(), Y.flatten()

#### Compute and plot fma
def plot_fma(xn,pxn,yn,pyn,x0,y0,title=None,half=500, xlim=[0.3,0.316], ylim=[0.305,0.323], save_to=None):
  qx_tot = []
  qy_tot = []
  qx_tot2 = []
  qy_tot2 = []
  npart = xn.shape[1]
  nturns = xn.shape[0]
  for i in range (npart):
    df = pd.DataFrame({'x': xn[0:half,i], 'px': pxn[0:half,i], 'y':yn[0:half,i], 'py': pyn[0:half,i]})
    df2 = pd.DataFrame({'x': xn[nturns-half:nturns,i], 'px': pxn[nturns-half:nturns,i], 'y':yn[nturns-half:nturns,i], 'py': pyn[nturns-half:nturns,i]})
    qx,amp= naff_cpp(df,[1], plane='x',fmax=1, write_file=False)
    qy,amp= naff_cpp(df,[1], plane='y',fmax=1, write_file=False)
    qx_tot.append(abs(qx[0]))
    qy_tot.append(abs(qy[0]))
    qx2,amp2= naff_cpp(df2,[1], plane='x',fmax=1, write_file=False)
    qy2,amp2= naff_cpp(df2,[1], plane='y',fmax=1, write_file=False)
    qx_tot2.append(abs(qx2[0]))
    qy_tot2.append(abs(qy2[0]))
  diff_x = [abs(x-y) for x,y in zip(qx_tot,qx_tot2)]
  diff_y = [abs(x-y) for x,y in zip(qy_tot,qy_tot2)]
  diff_x = [x**2 for x in diff_x]
  diff_y = [x**2 for x in diff_y]
  c1 = [np.log10(np.sqrt(x + y)) for x,y in zip(diff_x, diff_y) ]
  plt.xlabel(r'$ \rm Q_x$', fontsize=13)
  plt.ylabel(r'$ \rm Q_y$', fontsize=13)
  plt.xlim(xlim)
  plt.ylim(ylim)
  plot_res_upto_order(15,c1='darkgrey', c2='darkgrey', c3='r',annotate=False)
  plt.scatter(qx_tot, qy_tot, edgecolors=None,c=c1,s=2,cmap='jet', vmin=-7,vmax=-3)
  if title:
    plt.title(title)
  plt.tight_layout()
  if save_to:
    plt.savefig('fma/6d/newfma_%s_%s_%s.png' %(phi,sigmaz,N_slices))
    plt.close()
  else:
    plt.show()
  return pd.DataFrame({'qx_tot': qx_tot,"qy_tot": qy_tot,"c1": c1, "x0":x0, "y0":y0 })

'''
#### Plot inital phase space along with histogram, weight flag is to convert uniform distribution to gaussian
def plot_init_hist(x,y,weight, weights_flag=False):
  nullfmt = NullFormatter()
  left, width = 0.1, 0.65
  bottom, height = 0.1, 0.65
  bottom_h = left_h = left + width + 0.02
  rect_scatter = [left, bottom, width, height]
  rect_histx = [left, bottom_h, width, 0.2]
  rect_histy = [left_h, bottom, 0.2, height]
  plt.figure(1, figsize=(7, 7))
  axScatter = plt.axes(rect_scatter)
  axHistx = plt.axes(rect_histx)
  axHisty = plt.axes(rect_histy)
  axHistx.xaxis.set_major_formatter(nullfmt)
  axHisty.yaxis.set_major_formatter(nullfmt)
  axScatter.scatter(x,y, c=weight, cmap='jet', s=2, marker='o')
  axScatter.set_xlabel(r'$x \ [\sigma]$',fontsize=14)
  axScatter.set_ylabel(r'$y \ [\sigma]$',fontsize=14)
  binwidth = 0.25
  xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
  lim = (int(xymax/binwidth) + 1) * binwidth
  axScatter.set_xlim((0.0, 6.2))
  axScatter.set_ylim((0.0, 6.2))
  bins=100
  if weights_flag:
    axHistx.hist(x, bins=bins,edgecolor='k',color='r',weights=weight)
    axHisty.hist(y, bins=bins, orientation='horizontal',edgecolor='k',color='g',weights=weight)
  else:
    axHistx.hist(x, bins=bins,edgecolor='k',color='r')
    axHisty.hist(y, bins=bins, orientation='horizontal',edgecolor='k',color='g')
  axHistx.set_xlim(axScatter.get_xlim())
  axHisty.set_ylim(axScatter.get_ylim())
  axHistx.grid()
  axHisty.grid()
  plt.show()
'''

##### Save data in highest protocol
def save_data(df,result_folder):
    print "Saving tunes..."
    df.to_pickle(result_folder, protocol=-1)

##### Plots computed tunes and configuration space 
def plot_fma_from_file(qx,qy,c1,x0,y0,title=None,xlim=[0.3,0.316], ylim=[0.305,0.323], save_to=None, vmin=-7, vmax=-3, sigma_x =[0.1,6.1], sigma_y = [0.1,6.1]):
  fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12,6))
  if title:
    plt.suptitle(title)
  plt.sca(axes[0])
  plt.xlim(xlim)
  plt.ylim(ylim)
  plot_res_upto_order(15,c1 = 'darkgrey', c2='darkgrey')
  plt.scatter(qx, qy, edgecolors=None,c=c1,s=2,cmap='jet', vmin=vmin,vmax=vmax)
  plt.xlabel(r"$ \rm Q_x$")
  plt.ylabel(r"$\rm Q_y$")
  plt.sca(axes[1])
  plt.scatter(x0, y0, edgecolors=None,c=c1,s=2,cmap='jet', vmin=vmin,vmax=vmax)
  plt.xlabel(r"$ \rm \sigma_x$")
  plt.ylabel(r"$\rm \sigma_y$")
  plt.xlim(sigma_x)
  plt.ylim(sigma_y)
  plt.grid()
  cbar = plt.colorbar()
  cbar.set_label(r'$ \rm log(\sqrt{(Q_{x2}-Q_{x1})^2 + (Q_{y2}-Q_{y1})^2})$', rotation=270,labelpad=30)
  plt.subplots_adjust(left=0.09, bottom=0.12, right=0.95, top=0.9, wspace=None, hspace=None)
  if save_to:
    plt.savefig(save_to)
    plt.close()
  else:
    plt.show()

