intensity_pbun = 1e11
energy_GeV = 7000.
nemittx = 2e-6
nemitty = 2.5e-6

tune_x = 0.31
tune_y = 0.32
beta_s = 0.40

x_part = 1.
px_part = 0.
y_part = 2.
py_part = 0.
nturns  = 1024

include_beambeam = True



with open('mad_auto.madx', 'w') as fmad:
  fmad.write("beam, particle=proton, npart = %.2fe11, energy=%.2f, exn=%e, eyn=%e;\n\n"%(intensity_pbun/1e11, energy_GeV, nemittx, nemitty))
  fmad.write("tune_x = %.4f;\ntune_y =  %.4f;\nbeta_s = %.4f;\n"%(tune_x, tune_y, beta_s))
  
  fmad.write("""
one_turn: matrix, 
  rm11= cos(2*pi*tune_x),        rm12=sin(2*pi*tune_x)*beta_s,
  rm21=-sin(2*pi*tune_x)/beta_s, rm22=cos(2*pi*tune_x),
  rm33= cos(2*pi*tune_y),        rm34=sin(2*pi*tune_y)*beta_s,
  rm43=-sin(2*pi*tune_y)/beta_s, rm44=cos(2*pi*tune_y)  

;

linmap:   line=(one_turn);
""")

 #Insert beam beam here
  fmad.write("""
myring: line=(linmap);
use,period=myring;


""")
  
  #Track
  fmad.write("track, dump;\n")
  fmad.write("start, x= %e, px=%e, y=%e, py=%e;\n"%(x_part, px_part, y_part, py_part))
  
  fmad.write("run,turns=%d;\nendtrack;\n"%nturns)


import metaclass 

ob = metaclass.twiss('track.obs0001.p0001')

import harmonic_analysis as ha
extract_tune = lambda signal: ha.HarmonicAnalysis(signal).laskar_method(num_harmonics=1)[0][0]

tune  = extract_tune(ob.X)
