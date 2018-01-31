#!/usr/bin/env python
import numpy as np

import sixtracktools
import sixtracklib


six = sixtracktools.SixTrackInput('.')
line, rest, iconv = six.expand_struct()
names,types,args=zip(*line)
idx=dict( (nn,ii) for ii,nn in enumerate(six.struct) if not 'BLOC' in nn)
names2=np.array(names)[iconv]


sixtrackbeam=sixtracktools.SixDump3('dump3.dat')

block=sixtracklib.cBlock.from_line(line)
bref=sixtracklib.cBeam.from_full_beam(sixtrackbeam.get_full_beam())
bref=bref.reshape(-1,2)

import time
import sys

def mkbench(npart,nturn):
  block=sixtracklib.cBlock.from_line(line)
  cbeam=bref.copy().reshape(-1)[:npart]
  beg=time.time()
  #block.track(cbeam,nturn=nturn,turnbyturn=True)
  # desiable some diagnostics
  block.track(cbeam,nturn=nturn)
  end=time.time()
  perfcpu=((end - beg)/(npart*nturn))*1e6
  print("CPU part %4d, turn %4d: %10.3f usec/part*turn (total time: %10.3f s)"%(
    npart,nturn,perfcpu, end - beg))

  return end-beg,npart,nturn,perfcpu

out=open(time.strftime("bench_%Y%m%dT%H%m%S.txt"),'w')
#for npart in [100000,200000,500000,1000000]:
for npart in [1000, 2000, 5000, 10000]:
    #for nturn in [1,2,5,10]:
    for nturn in [10]:
        st,npart,nturn,perfcpu=mkbench(npart,nturn)
        fmt="%5d %5d %10.3f\n"
        out.write(fmt%(npart,nturn,perfcpu))

