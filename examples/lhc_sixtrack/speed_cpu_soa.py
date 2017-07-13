#!/usr/bin/env python

import time
import sixtracktools
import sixtracklib


def mkbench(npart, nturn, block, bref_soa):
    cbeam = bref_soa[:npart]
    beg = time.time()
    # block.track(cbeam,nturn=nturn,turnbyturn=True)
    # desiable some diagnostics
    block.track(cbeam, nturn=nturn, vec=True)
    end = time.time()
    perfcpu = ((end - beg)/(npart*nturn))*1e6
    print("CPU part %4d, turn %4d: %10.3f usec/part*turn (total time: %10.3f s)"
          % (npart, nturn, perfcpu, end - beg))
    return end-beg, npart, nturn, perfcpu


def main():
    six = sixtracktools.SixTrackInput('.')
    line, _, _ = six.expand_struct()
    sixtrackbeam = sixtracktools.SixDump3('dump3.dat')
    block = sixtracklib.cBlock.from_line(line)
    bref_soa = sixtracklib.cBeam_SoA.from_full_beam(sixtrackbeam.get_full_beam())
    #  bref = bref.reshape(-1, 2)

    out = open(time.strftime("bench_%Y%m%dT%H%m%S.txt"), 'w')
    #  for npart in [100000,200000,500000,1000000]:
    for npart in [1000, 2000, 5000, 10000]:
        # for nturn in [1,2,5,10]:
        for nturn in [10]:
            _, npart, nturn, perfcpu = mkbench(npart, nturn, block, bref_soa)
            fmt = "%5d %5d %10.3f\n"
            out.write(fmt % (npart, nturn, perfcpu))

if __name__ == "__main__":
    main()
