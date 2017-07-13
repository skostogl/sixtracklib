import ctypes
import numpy as np

particle_t=np.dtype([
         ('partid' ,'int32'),
         ('elemid' ,'int32'),
         ('turn'   ,'int32'),
         ('state'  ,'int32'),
         ('s'      ,float),
         ('x'      ,float),
         ('px'     ,float),
         ('y'      ,float),
         ('py'     ,float),
         ('sigma'  ,float),
         ('psigma' ,float),
         ('chi'    ,float),
         ('delta'  ,float),
         ('rpp'    ,float),
         ('rvv'    ,float),
         ('beta'   ,float),
         ('gamma'  ,float),
         ('mass0'  ,float),
         ('charge0',float),
         ('charge' ,float),
         ('beta0'  ,float),
         ('gamma0' ,float),
         ('p0c'    ,float),
         ])


class Particles(object):

    def __init__(self, npart, pdict=None):
        self.np = npart
        if pdict is None:
            self.partid = np.zeros(npart, dtype=np.int32, order='C')
            self.elemid = np.zeros(npart, dtype=np.int32, order='C')
            self.turn = np.zeros(npart, dtype=np.int32, order='C')
            self.state = np.zeros(npart, dtype=np.int32, order='C')
            self.s = np.zeros(npart, dtype=np.float64, order='C')
            self.x = np.zeros(npart, dtype=np.float64, order='C')
            self.y = np.zeros(npart, dtype=np.float64, order='C')
            self.px = np.zeros(npart, dtype=np.float64, order='C')
            self.py = np.zeros(npart, dtype=np.float64, order='C')
            self.sigma = np.zeros(npart, dtype=np.float64, order='C')
            self.psigma = np.zeros(npart, dtype=np.float64, order='C')
            self.chi = np.zeros(npart, dtype=np.float64, order='C')
            self.delta = np.zeros(npart, dtype=np.float64, order='C')
            self.rpp = np.zeros(npart, dtype=np.float64, order='C')
            self.rvv = np.zeros(npart, dtype=np.float64, order='C')
            self.beta = np.zeros(npart, dtype=np.float64, order='C')
            self.beta0 = np.zeros(npart, dtype=np.float64, order='C')
            self.gamma = np.zeros(npart, dtype=np.float64, order='C')
            self.gamma0 = np.zeros(npart, dtype=np.float64, order='C')
            self.mass0 = np.zeros(npart, dtype=np.float64, order='C')
            self.charge = np.zeros(npart, dtype=np.float64, order='C')
            self.charge0 = np.zeros(npart, dtype=np.float64, order='C')
            self.p0c = np.zeros(npart, dtype=np.float64, order='C')
        else:
            for k, v in pdict.iteritems():
                self.__setattr__(k, v)

    def __getitem__(self, key):
        new = {}
        for nn, prop, in self.__dict__.iteritems():
            if nn == 'np': continue  # not an array
            new[nn] = prop.copy()[key]
        newp = Particles(len(new['partid']), pdict=new)
        return newp

    def copy(self):
        new = Particles(self.np)
        for nn, prop, in self.__dict__.iteritems():
            if nn == 'np': continue  # not an array
            new.__setattr__(nn, prop.copy())
        return new


class Particles_ctypes(ctypes.Structure):
    _fields_ = [('partid',  ctypes.POINTER(ctypes.c_int)),
                ('elemid',  ctypes.POINTER(ctypes.c_int)),
                ('turn',    ctypes.POINTER(ctypes.c_int)),
                ('state',   ctypes.POINTER(ctypes.c_int)),
                ('s',       ctypes.POINTER(ctypes.c_double)),
                ('x',       ctypes.POINTER(ctypes.c_double)),
                ('px',      ctypes.POINTER(ctypes.c_double)),
                ('y',       ctypes.POINTER(ctypes.c_double)),
                ('py',      ctypes.POINTER(ctypes.c_double)),
                ('sigma',   ctypes.POINTER(ctypes.c_double)),
                ('psigma',  ctypes.POINTER(ctypes.c_double)),
                ('chi',     ctypes.POINTER(ctypes.c_double)),
                ('delta',   ctypes.POINTER(ctypes.c_double)),
                ('rpp',     ctypes.POINTER(ctypes.c_double)),
                ('rvv',     ctypes.POINTER(ctypes.c_double)),
                ('beta',    ctypes.POINTER(ctypes.c_double)),
                ('gamma',   ctypes.POINTER(ctypes.c_double)),
                ('mass0',   ctypes.POINTER(ctypes.c_double)),
                ('charge0', ctypes.POINTER(ctypes.c_double)),
                ('charge',  ctypes.POINTER(ctypes.c_double)),
                ('beta0',   ctypes.POINTER(ctypes.c_double)),
                ('gamma0',  ctypes.POINTER(ctypes.c_double)),
                ('p0c',     ctypes.POINTER(ctypes.c_double))]


class cBeam_ctypes(ctypes.Structure):
    _fields_ = [("npart",     ctypes.c_uint64),
                ("particles", ctypes.c_void_p)]


class cBeam_SoA_ctypes(ctypes.Structure):
    _fields_ = [("npart",     ctypes.c_uint64),
                ("particles", ctypes.POINTER(Particles_ctypes))]


class cBeam(object):

    clight = 299792458
    pi = 3.141592653589793238
    pcharge = 1.602176565e-19
    echarge = -pcharge
    emass = 0.510998928e6
    pmass = 938.272046e6
    epsilon0 = 8.854187817e-12
    mu0 = 4e-7*pi
    eradius = pcharge**2/(4*pi*epsilon0*emass*clight**2)
    pradius = pcharge**2/(4*pi*epsilon0*pmass*clight**2)
    anumber = 6.022140857e23
    kboltz = 8.6173303e-5 #ev K^-1 #1.38064852e-23 #   JK^-1

    @classmethod
    def from_full_beam(cls,beam):
        npart=len(beam['x'])
        particles=np.zeros(npart,particle_t)
        for nn in particle_t.names:
             particles[nn]=beam[nn]
        return cls(particles=particles)

    ptau =property(lambda p: (p.psigma*p.beta0))
    pc =property(lambda p: (p.beta*p.gamma*p.mass0))
    energy =property(lambda p: (p.gamma*p.mass0))
    energy0=property(lambda p: (p.gamma0*p.mass0))

    def __init__(self,npart=None,mass0=None,p0c=450,q0=1.0,particles=None):
        if mass0 is None:
            mass0 = cBeam.pmass
        if particles is None:
            self.npart=npart
            self.particles=np.zeros(npart,particle_t)
            self.particles['mass0']=mass0
            energy0=np.sqrt(p0c**2+mass0**2)
            gamma0=energy0/mass0
            beta0=p0c/mass0/gamma0
            chi=1.
            self.particles['partid']=np.arange(npart)
            self.particles['chi']=chi
            self.particles['beta0']=beta0
            self.particles['gamma0']=gamma0
            self.particles['p0c']=p0c
            self.particles['rvv']=1.
            self.particles['rpp']=1.
        else:
            self.particles=particles.view(particle_t)
            self.npart=len(self.particles)

    def ctypes(self):
        cdata=cBeam_ctypes(self.npart,self.particles.ctypes.data)
        return ctypes.pointer(cdata)

    def copy(self):
        return self.__class__(particles=self.particles.copy())

    def __getitem__(self,kk):
        particles=self.particles.copy().__getitem__(kk)
        return self.__class__(particles=particles)

    def get_size(self):
        return self.npart*particle_t.itemsize/8

    def __getattr__(self,kk):
        return self.particles[kk]

    def __dir__(self):
        return sorted(particle_t.names)

    def compare(self, ref, exclude=['s', 'elemid'], include=[], verbose=True):
        npart = self.particles.np
        if npart == self.particles.np:
            names = list(particle_t.names)
            for nn in exclude:
                names.remove(nn)
            for nn in include:
                names.append(nn)
            general = 0
            partn = 1
            fmts = "%-12s: %-14s %-14s %-14s %-14s"
            fmtg = "%-12s:  global diff  %14.6e"
            lgd = ('Variable', 'Reference', 'Value', 'Difference', 'Relative Diff')
            lgds = True
            fmt = fmts.replace('-14s', '14.6e')
            for pval, pref in zip(self.particles.flatten(), ref.particles.flatten()):
                pdiff = 0
                for nn in names:
                    val = pval[nn]
                    ref = pref[nn]
                    diff = ref-val
                    if abs(diff) > 0:
                        if abs(ref) > 0:
                            rdiff = diff/ref
                        else:
                            rdiff = diff
                        if lgds:
                            if verbose: print(fmts%lgd); lgds=False
                            if verbose: print(fmt%(nn,ref,val,diff,rdiff))
                        pdiff += rdiff**2
                if pdiff > 0:
                    pl = 'Part %d/%d'%(partn,npart)
                    if verbose: print(fmtg%(pl,np.sqrt(pdiff)))
                    general += pdiff
                partn += 1
            return np.sqrt(general)
        else:
          raise ValueError("Shape ref not compatible")

    def shape(self):
        return self.particles.shape

    def reshape(self,*args):
        return self.__class__(particles=self.particles.reshape(*args))


# TODO cBeam and cBeam_SoA should inherit from a common ancestor

class cBeam_SoA(object):

    clight = 299792458
    pi = 3.141592653589793238
    pcharge = 1.602176565e-19
    echarge = -pcharge
    emass = 0.510998928e6
    pmass = 938.272046e6
    epsilon0 = 8.854187817e-12
    mu0 = 4e-7*pi
    eradius = pcharge**2/(4*pi*epsilon0*emass*clight**2)
    pradius = pcharge**2/(4*pi*epsilon0*pmass*clight**2)
    anumber = 6.022140857e23
    kboltz = 8.6173303e-5  # ev K^-1 #1.38064852e-23 #   JK^-1

    @classmethod
    def from_full_beam(cls, beam):
        npart = len(beam['x'])
        particles = Particles(npart)
        for nn, prop, in particles.__dict__.iteritems():
            if nn == 'np': continue  # np doesn't exist in beam
            prop[:] = beam[nn]
        return cls(particles=particles)

    ptau = property(lambda p: (p.psigma*p.beta0))
    pc = property(lambda p: (p.beta*p.gamma*p.mass0))
    energy = property(lambda p: (p.gamma*p.mass0))
    energy0 = property(lambda p: (p.gamma0*p.mass0))

    def __init__(self, npart=None, mass0=None,
                 p0c=450, q0=1.0, particles=None):
        if mass0 is None:
            mass0 = cBeam.pmass
        if particles is None:
            self.npart = npart
            self.particles = Particles(npart)
            self.particles.mass0[:] = mass0
            energy0 = np.sqrt(p0c**2+mass0**2)
            gamma0 = energy0/mass0
            beta0 = p0c/mass0/gamma0
            chi = 1.
            self.particles.partid[:] = np.arange(npart)
            self.particles.chi[:] = chi
            self.particles.beta0[:] = beta0
            self.particles.gamma0[:] = gamma0
            self.particles.p0c[:] = p0c
            self.particles.rvv[:] = 1.
            self.particles.rpp[:] = 1.
        else:
            self.particles = particles
            self.npart = self.particles.np

    def ctypes(self):
        # NB: added np.ascontiguousarray to ensure arrays are not fragmented
        #     when passing into C.
        cdata = cBeam_SoA_ctypes(self.npart, ctypes.pointer(Particles_ctypes(
            *[np.ctypeslib.as_ctypes(
                np.ascontiguousarray(self.particles.__getattribute__(f)))
              for f, _ in Particles_ctypes._fields_]
        )))
        return ctypes.pointer(cdata)

    def copy(self):
        return self.__class__(particles=self.particles.copy())

    def __getitem__(self, kk):
        particles = self.particles.copy()[kk]
        return self.__class__(particles=particles)

    # FIXME: check this and implement correctly
    def get_size(self):
        return self.npart*particle_t.itemsize/8

    def __getattr__(self, kk):
        return self.particles[kk]

    def __dir__(self):
        return sorted(particle_t.names)

    def compare(self, ref, exclude=['s', 'elemid'], include=[], verbose=True):
        npart = self.particles.np
        if npart == self.particles.np:
            names = list(particle_t.names)
            for nn in exclude:
                names.remove(nn)
            for nn in include:
                names.append(nn)
            general = 0
            partn = 1
            fmts = "%-12s: %-14s %-14s %-14s %-14s"
            fmtg = "%-12s:  global diff  %14.6e"
            lgd = ('Variable', 'Reference', 'Value', 'Difference', 'Relative Diff')
            lgds = True
            fmt = fmts.replace('-14s', '14.6e')
            for pval, pref in zip(self.particles.flatten(), ref.particles.flatten()):
                pdiff = 0
                for nn in names:
                    val = pval[nn]
                    ref = pref[nn]
                    diff = ref-val
                    if abs(diff) > 0:
                        if abs(ref) > 0:
                            rdiff = diff/ref
                        else:
                            rdiff = diff
                        if lgds:
                            if verbose: print(fmts%lgd); lgds=False
                            if verbose: print(fmt%(nn,ref,val,diff,rdiff))
                        pdiff += rdiff**2
                if pdiff > 0:
                    pl = 'Part %d/%d'%(partn,npart)
                    if verbose: print(fmtg%(pl,np.sqrt(pdiff)))
                    general += pdiff
                partn += 1
            return np.sqrt(general)
        else:
          raise ValueError("Shape ref not compatible")

    def shape(self):
        return self.particles.shape

    def reshape(self, *args):
        return self.__class__(particles=self.particles.reshape(*args))
