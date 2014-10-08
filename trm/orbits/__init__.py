#!/usr/bin/env python

"""
a module to calculate elliptical orbits. Includes a class to define orbits, and 
various routines to deal with elliptical orbits as well.

Orbit basics
============

Useful things it is hard to remember with orbits: important planes to imagine
are the orbital plane and the plane of the sky. They intersect along a line
called the 'line of nodes'. At each end of this line orbiting objects either
come towards or away from us. The one where they go away is called the
'ascending node'.  The position angle of the ascending node is one of the
angles that defines the orientation of the orbit on the sky. It matters for
visual binaries but is indeterminate for others. Also for visual binaries, the
orbital inclination ranges from 0 to 180. 0 to 90 is reserved for binaries
that rotate counter-clockwise on the sky. The 'argument of periapsis' (omega)
is the angle measured from the ascending node to the periastron in the
direction of the orbit.

The periastron becomes ill-defined for small eccentricities thus the routines 
have some options to reference epoch to the ascending node instead as this helps 
MCMC runs.

Some routines construct orbits from a dictoinary of parameters. This should be defined
as follows (this is a complete lsit and not all are needed for all routines):

Ephemeris parameters:
 
 t0      -- zero point of binary ephemeris, days
 period  -- orbital period, days
 quad    -- quadratic term, days (assumed = 0 if not set)

General parameters:

 coord   -- coordinate system, 'Jacobi', 'Astro' or 'Marsh'
 tstart  -- reference time at which the rvms are computed (days)
 mass0   -- mass of object 0 (the binary for eclipsers) (solar)
 rint0   -- interaction radius for object 0 (AU)
 integ   -- True to use N-body integration, False for Kepler orbits
 gamma   -- systemic velocity (AU/day)

Planet parameters (the '1' becomes '2' for planet 2 etc):

 mass1    -- mass (solar)
 a1       -- semi-major axis (AU)
 epoch1   -- epoch at periastron if eperi1 is True, at ascending node if False
 e1       -- eccentricity if eomega1 is True
 omega1   -- argument of periastron (rads) if eomega1 is True
 sqeco1   -- sqrt(e)*cos(omega) if eomega1 is False
 sqeso1   -- sqrt(e)*sin(omega) if eomega1 is False
 iangle1  -- orbital inclination (rads)
 Omega1   -- longitude of ascending node (rads)
 eperi1   -- flag, see above
 eomega1  -- flag, see above

The module provides three (approximate) parameterisations of N-body
orbits. One uses 'Jacobi' coordinates, one astrocentric ('Astro') coordinates,
the other Jacobi coordinates with a small change of the relationship between
semi-major axis and orbital frequency designed to reflect a hierarchical
((star,planet1),planet2) ordering. These last are called 'Marsh' for
short. 'Jacobi' coordinates are hierarchical, taking a set of objects in
groups, first m0, then m0+m1, then m0+m1+m2 etc. Successive coordinate vectors
point from the centre of mass of one group to the next object to be included.
e.g. the first points from 0 to 1, the second points from the CoM of (0+1) to
2, etc. 'k' factors referred to in the documentation are calculated as
m1/(m0+m1), m2/(m0+m1+m2) etc and allow one to get the motion of object 0 in
reflex to the others.

'Astro' coordinates are simply referenced to the centre of mass of object 0.
'k' is calculated as m1/(m0+m1), m2/(m0+m2), in the case of ptolorb, but k =
m1/m0 in the case of prorvm (internal use only) since the vectors are scaled
to barycentric by factors of m0/m0+m1 etc. They are the same as Jacobi in the
case of 2 bodies.

Experiments I have carried out suggest relatively little to choose between
them.
"""

from _orbits import *
import math as m
import copy
import numpy as np
import trm.subs as subs
from trm.subs import Vec3

class Orbit(object):

    """
    Class to contain parameters needed to define an elliptical orbit. It
    is designed with MCMC runs in mind, and has a couple of options to facilitate
    a better choice of parameterisation.

    These are the attributes:

      a        -- semi-major axis, (AU)
      mu       -- gravitational parameter = (a**3)*(n**2) where n = angular frequency,
                  scaled by the solar value (solar masses)
      epoch    -- epoch of periastron if eperi is True, of ascending node if not, (days).
      e        -- eccentricity, if eomega is True. 
      omega    -- argument of periastron, if eomega is True, (radians).
      iangle   -- orbital inclination, (radians).
      Omega    -- longitude of ascending node, (radians).
      eperi    -- True/False flag, see epoch.
    """

    def __init__(self, pd, num, eomega):
        """
        Constructor from a dictionary, with option to supply an add-on number 
        'num' to distinguish multiple orbits (set to None if you don't want
        to add a number). Checks for existence of various parameters but 
        not much else. 

        pd     -- dictionary of parameters.
        num    -- orbit number, can be 'None'
        eomega -- If True, e and omega will be set directly, else parameter sqeco and sqeso
                  representing sqrt(e)*cos(omega) and sqrt(e)*sin(omega) will be looked for.
        """

        self.num = num
        if num is not None:
            sn = str(num)
        else:
            sn = ''

        self.a      = pd['a'+sn]
        self.mu     = pd['mu'+sn]
        self.epoch  = pd['epoch'+sn]
        if eomega:
            e     = pd['e'+sn]
            omega = pd['omega'+sn]
        else:
            e     = pd['sqeco'+sn]**2 + pd['sqeso'+sn]**2
            omega = m.atan2(pd['sqeso'+sn], pd['sqeco'+sn])            
        self.e      = e
        self.omega  = omega
        self.Omega  = pd['Omega'+sn]
        self.epoch  = pd['epoch'+sn]
        self.iangle = pd['iangle'+sn]
        self.eperi  = pd['eperi'+sn]
        
    def get(self, eomega):
        """
        Returns values as a numpy array, excluding the flags.

        On return the following are packed into a numpy array, in two
        possible different ways:

        if eomega:

          a,mu,epoch,e,omega,iangle,Omega 

        else:

          a,period,epoch,sqrt(e)*cos(omega),sqrt(e)*sin(omega),iangle,Omega 

        """
        if eomega:
            return np.array([self.a,self.mu,self.epoch,self.e,self.omega,self.iangle,self.Omega])
        else:
            sqeco = m.sqrt(self.e)*cos(self.omega)
            sqeso = m.sqrt(self.e)*sin(self.omega)
            return np.array([self.a,self.mu,self.epoch,sqeco,sqeso,self.iangle,self.Omega])

    def set(self, pa, eomega):
        """
        Sets orbit elements from a parameter array (but not the flags). 
        See the get function for the order expected.
        """
        if eomega:
            self.a,self.mu,self.epoch,self.e,self.omega,self.iangle,self.Omega = pa
        else:
            self.a,self.mu,self.epoch,sqeco,sqeso,self.iangle,self.Omega = pa
            self.e     = sqeco**2+sqeso**2
            self.omega = m.atan2(sqeso, sqeco) 

    def __str__(self):
        return 'a=' + str(self.a) + ' AU, mu=' + str(self.mu) + \
            ' msun, epoch=' + str(self.epoch) + ' days, e=' + str(self.e) + \
            ', omega=' + str(self.omega) + ' rads, iangle=' + str(self.iangle) + \
            ' rads, Omega=' + str(self.Omega) + ' rads'

    def __repr__(self):
        return self.__str__()

    @property
    def n(self):
        """
        angular frequency corresponding to the Orbit, in radians/day.
        """
        return comp_n(self.mu,self.a)

    @property
    def p(self):
        """
        orbital period corresponding to the Orbit, in days.
        """
        return 2.*np.pi/self.n

    def mean(self, t):
        """
        Compute mean anomaly of the Orbit at a particular time or times.
        
        t -- time to compute mean anomaly (days). Can be a numpy array.
        """

        if self.eperi:
            mean0  = 0.
        else:
            # compute phase offset from periastron to ascending node
            mean0 = true2mean(-self.omega, self.e)

        return self.n*(t-self.epoch) + mean0

    def true(self, t, acc=1.e-12):
        """
        Compute true anomaly of the Orbit at a particular time or times.
        
        t -- time to compute true anomaly for (days). Can be a numpy array.
        """
        
        manom = self.mean(t)
        return mean2true(manom, self.e, acc)

    def torv(self, t, acc=1.e-12):
        """
        Compute position and velocity of object at a particular time.
        
        t -- time to use (days).
        """
        tanom = self.true(t, acc)
        return kepl2cart(self.a,self.iangle,self.e,self.omega,self.Omega,self.mu,tanom)

class Norbit(Orbit):
    """
    Class for N-body orbits, sub-classing Orbit. This is a light 
    layer on top of the Orbit class which adds a parameter k which 
    is needed when calculating reflex motion of a 'central object' 
    due to one or more light companions. 'k' depends upon the masses
    and the exact coordinates being used (Jacobi or barycentric).

    See Orbit for a description of attributes to which this adds 'k'
    """

    def __init__(self, pd, num, eomega, k):
        Orbit.__init__(self, pd, num, eomega)
        self.k = k

    def __str__(self):
        return Orbit.__str__(self) + ', k=' + str(self.k)

def t2mean(t, epoch, period):
    """
    Computes orbital phase (rads) aka
    'mean anomaly'

    t      -- time to compute phase at (can be an array)
    epoch  -- reference for zero phase
    period -- orbital period
    """
    return 2.*np.pi*(t-epoch)/period

def true2ecc(true, e):
    """
    Returns the eccentric anomaly E given the true anomaly T using the relation

    tan(E/2) = [sqrt(1-e)/sqrt(1+e)] * tan(T/2)

    true  -- true anomaly, can be an array (radians)
    e     -- the eccentricity

    Returns the eccentric anomalies in radians
    """
    return 2.*np.arctan2(np.sqrt(1.-e)*np.sin(true/2.), np.sqrt(1.+e)*np.cos(true/2.)) 

def ecc2true(ecc, e):
    """
    Returns the true anomaly T given the eccentric anomaly E using the relation

    tan(T/2) = [sqrt(1+e)/sqrt(1-e)] * tan(E/2)

    ecc   -- eccentric anomaly, can be an array (radians)
    e     -- the eccentricity

    Returns the true anomalies in radians
    """
    return 2.*np.arctan2(np.sqrt(1.+e)*np.sin(ecc/2.), np.sqrt(1.-e)*np.cos(ecc/2.)) 

def mean2true(manom, e, acc):
    """
    Converts from the mean to the true anomaly

    manom  -- mean anomoly in radians (=2.*pi times usual orbital phase measured
              from periastron). Can be a numpy array.
    e      -- orbital eccentricity < 1
    acc    -- accuracy in radians for the solution. 

    Returns the true anomalies in radians
    """

    return ecc2true(mean2ecc(manom, e, acc), e)

def ecc2mean(ecc, e):
    """
    Returns the mean anomaly M from the eccentric anomaly E

    ecc   -- eccentric anomaly, can be an array (radians)
    e     -- the eccentricity

    Returns the mean anomalies in radians
    """
    return ecc - e*np.sin(ecc)

def true2mean(true, e):
    """
    Convert from true to mean anomaly

    true  -- the true anomaly, radians
    e     -- eccentricity

    Returns the mean anomalies in radians
    """
    return ecc2mean(true2ecc(true, e),e)

def comp_n(mu, a):
    """
    Computes orbital angular frequency given a gravitational
    parameter mu and semi-major axis a. This is to ensure
    consistency.

    mu    -- the quantity appearing in Kepler's 3rd law n = k*sqrt(mu/a**3)
             where k is Gauss' grav constant, hence mu has units
             of solar masses.
    a     -- semi-major axis
    """
    return subs.KGAUSS*m.sqrt(mu/a**3)

def tdelay(t, orbits, acc=1.e-12):
    """
    Returns predicted light-travel time delays from an object, 
    e.g. a binary, in seconds given sets of orbital elements 
    of perturbing bodies assuming pure Keplerian orbits.

    t      -- time or times (days), can be a numpy array.
    orbits -- list of Norbit objects defining the orbits. 
    acc    -- accuracy of computation of true anomalies. 

    Returns light travel-time delays in seconds.
    """

    # handle just one orbit
    if isinstance(orbits, Norbit):
        orbits = [orbits]
    
    # Three minus signs combine to give the final sign: z points towards 
    # Earth which gives one minus sign. Another comes because we are 
    # considering reflex motion. A third comes because of the way omega
    # is defined relative to the ascending node.
    tdel = t - t
    for orb in orbits:
        tanom = orb.true(t, acc)
        tdel -= subs.AU*orb.k*orb.a*(1-orb.e**2)/subs.C* \
            m.sin(orb.iangle)/(1.+orb.e*np.cos(tanom))*np.sin(tanom+orb.omega)

    return tdel

def rv(t, orbits, acc=1.e-12):
    """
    Returns predicted radial velocities of an object subject to perturbations from
    a set of extra bodies (in km/s), assuming Keplerian orbits

    t      -- time or times (days), can be a numpy array.
    orbits -- list of Norbit objects defining the orbits. 
    acc    -- accuracy of computation of true anomalies. 

    Returns radial velocities in km/s
    """

    # handle just one orbit
    if isinstance(orbits, Norbit):
        orbits = [orbits]
    
    # Three minus signs combine to give the final sign: z points towards 
    # Earth which gives one minus sign. Another comes because we are 
    # considering reflex motion. A third comes because of the way omega
    # is defined relative to the ascending node.
    rvs = t - t
    for orb in orbits:
        tanom = orb.true(t, acc)
        rvs  -= orb.n*subs.AU/subs.DAY/1000.*orb.k*orb.a/np.sqrt(1-orb.e**2)* \
            m.sin(orb.iangle)*(np.cos(tanom+orb.omega)+orb.e*np.cos(orb.omega))

    return rvs

def lprior(p):
    """
    Computes ln(prior probabilty) of a model. 
    Jeffries on the separations and masses.
    """

    # Distinguish RV from timing models from the presence of 'a'
    if 'a' in p:
        lp = -m.log(p['a'])
    else:
        # compute number of orbits
        norb,neph = norbeph(p)
        lp = 0.
        for i in range(norb):
            si = str(i+1)
            pname = 'a' + si
            lp -= m.log(p[pname])
            pname = 'mass' + si
            lp -= m.log(p[pname])
            
    return lp

def norbeph(p):
    """
    Returns number of orbits and ephemeris coefficients in a dictionary 'p' 
    which has entries like p['quad'], p['period1'], etc.

    Returns (norb,neph)
    """

    # compute number of orbits and number of ephemeris terms

    neph = 2
    if 'quad' in p:
        neph += 1

    norb = 0
    while True:
        ap = 'a' + str(norb+1)
        if ap not in p:
            break
        norb += 1

    return (norb,neph)

class Rvm(object):
    """
    For a particle stores the following attributes

      r    -- position vector (AU)
      v    -- velocity vector (AU/day)
      m    -- mass (solar), 
      ls   -- length scale to use for integrator (AU)
      vs   -- velocity scale to use for integrator (AU/day)
      ri   -- interaction radius (AU)
    """

    def __init__(self, r, v, m, ls, vs, ri):
        """
        r   -- position vector (AU, Vec3)
        v   -- velocity vector (AU/day Vec3)
        m   -- mass (solar)
        ls  -- length scale (AU)
        vs  -- velocity scale (AU/day)
        ri  -- interaction radius (AU)
        """
        self.r  = copy.copy(r)
        self.v  = copy.copy(v)
        self.m  = m
        self.ls = ls
        self.vs = vs
        self.ri = ri

    def __str__(self):
        string = '(' + str(self.r) + ', ' + str(self.v) + ', ' + \
            str(self.m) + ', ' + str(self.ls) + ', ' + str(self.vs) + ', ' + str(self.ri) + ')'
        return string

    def __repr__(self):
        return self.__str__()


def ptolorb(pd):
    """
    Converts from a dictionary of parameters pd to a list of Norbits
    suitable for computing as Kepler ellipses. Note that the exact
    interpretation of these depends upon the coordinates ('coord')
    being used, 'Jacobi', 'Astro' or 'Marsh'

    pd     -- Dictionary of parameters. See intro for a complete list. Not all
              are required for this routine to operate, but the complete list
              won't hurt.

    Returns a list of Norbits for objects 1, 2, 3 etc, but not object 0. The
    latter is assumed to be computed by reflex. For Astro, k = mass of
    planet/mass central object; see write-up of Jacobi coords for definition
    of k in that case.
    """

    norbit, nephem = norbeph(pd)
 
    mass0  = pd['mass0']
    coord  = pd['coord']
    
    msum   = mass0
    if coord == 'Astro':
        for i in range(1,norbit+1):
            stri = str(i)
            msum += pd['mass'+stri]

    pdc    = pd.copy()

    # 'n' in what follows is the conventional symbol for the angular frequency
    # of an orbit.
    lorb   = []
    for i in range(1,norbit+1):
        stri = str(i)

        # compute angular frequency -- some subtlety to the formula
        # used for the Jacobi coordinates. See my document on this.
        # The 1+k factor appears in the barycentric case because the
        # the semi-major axis is assumed to be relative to the 
        # barycentre, not the other star.
        a      = pd['a'+stri]
        mass   = pd['mass'+stri]
        if coord == 'Jacobi':
            msum += mass
            k     = mass/msum
            mu    = mass0/(1-k)
        elif coord == 'Marsh':
            msum += mass
            k     = mass/msum
            mu    = msum
        elif coord == 'Astro':
            mu    = mass0+mass
            k     = mass/msum
        else:
            raise Exception('coord = ' + coord + ' not recognised.')

        pdc['mu'+stri] = mu
        lorb.append(Norbit(pdc, i, pd['eomega'+stri], k))

    return lorb

def ptolrvm(pd):
    """
    Converts from a dictionary of parameters p to a list of rvms suitable
    for N-body integration

    pd     -- Dictionary of parameters. See intro for a complete list. Not all
              are required for this routine to operate, but the complete list
              won't hurt.

              Note that the motion of the body 0 is computed from reflex
              relative to the others.  Note also that the periods of the
              planets are *NOT* specified since they are computed consistently
              with the semi-major axes and the masses. The equation used
              depends upon the 'coord' being used.

    Returns a list of rvm objects in order 0, 1, 2 etc from which N-body
    integrations can proceed. Note that these are always in barycentric
    coordinates, whatever the input orbital elements may be referenced
    relative to and so can be used for N-body integration immediately. The rvm
    objects apply at the reference time 'tref' supplied in the dictionary pd
    which is taken to represent the time at which the orbital elements are
    correct.
    """

    norbit, nephem = norbeph(pd) 

    # tref is the time at which the orbital elements are assumed to apply
    mass0  = pd['mass0']
    tref   = pd['tstart']
    coord  = pd['coord']
    
    r0, v0 = Vec3(), Vec3()
    msum   = mass0
    if coord == 'Astro':
        for i in range(1,norbit+1):
            stri = str(i)
            msum += pd['mass'+stri]
    lrvm   = []
    pdc    = pd.copy()

    # 'n' in what follows is the conventional symbol for the angular frequency
    # of an orbit.
    ks     = []
    for i in range(1,norbit+1):
        stri = str(i)

        # compute angular frequency
        a    = pd['a'+stri]
        mass = pd['mass' + stri]
        if coord == 'Jacobi':
            msum += mass
            k     = mass/msum
            mu    = mass0/(1-k)
        elif coord == 'Marsh':
            msum += mass
            k     = mass/msum
            mu    = msum
        elif coord == 'Astro':
            mu    = mass0+mass
            k     = mass/msum
        else:
            raise Exception('Unrecognised coordinates in ptolrvm')

        n = comp_n(mu,a)
        pdc['mu'+stri] = mu

        orb  = Orbit(pdc,i,pdc['eomega'+stri])
        r,v  = orb.torv(tref)

        # accumulate reflex sums (automatically barycentric)
        r0  -=  k*r
        v0  -=  k*v

        # store in Rvm list, store k values
        lrvm.append(Rvm(r, v, mass, a, n*a, pdc['rint'+stri]))
        ks.append(k)

    if coord == 'Jacobi' or coord == 'Marsh':
        # Need to convert the Jacobi coordinates to barycentric ones
        # for N-body work. Work through rvm list in reverse order:
        rsum, vsum = Vec3(), Vec3()
        for i in range(len(ks)-1,-1,-1):
            rsum      += ks[i]*lrvm[i].r
            vsum      += ks[i]*lrvm[i].v
            lrvm[i].r -= rsum
            lrvm[i].v -= vsum

    elif coord == 'Astro':
        # to get from astro to barycentric simply add r0, v0
        for i in range(len(ks)):
            lrvm[i].r += r0
            lrvm[i].v += v0

    # Create and insert the zeroth object Rvm and return
    rvm0 = Rvm(r0, v0, mass0, r0.norm(), v0.norm(), pdc['rint0'])
    lrvm.insert(0,rvm0)
    return lrvm

def etimes(pd, cycle, acc=1.e-12, nmax=10000, stoerm=True):
    """
    Computes eclipse times for a given set of ephemeris and third body
    parameters and cycle numbers.

    This works in two possible ways according to whether pd['integ'] is True
    or False

      1) The orbits (entered via p) are computed as straight ellipses using
         Kepler's laws alone. (integ = False)

      2) The orbits are assumed to be due to third-bodies. The parameters
         entered are translated into orbits of the third bodies and the binary
         and then integrated using N-body Newtonian physics. In this case a
         start time and a binary star total mass are required in pd.

    pd     -- dictionary of parameters. See intro.

    cycle  -- eclipse cycle numbers

    acc    -- accuracy parameter (for orbit integrations only)
    
    nmax   -- maximum number of steps to take

    stoerm -- use Stoermer's rule for integrating 2nd-order conservative
              ODEs. 2x faster than the alternative which is maintained for
              comparison purposes.
    """

    # compute number of orbits and number of ephemeris terms
    norbit, nephem = norbeph(pd) 

    tmodel = pd['t0'] + pd['period']*cycle
    if nephem == 3: 
        tmodel += pd['quad']*cycle**2
    elif nephem != 2:
        print 'ERROR: invalid ephemeris length =',nephem,'encountered in etimes'
        exit(1)

    if pd['integ']:
        lrvm = ptolrvm(pd)
        # integrate
        ttime,ntest,eratio,npoint,ierr,tnext,nstore,arr = \
            integrate(lrvm, tmodel-pd['tstart'], acc=acc, nmax=nmax, stoerm=stoerm)

        tmodel -= subs.AU*arr[:,3]/subs.C/subs.DAY

    else:
        orbs    = ptolorb(pd)
        tmodel += tdelay(tmodel, orbs, acc)/subs.DAY

    return tmodel

def read_model(mfile):
    """
    Reads in an orbit model file from file called mfile

    Returns (head, vvar, pars, model)

    head  --- list of header strings
    vvar  --- list of variable parameter names
    pars  --- list of all parameter names
    model --- dictionary of model values keyed by name
    """

    with open(mfile) as fm:
        model = {}
        vvar  = []
        head  = []
        pars  = []
        for line in fm:
            if line.startswith('#'):
                head.append(line)
            elif line.find('=') > -1:
                lsplit = line.split()
                if len(lsplit) == 3:
                    name,eq,value = lsplit
                    if name != 'coord' and name != 'mass0' and name != 'integ' \
                            and not name.startswith('eperi') and not name.startswith('eomega') \
                            and not name.startswith('rint') and name != 'tstart':
                        print 'Model line = "' + line + '" has invalid format (missing "v" or "f")'
                        exit(1)
                    vorf = 'f'
                elif len(lsplit) == 4:
                    name,eq,value,vorf = lsplit
                else:
                    print 'Model line = "' + line + '" has invalid format.'
                    exit(1)

                if eq != '=' or (vorf != 'v' and vorf != 'f' and vorf != 'V' and vorf != 'F'):
                    print 'Model line = "' + line + '" has invalid format.'
                    exit(1)

                model[name] = value
                pars.append(name)
                if vorf == 'v' or vorf == 'V':
                    vvar.append(name)

    model = fixpar(model)
    return (head, vvar, pars, model)

def ephmod(cycle, t0, period, quad=0.):
    """
    Compute times given ephemeris
    
    cycle  -- cycle numbers
    t0     -- zeropoint of ephemeris
    period -- period
    quad   -- quadratic term
    """
    return t0 + cycle*(period + quad*cycle)

def fixpar(pzero, vvars=None):
    """
    Correct dictionary of orbital parameters. Makes them floating point if need be
    and updates from lists of variable names and values. 

    pzero  -- dictionary of parameter name / values
    vvars  -- list of name / value pairs of parameters to be updated (e.g. variable ones)
              Could be made from separate lists with e.g. zip(names, values)

    Returns a dictionary similar to pzero with variables updated by p and converted to right type.
    """

    ptry = dict(pzero)
    if vvars is not None:
        for vn,vv in vvars:
            if vn in ptry:
                ptry[vn] = vv

    # Convert strings to floats, except for special cases
    for vn,vv in ptry.iteritems():
        if vn.startswith('eperi') or vn == 'integ' or vn.startswith('eomega'):
            if vv == 'True':
                ptry[vn] = True
            elif vv == 'False':
                ptry[vn] = False
            else:
                ptry[vn] = bool(int(vv))
        elif not vn.startswith('coord'):
            try:
                ptry[vn] = float(vv)
            except:
                if vv == 'PIBY2':
                    ptry[vn] = np.pi/2.
                elif vv == 'PI':
                    ptry[vn] = np.pi
                elif vv == 'TWOPI':
                    ptry[vn] = 2.*np.pi
    if 'coord' in pzero and pzero['coord'] != 'Jacobi' and \
            pzero['coord'] != 'Marsh' and pzero['coord'] != 'Astro':
        raise Exception('Do not understand coordinates = ' + pzero['coord'])
    return ptry

def check_tpar(pd):
    """
    Check parameter names of a model. Raises Exception if there are problems, otherwise
    does nothing.
    
    pd -- dictionary of parameter names and values
    """

    print 'Checking timing model parameters'
    if 't0' not in pd:
        raise Exception('Could not find t0 in model')
    if 'period' not in pd:
        raise Exception('Could not find period in model')
    if 'integ' not in pd:
        raise Exception('Could not find integ in model')
    if 'tstart' not in pd:
        raise Exception('Could not find tstart in model')
    if 'coord' not in pd:
        raise Exception('Could not find coord in model')
    if 'mass0' not in pd:
        raise Exception('Could not find mass0 in model')
    if 'rint0' not in pd:
        raise Exception('Could not find mass0 in model')

    norb, neph = norbeph(pd) 

    if neph == 2:
        print 'Found ephemeris coefficients t0 and period.'
    elif neph == 3:
        if 'quad' not in pd:
            raise Exception('3 ephemeris coefficients but could not find quad in model')
        print 'Found ephemeris coefficients t0, period and quad'
    else:
        raise  Exception('Number of ephemeris coefficients = ' + str(neph) + ' is invalid.')
    
    # check orbits
    for i in range(1,norb+1):
        si  = str(i)

        eomega = 'eomega'+si
        if eomega not in pd:
            raise Exception('Could not find ' + eomega + ' in timing model')
        if isinstance(pd[eomega],str):
            eomega = pd[eomega] == True
        else:
            eomega = pd[eomega]

        if eomega:
            names = ('a'+si,'epoch'+si,'Omega'+si,'e'+si,'omega'+si,'mass'+si,'eperi'+si,'rint'+si)
        else:
            names = ('a'+si,'epoch'+si,'Omega'+si,'sqeco'+si,'sqeso'+si,'mass'+si,'eperi'+si,'rint'+si)

        for name in names:
            if name not in pd:
                raise Exception('Could not find ' + name + ' in timing model')

def check_tmod(p):
    """
    Check values of a model. Returns True if all ok.
    Max eccentricity = 0.999 to eliminate problems
    solving Kepler's equation

    p -- dictionary of parameter names and values
    """
    norb, neph = norbeph(p)
    if p['period'] <= 0.:
        return False
    for i in range(1,norb+1):
        si = str(i)
        if p['a' + si] <= 0. or p['mass' + si] <= 0.:
            return False
        if p['eomega' + si]:
            e = 'e' + si
            if p[e] < 0. or p[e] >= 0.999:
                return False
        else:
            if p['sqeco' + si]**2 + p['sqeso' + si]**2 >= 0.999:
                return False

    return True

def check_rvmod(p):
    """
    Check values of a model. Returns True if all ok,

    p -- dictionary of parameter names and values
    """
    if p['period'] <= 0. or p['a'] < 0.:
        return False
    if p['type'] == 1:
        if p['e'] < 0. or p['e'] >= 1:
            return False
    else:
        if p['esx']**2+p['esy']**2 >= 0.999:
            return False

    return True

def eccres(p, cycle, time, etime, acc=1.e-12, nmax=10000):
    """
    Generates RMS normalised residuals for eclipse timing data.
    Needed for scipy.optimize.leastsq. Only returns ones with
    positive errors. 

    p      -- dictionary of parameters
    cycle  -- cycle numbers
    time   -- equivalent times
    etime  -- equivalent errors
    """

    ok  = etime > 0
    res = (time[ok] - etimes(p, cycle[ok], acc, nmax))/etime[ok]
    return res

def tchisq(p, cycle, time, etime, acc=1.e-12, nmax=10000):
    """
    Computes chi**2 for MCMC. 

    p      -- the full model in dictionary form
    cycle  -- cycle numbers of the eclipses
    time   -- corresponding times
    etime  -- uncertainties in the times
    """
    resids = eccres(p, cycle, time, etime, acc, nmax)
    return (resids*resids).sum()

def rvp(pd, times):
    """
    Computes radial velocities from a model contained in the dictionary pd.
    Returns numpy.ndarray of velocities. 

    pd      -- dictionary of parameters (full).
    times   -- observation times
    """

    # compute number of orbits and number of ephemeris terms
    norbit, nephem = norbeph(pd) 

    if pd['integ']:
        lrvm = ptolrvm(pd)
        # integrate
        ttime,ntest,eratio,npoint,ierr,tnext,nstore,arr = \
            integrate(lrvm, tmodel-pd['tstart'], acc=acc, nmax=nmax, stoerm=stoerm)

        rvs = -subs.AU*arr[:,6]/subs.DAY

    else:
        orbs = ptolorb(pd)
        rvs  = rv(times, orbs)

    return rvs+pd['gamma']

def rvchisq(pd, times, vels, verrs):
    """
    Computes chi**2 for given set of parameters given a set of velocity measurements

    pd      -- dictionary of parameters (full)
    times   -- observation times
    vels    -- radial velocities km/s
    verrs   -- velocity errors, km/s
    """

    return (((vels-rvp(pd, times))/verrs)**2).sum()

def setcen(lrvm):
    """
    Sets a list of rvm objects to the centre of mass frame. i.e. both the position
    and velocity components are adjusted to the centre of mass frame
    """

    rm = Vec3()
    vm = Vec3()
    sm = 0.
    for rvm in lrvm:
        rm += rvm.m*rvm.r
        vm += rvm.m*rvm.r
        sm += rvm.m

    rm /= sm
    vm /= sm

    for rvm in lrvm:
        rvm.r -= rm
        rvm.v -= vm

def kepl2cart(a,i,e,omega,Omega,mu,tanom):
    """
    Computes Cartesian position and velocity given orbital elements in the form:

    a        -- semi-major axis (AU)
    i        -- orbital inclination (rads)
    e        -- eccentricity
    omega    -- argument of periastron (rads)
    Omega    -- longitude of ascending node (rads)
    mu       -- gravitational parameter, see e.g. Orbit (solar masses)
    tanom    -- true anomaly (rads)

    Returns r,v

    r  -- Vec3 position (AU)
    v  -- Vec3 velocity (AU/day)
    """
    
    # cos and sin of omega + tanom 
    cospt = m.cos(omega+tanom)    
    sinpt = m.sin(omega+tanom)

    # cos and sin of tanom 
    cost = m.cos(tanom)    
    sint = m.sin(tanom)

    # cos and sin of omega 
    cosp = m.cos(omega)    
    sinp = m.sin(omega)

    # cos and sin of Omega 
    coso = m.cos(Omega)    
    sino = m.sin(Omega)

    # cos and sin of inclination
    cosi = m.cos(i)    
    sini = m.sin(i)

    l = a*(1-e**2)

    # position coordinates first
    LSCALE = l/(1+e*cost)

    x = -LSCALE*(cospt*sino+sinpt*cosi*coso)
    y =  LSCALE*(cospt*coso-sinpt*cosi*sino)
    z = -LSCALE*sinpt*sini

    # velocity coordinates second
    n      = comp_n(mu, a)
    VSCALE = n*l/(1-e**2)**1.5

    vx =  VSCALE*((sinpt+e*sinp)*sino-(cospt+e*cosp)*cosi*coso)
    vy =  VSCALE*(-(sinpt+e*sinp)*coso-(cospt+e*cosp)*cosi*sino)
    vz = -VSCALE*(cospt+e*cosp)*sini

    return (Vec3(x,y,z), Vec3(vx,vy,vz))

def cart2kepl(r,v,mu):
    """
    Computes Keplerian elements given vectors r and v and the
    relationship between semi-major axis and angular frequency.

    r        -- Vec3 position relative to focus (AU)
    v        -- Vec3 velocity relative to focus (AU/day)
    mu       -- Gravitational parameter, see Orbit (solar masses)

    Returns (a,i,e,omega,Omega,tanom)

    a        -- semi-major axis (AU)
    i        -- orbital inclination (rads)
    e        -- eccentricity
    omega    -- argument of periastron (rads)
    Omega    -- longitude of ascending node (rads)
    tanom    -- true anomaly (rads)

    If e >=1 it returns l instead of a and -1 for the period.
    """
    
    # in what follows, quantities ending with 'v' are vectors. 
    # unit indicates a unit vector
    rn     = r.norm()
    vn     = v.norm()
    runitv = r.unit()

    # radial and tangential cpts of velocity
    vr     = subs.dot(runitv, v)
    vtv    = v - vr*runitv
    vt     = vtv.norm()

    # unit vector in theta direction
    tunitv = vtv.unit()

    l   = (rn*vt/subs.KGAUSS)**2/mu
    vc  = subs.KGAUSS*m.sqrt(mu/l)
    if vr == 0. and vt == vc:
        e     = 0.
        tanom = 0.
        cosnu = 1.
        sinnu = 0.
    else:
        ecos   = vt/vc - 1.
        esin   = vr/vc
        e      = m.sqrt(ecos**2+esin**2)
        tanom  = m.atan2(esin, ecos)
        cosnu  = ecos / e
        sinnu  = esin / e

    punitv    = cosnu*runitv - sinnu*tunitv
    hunitv    = subs.cross(r, v).unit()
    xunitv    = Vec3(1.,0.,0.)
    yunitv    = Vec3(0.,1.,0.)
    zunitv    = Vec3(0.,0.,1.)
    nunitv    = subs.cross(hunitv,zunitv).unit()
    i       = m.acos(subs.dot(hunitv,zunitv))
    tvec    = subs.cross(hunitv,nunitv)
    omega   = m.atan2(subs.dot(punitv,tvec),subs.dot(punitv,nunitv))
    if omega < 0:
        omega += 2.*np.pi
    Omega = m.atan2(-subs.dot(xunitv,nunitv),subs.dot(yunitv,nunitv))
    if e < 1.:
        a       = l/(1.-e**2)
        return (a,i,e,omega,Omega,tanom)
    else:
        return (l,i,e,omega,Omega,tanom)

def amthird(asemi, mass, period):
    """
    This returns the semi-major axis and mass of a gravitating body 
    given the semi-major axis, mass and orbital period of an object 
    in orbit around it.

    asemi  -- semi-major axis of object (AU)
    mass   -- mass of object (solar)
    period -- orbital period (years)

    Returns (at, mt) semi-major axis and mass of third body (AU, solar).

    asemi and period can be numpy arrays
    """

    mt   = period - period
    temp = asemi**3/period**2
    while True:
        mtold = mt
        mt    = (temp*(mass+mt)**2)**(1./3.)
        if np.all(np.abs(mt-mtold) < 1.e-15*mt):
            break

    return ((mass/mt)*asemi, mt)

def load_times(fname):
    """
    Loads times from a timing file given a file name fname.

    Returns cycle,time,etime,offset,soften where
    
    cycle  -- array of cycle numbers, offset if set.
    time   -- the corresponding times
    etime  -- the uncertainties in the times, softened if set.
    offset -- the cycle number offset applied
    soften -- the softening factor added in quadrature to give etime.
    """

    data    = np.loadtxt(fname)

    # check for a cycle number offset and softening factor
    fp = open(fname)
    offset = 0
    soften = 0.
    for line in fp:
        if line.startswith('# offset = '):
            offset = int(line[11:])
        if line.startswith('# soften = '):
            soften = float(line[11:])
    fp.close()
    print 'Cycle number offset =',offset,'will be applied.'
    print 'Softening factor    =',soften,'will be applied.'
    cycle   = data[:,0] - offset
    time    = data[:,1]
    etime   = np.sqrt(soften**2 + data[:,2]**2)
    asort   = np.argsort(cycle)
    return cycle[asort], time[asort], etime[asort], offset, soften

class Times(object):
    """"
    Represents eclipse times from a file.

    Attributes:

    data    -- cycle numbers, eclipse times and errors
    offset  -- offset to be applied to cycle numbers
    soften  -- softening factor to ameliorate uncertainties
    """

    def __init__(self, fname):
        self.data    = np.loadtxt(fname)

        # check for a cycle number offset and softening factor
        self.offset = 0
        self.soften = 0.
        with open(fname) as fp:
            for line in fp:
                if line.startswith('# offset = '):
                    self.offset = int(line[11:])
                if line.startswith('# soften = '):
                    self.soften = float(line[11:])

        print 'Cycle number offset =',self.offset,'will be applied.'
        print 'Softening factor    =',self.soften,'will be applied.'

    def __len__(self):
        """
        Defines a length (number of times)
        """
        return len(data)

    def get(self):
        """
        Returns cycle numbers, times and errors with corrections applied
        """
        return (self.data[:,0] - self.offset, self.data[:,1], np.sqrt(self.soften**2 + self.data[:,2]**2))

    def set_times(self, times):
        """
        Sets the times
        """
        self.data[:,1] = times

    def write(self, fname):
        """
        Dumps to a file.
        """
        if self.data[:,1].max() > 2400000.:
            fmt = '%17.9f'
        else:
            fmt = '%15.9f'
        with open(fname,'w') as fp:
            fp.write('# offset = %d\n' % (self.offset,))
            fp.write('# soften = %9.3e\n' % (self.soften,))
            np.savetxt(fp, self.data, '%d ' + fmt + ' %9.3e')

def trans_code(code):
    """
    Interprets return code from integrate and composes suitable message
    
    code -- code indicating collision, escape etc.
    """
    if code == 0:
        return 'integration was fully completed.'
    elif code < 10:
        return 'object ' + str(code) + ' had too much kinetic energy.'
    elif code < 1000:
        n  = code // 10
        n1 = n % 10
        n2 = (n - n1) // 10
        return 'objects ' + str(n1) + ' and ' + str(n2) + ' collided.'
    else:
        return 'object ' + str(code // 1000) + ' exceeded escape radius.'

        
def xyz(pd, times, acc=1.e-12, nmax=10000, stoerm=True):
    """
    Computes x,y,z coordinates for a set of bodies with orbital
    elements specified in pd

    pd     -- dictionary of parameters.

    times  -- the times in days for which the coordinates will be returned

    acc    -- accuracy parameter (for orbit integrations only)
    
    nmax   -- maximum number of steps to take

    stoerm -- use Stoermer's rule for integrating 2nd-order conservative
              ODEs. 2x faster than the alternative which is maintained for
              comparison purposes.

    Returns arrays of x,y,z for each body. Thus a 2 planet orbit will return
    9 numpy arrays x,y,z for the central object, planet 1, planet 2.
    Units of AU.
    """

    # compute number of orbits and number of ephemeris terms
    norbit, nephem = norbeph(pd) 

    mass0 = pd['mass0']
    if pd['integ']:
        # Newtonian
        lrvm = ptolrvm(pd)
        # integrate
        ttime,ntest,eratio,npoint,ierr,tnext,nstore,arr = \
            integrate(lrvm, times-pd['tstart'], acc=acc, nmax=nmax, stoerm=stoerm)
        ret = [arr[:,1],arr[:,2],arr[:,3]]
        for nb in range(norbit):
            ind = 6*(nb+1)
            ret += [arr[:,ind+1],arr[:,ind+2],arr[:,ind+3]]
    else:
        # Keplerian
        orbs    = ptolorb(pd)

        x0 = np.zeros_like(times)
        y0 = np.zeros_like(times)
        z0 = np.zeros_like(times)
        ret = []
        for orb in orbs:
            tanom = orb.true(times, acc)
            scale = orb.a*(1-orb.e**2)/(1.+orb.e*np.cos(tanom))
            cto   = np.cos(tanom+orb.omega)
            sto   = np.sin(tanom+orb.omega)
            x     = -scale*(m.sin(orb.Omega)*cto+m.cos(orb.Omega)*m.cos(orb.iangle)*sto)
            y     =  scale*(m.cos(orb.Omega)*cto-m.sin(orb.Omega)*m.cos(orb.iangle)*sto)
            z     = -scale*m.sin(orb.iangle)*sto
            ret  += [x,y,z]
            x0 -= orb.k*x
            y0 -= orb.k*y
            z0 -= orb.k*z

        ret = [x0,y0,z0] + ret

    return ret
