#!/usr/bin/env python

"""
helps carry out MCMC work by defining a few classes that are supposed to be of
generic utility. It is mainly designed at the moment around light curve models
from the 'lcurve' suite of C++ routines, acting as a wrapper, but it has been
used for very different cases as well.

Classes
=======

Chain     -- contains the results of an MCMC run, and the information needed to re-start iterations if needed.
Jump      -- class to help in generation of new models. Basically a container for a covariance matrix.
Lcurve    -- contains an lcurve model
MCMCError -- exception class for the module


Functions
========

merge_chains -- combines chains together into one.

Scripts
=======

Most users will only be interested in the scripts (all end in .py):

cmin   -- prints minimum chi**2 and -2*ln(P) of chains
corr   -- plots auto-correlation for a parameter
covar  -- computes covariance resulting in a '.lmq' file
diags  -- convergence diagnostics (Gelman & Rubin 1992)
genlmq -- generates and updates initial jump distributions.
lcmcmc -- main script for running MCMC using the light curve fitter routine 'lroche'
masses -- ?
modlmq -- tempers a jump distribution
oned   -- plots 1D histogram
oneds  -- plots 1D histogram with confidence regions coloured
pall   -- plots p1 vs p2 diagrams for all parameters
series -- plots time-series of given parameter
twod   -- plots p1 vs p2 for any two parameters
twods  -- plots p1 vs p2 for any two parameters with colours

Usage: for documentation on usage of the above scripts, look at the help in lcmcmc.py (e.g. 'lcmcmc.py -h | more'),
then look at modlmq.py similarly.

"""

import os, os.path, subprocess, copy
import math as m
import numpy as np
from trm import subs, roche

# Look for lroche
if 'TRM_SOFTWARE' not in os.environ:
    print 'Environment variable "TRM_SOFTWARE" is not defined.'
    exit(1)

# the lroche command path
lroche = os.path.join(os.environ['TRM_SOFTWARE'], 'bin', 'lcurve', 'lroche')

if not os.path.exists(lroche):
    print 'Could not find',lroche
    exit(1)


def _chain_from_file(fname):
    """
    Carries out start up stuff for reading a chain from a file. This is 
    to allow some code re-use without getting too fancy.

    fname -- file name
    
    Returns (fobj,names,model,nstore,nwalker,stretch,method,chmin,jump)

    fobj     -- file object open for reading at the start of the data.
    names    -- names of stored items
    model    -- dictionary of the model parameters
    nstore   -- how often data are stored
    nwalker  -- number of walkers for affine method
    stretch  -- stretch factor for affine
    method   -- 'a' for affine, 'm' metropolis
    chmin    -- input minimum chi**2
    jump     -- jump distribution
    """

    fobj    = open(fname)

    # First read variables and nstore
    model   = {}
    namread = False
    nstread = False
    nwlread = False
    mtread  = False
    modread = False
    strread = False
    jfound  = False

    for line in fobj:
        if line[0:2] == '##':
            names = line[2:].split()
            names = [n.strip() for n in names]
            namread = True
        elif line.startswith('# nstore'):
            nstore = int(line[line.find('=')+1:])
            nstread = True
            modread = False
        elif line.startswith('# nwalker'):
            nwalker = int(line[line.find('=')+1:])
            nwlread = True
            modread = False
        elif line.startswith('# stretch'):
            stretch = float(line[line.find('=')+1:])
            strread = True
            modread = False
        elif line.startswith('# method'):
            method = line[line.find('=')+1:].strip()
            mtread = True
            modread = False
        elif line.startswith('# Minimum chi**2'):
            chmin  = float(line[line.find('=')+1:])
        elif line.startswith('# Initial model:'):
            modread = True
        elif modread and line.startswith('#') and '=' in line:
            key = line[2:line.find('=')].strip()
            val = line[line.find('=')+1:].strip()
            model[key] = val
        elif line.startswith('# Jump start'):
            jfound = True
            break

    if not namread or not nstread:
        raise MCMCError('_chain_from_file: failed to read either the parameter names or nstore')

    if not nwlread and not modread and not strread and not modread:
        print 'Old style file assumed'
        method  = 'm'
        nwalker = 1
        stretch = 2.

    if not jfound:
        raise MCMCError('_chain_from_file: failed to find start of jump distribution')

    # Then read the jump distribution
    try:
        jump = Jump(fobj)
    except MCMCError, err:
        jump = None
    return (fobj,names,model,nstore,nwalker,stretch,method,chmin,jump)

class Chain(object):
    
    """
    Represents the results of an MCMC run.

    Attributes:

    names   -- list of parameter names, possibly including extras such as chisq, lp, wdwarf. The first of the
               extras should always be 'chisq' to allow identification of the variable parameters.
    vals    -- the values, a 2D numpy array.
    jump    -- the jump distribution used if the Metroplis-Hastings method is being used (a Trial object)
    nstore  -- how often data is stored.
    chmin   -- minimum chi**2, used if the errors need re-scaling.
    model   -- dictionary defining model
    method  -- the method being used
    nwalker -- the number of walkers for the affine method
    stretch -- the stretch factor for the affine method.
    nmods   -- number of models stored. Not same as len(vals) because one can skip entries to save memory.
    """
    
    def __init__(self, *args):
        """
        Initialises a Chain. This return the minimum chi**2 if read from
        a file (-1 if started from scratch)

        One arguments:

         name  -- a file name containing a previously stored run

        Two arguments:

         name  -- a file name containing a previously stored run
         nskip -- number of point groups to skip (-n to just keep the final n groups)

        Three arguments:

         name  -- a file name containing a previously stored run
         nskip -- number of point groups to skip (-ve to keep the final nskip groups)
         nthin -- number to skip over per point accepted

        Eight arguments:

         model   -- dictionary defining all parameters and initial values. These can include both
                    variable and non-variable parameters and are simply included for reference
         nstore  -- how often models have been written (15 = 1 in every 15)
         method  -- 'm' for Metropolis-Hastings jumps using 'jump' or 'a' for the affine method.
         jump    -- the jump distribution
         nwalker -- number of "walkers" if using affine method.
         stretch -- stretch factor for affine method
         extras  -- extra parameters in addition to the variable parameters contained in jump
                    should start with 'chisq' for some routines to work.
         cmin    -- minimum chi**2. This used to scale the uncertainties.
        """

        if len(args) == 1 or len(args) == 2 or len(args) == 3:

            # Read headers
            fobj,self.names,self.model,self.nstore,self.nwalker,self.stretch,\
                self.method,self.chmin,self.jump = _chain_from_file(args[0])

            if len(args) == 1:
                nskip = 0
                nthin = 0  
            elif len(args) == 2:
                nskip = int(args[1])
                nthin = 0
            elif len(args) == 3:
                nskip = int(args[1])
                nthin = int(args[2])

            # Read the data
            data = []
            ns = 0
            for line in fobj:
                if not line[0:1] == '#':

                    if nskip >=0 and ns >= (nthin+1)*nskip*self.nwalker and (ns % (nthin+1)) == 0:
                        vals = line.split()
                        if len(vals) == len(self.names):
                            data.append(np.cast[float](vals))
                        else:
                            print 'Line = ' + line
                            print 'has',len(vals),'values compared to expected number =',len(self.names)
                            break

                    if nskip < 0 and (ns % (nthin+1)) == 0:
                        vals = line.split()
                        if len(vals) == len(self.names):
                            data.append(np.cast[float](vals))
                            if len(data) > -nskip*self.nwalker:
                                data.pop(0)
                        else:
                            print 'Line = ' + line
                            print 'has',len(vals),'values compared to expected number =',len(self.names)
                            break
                    ns += 1

            fobj.close()
            self.nmods = ns
            self.vals  = np.atleast_2d(np.array(data))

        elif len(args) == 8:

            self.model   = dict(args[0])
            self.nstore  = int(args[1])
            if args[2] == 'm' or args[2] == 'a':
                self.method  = args[2]
            else:
                raise MCMCError('Chain.__init__: could not understand method which must be "m" or "a"')

            self.nstore  = int(args[1])
            if args[3] is not None:
                self.jump    = copy.deepcopy(args[3])
            else:
                self.jump    = None
            self.nwalker = int(args[4])
            self.stretch = float(args[5])
            if self.jump is not None:
                self.names   = self.jump.par + list(args[6])
            else:
                self.names   = list(args[6])
            self.chmin   = float(args[7])
            self.vals    = None
            self.nmods   = 0
        else:
            raise MCMCError('Chain.__init__: could not interpret arguments of Chain constructor; expecting 1 or 8')

    def add_file(self, fname):
        """
        Adds in data from another file
        """
        fp = open(fname)

        for line in fp:
            if line[0:2] == '##':
                names = line[2:].split()
                names = [n.strip() for n in names]
                if names != self.names:
                    raise MCMCError('Chain.add_file: file ' + fname + ' contains a different set of parameters or parameters' \
                                    + ' in a different order from those already read.')
            elif line.startswith('# nstore = '):
                nstore = int(line[line.find('=')+1:].strip())
                if nstore != self.nstore:
                    raise MCMCError('Chain.add_file: nstore of file ' + fname + ' differs from that already loaded')

            elif line.startswith('# nwalker = '):
                nwalker = int(line[line.find('=')+1:].strip())
                if nwalker != self.nwalker:
                    raise MCMCError('Chain.add_file: nwalker of file ' + fname + ' differs from that already loaded')

            elif line.startswith('# method = '):
                method = line[line.find('=')+1:].strip()
                if method != self.method:
                    raise MCMCError('Chain.add_file: method of file ' + fname + ' differs from that already loaded')

            elif line.startswith('# stretch = '):
                stretch = float(line[line.find('=')+1:].strip())
                if stretch != self.stretch:
                    raise MCMCError('Chain.add_file: stretch of file ' + fname + ' differs from that already loaded')

            elif line.startswith('# Jump start'):
                break

        # Then read and check the jump distribution
        tjump = Jump(fp)
        if tjump != self.jump:
            raise MCMCError('Chain.add_file: jump distribution of file ' + fname + ' differs from that already loaded')

        # Read the data
        data = []
        for line in fp:
            if not line[0:1] == '#':
                vals = line.split()
                if len(vals) == len(self.names):
                    data.append(np.cast[float](vals))
                else:
                    break

        fp.close()

        # add onto end of current data
        if self.vals is not None:
            self.vals = np.concatenate((self.vals,np.array(data)),0)
            self.nmods += len(data)
        else:
            self.vals = np.atleast_2d(np.array(data))
            self.nmods = len(data)

    def add_line(self, values):
        """
        Adds in another line of data (which must be identical in format to currently stored lines,
        although only the number of items is checked)
        """
        # convert to floats
        vls = np.array(map(float, values))
        if len(vls) != len(self.names):           
            raise MCMCError('Chain.add_line: number of values supplied clashes with number of parameters in Chain')
        if self.vals is not None:
            self.vals = np.vstack((self.vals,vls))
        else:
            self.vals = np.atleast_2d(vls)
        self.nmods += 1

    def genjump(self):
        """
        Compute the Jump corresponding to the Chain.
        """

        # compute covariances
        lp     = len(self.jump.par)
        mean   = self.vals.mean(0)[:lp]
        dvals  = np.matrix(self.vals[:,:lp] - mean)
        tdvals = dvals.transpose()
        cov    = tdvals*dvals/len(dvals)
        sigma  = np.sqrt(np.diag(cov))
        corr   = ((cov/sigma).transpose()/sigma).transpose()
        return Jump(self.jump.par, sigma, corr)

    def mean(self):
        """
        Return means corresponding to the chain
        """

        lp     = len(self.jump.par)
        return self.vals.mean(0)[:lp]

    def rms(self):
        """
        Return RMSs corresponding to the chain
        """

        lp     = len(self.jump.par)
        return self.vals.std(0)[:lp]


    def cov(self):
        """
        Return covariance matrix corresponding to the chain
        """

        # compute covariances
        lp     = len(self.jump.par)
        mean   = self.vals.mean(0)[:lp]
        dvals  = np.matrix(self.vals[:,:lp] - mean)
        tdvals = dvals.transpose()
        cov    = tdvals*dvals/len(dvals)
        return cov

    def chop(self, n):
        """
        Chops n points from the start of a chain.
        
        Raises MCMCErrors if there are no values or too few values
        """
        if self.vals is None:
            raise MCMCError('Cannot chop ' + str(n) + ' points as there are no values at all at present')
        if n >= len(self):
            raise MCMCError('Cannot chop ' + str(n) + ' points as there are only ' + str(len(self)) + ' in the chain.')
        # The copy is here in order (ultimately) to save on memory
        self.vals = np.copy(self.vals[n:,:])

    def chop_burn(self, cburn):
        """
        Chops burn-in from start of a chain with burn-in defined
        by the initial period that the chisq exceeds the value cburn.

        Returns (True,n) where n is the number removed if cburn has been 
        reached, else (False,0). Raises an MCMCError if no chisq is found.
        """
        if 'chisq' in self.names:
            if self.vals is not None:
                nc = self.names.index('chisq')
                c  = self.vals[:,nc]
                if len(c[c < cburn]):
                    # find first instance of c < cburn
                    n = np.arange(len(c))[c < cburn][0]
                    self.vals = self.vals[n:,:]
                    return (True,n)
                else:
                    return (False,0)
            else:
                return (False,0)
        else:
            raise MCMCError('Chain.chop_burn: no chisq found')

    def cmin(self):
        """
        Returns minimum chi**2 stored in the chain. Returns 1e30 if there
        are no values yet.

        Raises an MCMCError if no parameter named 'chisq' is found.
        """
        if 'chisq' in self.names:
            if self.vals is not None:
                nc = self.names.index('chisq')
                return self.vals[:,nc].min()
            else:
                return 1.e30
        else:
            raise MCMCError('Chain.cmin: no chisq found')

    def pmax(self):
        """
        Returns maximum -2*ln(posterior prob) stored in the chain. Returns 1e30 if there
        are no values yet.

        Raises an MCMCError if no parameter named 'pprob' is found.
        """
        if 'pprob' in self.names:
            if self.vals is not None:
                nc = self.names.index('pprob')
                return self.vals[:,nc].max()
            else:
                return 1.e30
        elif 'lpost' in self.names:
            if self.vals is not None:
                nc = self.names.index('lpost')
                return self.vals[:,nc].max()
            else:
                return -1.e30
        else:
            raise MCMCError('Chain.pmin: no pprob or lpost found')


    def tofile(self, fname, fmt='%.16e'):
        """
        Dumps out a chain to a file. fmt is the formatting to be used for values.
        """
        optr = open(fname,'w')
        optr.write('# Written by Chain.dump\n')
        optr.write('#\n')
        optr.write('## ' + ' '.join(self.names) + '\n')
        optr.write('#\n')
        optr.write('# Initial model:\n#\n')
        for k,v in self.model.iteritems():
            if isinstance(v, list):
                optr.write('# ' + k + ' = ' + ' '.join([str(x) for x in v]) + '\n')
            else:
                optr.write('# ' + k + ' = ' + str(v) + '\n')
        optr.write('#\n')
        optr.write('# nstore  = ' + str(self.nstore) + '\n')
        optr.write('# method  = ' + self.method + '\n')
        optr.write('# nwalker = ' + str(self.nwalker) + '\n')
        optr.write('# stretch = ' + str(self.stretch) + '\n')
        optr.write('# Minimum chi**2 = ' + str(self.chmin) + '\n#\n')
        optr.write('# Jump start:\n#\n')
        if self.jump is not None:
            self.jump.tofobj(optr, '# ')
        optr.write('#\n# Jump end\n#\n')
        if self.vals is not None:
            np.savetxt(optr, self.vals, fmt)
        optr.close()

    def range(self, indx):
        """
        Returns range (min,max) of the parameter of index indx
        """
        return (self.vals[:,indx].min(),self.vals[:,indx].max())

    def twod(self, xind, x1, x2, nx, yind, y1, y2, ny, xb=0., yb=0., frac=None):
        """
        Returns a 2D image representing the distribution of two parameters
        of a Chain by accumulating them and blurring them. Can also compute
        contour levels marking particular fractions of points, e.g. 0.68, 0.95
        
        xind, x1, x2, nx -- index, range and number of pixels for the X axis. 
                            Leftmost pixel starts at x1, rightmost ends at x2.
        yind, y1, y2, ny -- index, range and number of pixels for the Y axis
        xb, yb           -- amount of blurring in X and Y (RMS in pixels).
        frac             -- fractions to determine equivalent contour levels.
        """

        # get values
        x = self.vals[:,xind]
        y = self.vals[:,yind]

        return subs.pts2cont(x, y, x1, x2, nx, y1, y2, ny, xb, yb, frac)

    def get(self, vname):
        """
        Returns numpy 1D array with values of the parameter called vname. In incorrect name will cause
        an error
        """
        return self.vals[:,self.names.index(vname)]

    def vpars(self):
        """
        Returns list of variable parameters, i.e. excluding dependent extras such as 'chisq'
        """
        return self.names[:self.names.index('chisq')]

    def nvpars(self):
        """
        Returns number of variable parameters, i.e. excluding dependent extras such as 'chisq'
        """
        return self.names.index('chisq')

    def chisq(self):
        """
        Returns array of chisq values. Just a convenience interface to 'get'
        """
        return self.get('chisq')

    def __len__(self):
        """
        Returns the length of the Chain purely in terms of the number of entries
        stored in memory.
        """
        if self.vals is not None:
            return len(self.vals)
        else:
            return 0

class Fchain(Chain):
    """
    A Chain object linked to a file to allow a continuous writing to disk.
    """
    def __init__(self, *args):
        """
        Initialises a Chain.

        One argument:

        fname  -- a file name containing a previously stored run

        Two arguments:

        fname  -- a file name containing a previously stored run
        nskip  -- number of point groups to skip when inputting old file
 
        Nine arguments:

        model   -- dictionary defining all parameters and initial values. These can include both
                  variable and non-variable parameters and are simply included for reference
        nstore  -- how often models have been written (15 = 1 in every 15)
        method  -- 'm' for Metropolis-Hastings jumps using 'jump' or 'a' for the affine method.
        jump    -- the jump distribution
        nwalker -- number of "walkers" if using affine method.
        stretch -- stretch factor for affine method
        extras  -- extra parameters in addition to the variable parameters contained in jump
        cmin    -- minimum chi**2. This used to scale the uncertainties.
        fname   -- name of file to associate with this Chain
        """

        if len(args) == 1:
            fname = args[0] 
            Chain.__init__(self, args[0])
        elif len(args) == 2:
            fname = args[0] 
            Chain.__init__(self, args[0], args[1])
        elif len(args) == 9:
            Chain.__init__(self, args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7])
            fname = args[8]
            self.tofile(fname)
        else:
            raise MCMCError('Chain.Fchain: invalid number of arguments')
        self.fp = open(fname, 'a')

    def add_line(self, values, fmt='%.16e', nomem=False):
        """
        Add a line and writes to disk
        """
        if not nomem:
            Chain.add_line(self, values)
        for val in values:
            self.fp.write((fmt + ' ') % (val,))
        self.fp.write('\n')
        self.fp.flush()

    def add_info(self, message):
        """
        Adds an information line
        """
        self.fp.write(message)
        self.fp.flush()

    def close(self):
        self.fp.close()

def gchain(chain):
    """
    generator of a chain. Allows iteration. Each call to it yields the
    next row of values as a dictionary
    """
    for i, vals in enumerate(chain.vals):
        yield i, dict(zip(chain.names, vals))

class IChain(object):
    
    """
    Represents the results of an MCMC run for input from a file as an iterator. This allows line-by-line
    access without using tons of memory for a large file. 'IChain' stands for 'Input Chain'. Every next
    line returns an entry number (starting at 0) and the corresponding values

    Attributes:

    names   -- list of parameter names, possibly including extras such as chisq, lp, wdwarf. The first of the
               extras should always be 'chisq' to allow identification of the variable parameters.
    jump    -- the jump distribution used if the Metroplis-Hastings method is being used (a Trial object)
    nstore  -- how often data is stored.
    chmin   -- minimum chi**2, used if the errors need re-scaling.
    model   -- dictionary defining model
    method  -- the method being used
    nwalker -- the number of walkers for the affine method
    stretch -- the stretch factor for the affine method.
    nline   -- entry number we are at the start of. Next call will return the values corresponding to these
               starts from 0
    """
    
    def __init__(self, fname):
        """
        Initialises an IChain. 

        Argument:

        fname  -- a file name containing a previously stored run
        """

        # Read headers
        self.fobj,self.names,self.model,self.nstore,self.nwalker,self.stretch,\
            self.method,self.chmin,self.jump = _chain_from_file(fname)
        self.nline = 0

    def __iter__(self):
        """
        Start of loop return to start of data. Go to start
        of file to ensure that we are positioned correctly.
        """
        self.fobj.seek(0)
        self.nline = 0
        return self

    def next(self):
        """
        Return next entry number and values as a numpy 1D array
        """
        self.nline += 1
        while True:
            line = self.fobj.readline()
            if line == '':
                raise StopIteration()
            elif not line.startswith('#'):
                break

        vals = line.split()
        if len(vals) == len(self.names):
            return (self.nline-1,np.cast[float](vals))
        else:
            raise Exception('Line = ' + line + ' has ' +  str(len(vals)) + ' values compared to expected number = ' + str(len(self.names)))

def merge_chains(chains):
    """
    Merges multiple chains contained in the input sequence together.
    Returns the combined sequence.
    """
    if len(chains) == 0:
        raise MCMCError('merge_chains: no Chains supplied.')
    if not isinstance(chains[0],Chain):
        raise MCMCError('merge_chains: first element of sequence is not a Chain.')

    tchain = Chain(chains[0].names)
    for chain in chains:
        if not isinstance(chain,Chain):
            raise MCMCError('merge_chains: encountered sequence elelemt that is not a Chain.')
        
        if chain.names != tchain.names:
            raise MCMCError('merge_chains: names conflict: ' + str(chain.names) + ' vs ' + str(tchain.names))
        
        if tchain.vals is None:
            tchain.vals = chain.vals.copy()
        else:
            tchain.vals = np.concatenate((tchain.vals, chain.vals))

    return tchain

class Lcurve(subs.Odict):
    
    """
    Class representing an Lcurve model. Initialised by reading
    in an Lcurve model file.
    """

    def __init__(self, model, prior=None):
        """
        model  -- the name of the model file to initialise this
        prior  -- a function to return -2*log(prior prob) to add to the
                  chi**2 during trials. This is a function that will be passed
                  a dictionary of parameter name/value pairs. Here is an example
                  applying a prior constraint on a parameter 'q':
                  def prior(pdict):
                     return ((pdict['q']-0.1)/0.05)**2
                  which imposes a gaussian prior mean 0.1, RMS=0.05
                  Note that the prior can be relative to its peak value so no normalisation
                  terms are needed.
        args   -- a tuple of extra arguments to pass through to 'prior' if necessary. These
                  come after name/value dictionary argument.
        """
        super(subs.Odict, self).__init__(self)
        fptr = open(model)
        for line in fptr:
            eq = line.find('=') 
            if eq > -1:
                self[line[:eq].strip()] = line[eq+1:].strip().split()
        fptr.close()
        self.prior = prior

    def set(self, p):
        """
        Sets the model to a dictionary of parameter values
        """
        for key,value in p.iteritems():
            if key in self:
                self[key][0] = value
            else:
                raise Exception('Lcurve.set: could not find parameter = ' + key + ' in current model.')

    def vcheck(self, name, vmin, vmax):
        """Checks that a parameter of a given name is in range"""
        v = float(self[name][0])
        return (v >= vmin and v <= vmax)

    def ok(self):
        """
        Checks for silly parameter values. Runnning this before chisq could
        save time and disaster
        """
        return self.vcheck('q', 0.001, 100.) and self.vcheck('iangle', 0., 90.) and \
            self.vcheck('r1',0.,1.) and (float(self['r2'][0]) <= 0. or \
                                             self.vcheck('r2',0.,1-roche.xl1(self['q'][0]))) and \
            self.vcheck('cphi3',0.,0.25) and self.vcheck('cphi4',float(self['cphi3'][0]),0.25) and \
            self.vcheck('t1',0.,1.e6) and self.vcheck('t2',-1.e6,1.e6) and \
            self.vcheck('ldc1_1',-20.,20.) and self.vcheck('ldc1_2',-20.,20.) and \
            self.vcheck('ldc1_3',-20.,20.) and self.vcheck('ldc1_4',-20.,20.) and \
            self.vcheck('ldc2_1',-20.,20.) and self.vcheck('ldc2_2',-20.,20.) and \
            self.vcheck('ldc2_3',-20.,20.) and self.vcheck('ldc2_3',-20.,20.) and \
            self.vcheck('period',0.,100.) and self.vcheck('gravity_dark1',0.,1.) and \
            self.vcheck('gravity_dark2',0.,1.) and \
            (self['add_disc'][0] == '0' or \
                 ((float(self['rdisc1'][0]) <= 0. or self.vcheck('rdisc1',float(self['r1'][0]),1.0)) \
                      and self.vcheck('rdisc2',float(self['rdisc1'][0]),1.) and self.vcheck('temp_disc',0.,1.e6) \
                      and self.vcheck('texp_disc',-10.,10.))) and \
                      (self['add_spot'][0] == '0' or (self.vcheck('radius_spot',0.,1.) and \
                                                       self.vcheck('length_spot',0.,10.) and \
                                                       self.vcheck('expon_spot',0.,30.) and self.vcheck('temp_spot',0.,1.e6) and \
                                                       self.vcheck('cfrac_spot',0.,1.) and self.vcheck('epow_spot',0.,10.)))
               
    def write(self,fname=None):
        """
        Writes model to a file, or a temporary file if fname=None. It returns the name of the file.
        You should delete the temporary file at some point.
        """
        
        import tempfile
        if not fname:
            (fd,fname) = tempfile.mkstemp()
            fptr = os.fdopen(fd, 'w')
        else:
            fptr = open(fname,'w')
                
        for (key,value) in self.iteritems():
            fptr.write('%-20s =' % key)
            for item in value:
                fptr.write(' ' + str(item))
            fptr.write('\n')
        fptr.close()
        return fname

    def lnprior(self):
        """
        Returns ln(prior prob) for current light curve model
        which contains a prior probability function
        """

        if self.prior is None:
            lp = 0.
        else:
            vdict = {}
            for (key,val) in self.iteritems():
                if len(val) > 1:
                    vdict[key] = float(val[0])
            lp = self.prior(vdict)
        return lp

    def chisqlp(self, data):
        """
        Computes chi**2 of the current model against a given 
        data file. This first writes out the current model to 
        a temporary disc file which it then passes to a call 
        to lroche along with the name of the data file. It returns 
        (chisq,wnok,wdwarf) where chisq and wnok are the weighted 
        chi**2 and number of data points while wdwarf
        is the white dwarf's contribution as reported by lroche.
        If something goes wrong, it comes back with 'None'
        in each of these values.
        """

        data = str(data)
        if not os.path.isfile(data):
            raise MCMCError(data + ' does not appear to be a file')

        # write the model to a temporary file
        tfile = self.write()

        # build the command, run it, read the results.
        args = (lroche, 'model=' + tfile, 'data=' + data, 'noise=0','scale=yes', 'seed=12345', 'nfile=0', 'device=null', 'nodefs')
        sout, serr = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        output = sout.split('\n')

        # delete the temporary file
        os.remove(tfile)

        # interpret the output
        subout = [out for out in output if out.startswith('Weighted')]

        if len(subout) == 1:
            eq    = subout[0].find('=')
            comma = subout[0].find(',')
            chisq = float(subout[0][eq+2:comma])
            wnok  = float(subout[0][subout[0].find('wnok =')+6:].strip())
            subout = [out for out in output if out.startswith('White dwarf')]
            if len(subout) == 1:
                eq     = subout[0].find('=')
                wdwarf = float(subout[0][eq+2:])
            else:
                print "Can't find white dwarf's contribution. lroche out of date"
                chisq = wnok = lp = wdwarf = None

        else:
            chisq = wnok = wdwarf = None
            raise MCMCError('Output from lroche failure: ' + sout + '\n' + serr)


        return (chisq,wnok,wdwarf)

    def var(self):
        """
        Returns dictionary of current values of variable parameters
        keyed on the parameter name
        """
        vpar = {}
        use_radii = (self['use_radii'][0] == '1')
        for (key,value) in self.iteritems():
            if len(value) > 1 and value[3] == '1':
                if (use_radii and (key == 'r1' or key == 'r2')) or \
                        (not use_radii and (key == 'cphi3' or key == 'cphi4')) or \
                        (key != 'r1' and key != 'r2' and key != 'cphi3' and key != 'cphi4'):
                    vpar[key] = float(value[0])

        return vpar

    def vvar(self):
        """
        Returns names of variable parameters
        """
        vnam = []
        use_radii = (self['use_radii'][0] == '1')
        for (key,value) in self.iteritems():
            if len(value) > 1 and value[3] == '1':
                if (use_radii and (key == 'r1' or key == 'r2')) or \
                        (not use_radii and (key == 'cphi3' or key == 'cphi4')) or \
                        (key != 'r1' and key != 'r2' and key != 'cphi3' and key != 'cphi4'):
                    vnam.append(key)

        return vnam

    def nvar(self):
        """
        Returns the number of variable parameters
        """
        nvar = 0
        use_radii = (self['use_radii'][0] == '1')
        for (key,value) in self.iteritems():
            if len(value) > 1 and value[3] == '1':
                if (use_radii and (key == 'r1' or key == 'r2')) or \
                        (not use_radii and (key == 'cphi3' or key == 'cphi4')) or \
                        (key != 'r1' and key != 'r2' and key != 'cphi3' and key != 'cphi4'):
                    nvar += 1
        return nvar
        
    def trial(self, trl, data, cold, lpold, cmin, tmodel, wdold, plarge=0., slarge=3.):
        """
        This tries out an MCMC jump. Returns (new,cnew,lp,wdwarf)
        new = flag indicating whether the jump was accepted,
        cnew is the new chi**2 (= cold if the jump failed), lp and 
        wdwarf are the prior probability and white dwarf contributions
        of the new model. The Lcurve is modified if the trial was a 
        success, otherwise it is unchanged. Note that a trial that 
        jumps into a forbidden region of parameter space is treated 
        exactly as if it were a failed trial.

        Arguments:

        trl    --- Jump object to define covariances
        data   --- data file (sent to lroche)
        cold   --- previous chi**2
        lpold  --- previous -2*log(prior prob)
        cmin   --- minimum value of chi**2: this is a number to enable scaling
                   of errors to get chi**2/dof = 1 and should remain fixed.
        tmodel --- file to store model prior to running 'lroche'
        wdold  --- previous value of white dwarf's contribution
        plarge --- probability of taking an extra large jump
        slarge -- scale factor of extra large jumps (>1)
        """

        data = str(data)
        vold = self.var()
        if np.random.uniform() < plarge:
            ltrl = slarge*trl
            vnew = ltrl.trial(vold)
        else:
            vnew = trl.trial(vold)

        # modify variable values
        for (key,value) in vnew.iteritems():
            self[key][0] = value

        new = False
        if self.ok():
            # save model and then compute chi**2
            if tmodel:
                self.write(tmodel)
            (cnew,lp,wnok,wdwarf) = self.chisqlp(data)
            if cnew is not None:
                # Here is the crucial part of the MCMC method.
                # diff = -2*log(posterior prob new/posterior prob old). 
                # Jump immediately if this is  < 0, i.e. new is more probable
                # If diff > 0, we jump on a probablistic basis. wnok and cmin
                # are used here to scale the chi**2
                diff = (cnew-cold)*wnok/cmin + (lp-lpold)
                if diff <= 0. or np.random.uniform() < m.exp(-diff/2.):
                    # success
                    new = True

        if not new:
            # Failed, go back to pre-trial values
            for (key,value) in vold.iteritems():
                self[key][0] = value
            cnew   = cold
            lp     = lpold
            wdwarf = wdold
 
        return (new,cnew,lp,wdwarf)
    
class Jump(object):
    """
    Object to store covariances used to generate trial values of variables. 

    Attributes:

    par    -- parameter names.
    sig    -- standard deviations, 1D numpy array
    corr   -- correlation coefficients, 2D numpy array
    """

    def __init__(self, *args):
        """
        Initialises a Jump object.

        args   -- If one argument, then it can be the name of a file with covariance 
                  information or a file object set at the start of such information
                  as written to a Chain log. It can cope with output from levmarq and
                  genlmq.py (which are slightly different from each other).

                  If two arguments: parameter names and standard deviations, (corrleation matrix
                  will be diagonal)

                  If three arguments: parameter names, standard deviations and correlation matrix
        """

        if len(args) == 1 and (isinstance(args[0], str) or isinstance(args[0], file)):
            # Read from an lmq file (specified by a name string) or a file object 
            # opened at the right place. In latter case it is expected that each line
            # starts with '# '
            if isinstance(args[0], str):
                fptr = open(args[0])
                whole = fptr.readlines()
                fptr.close()
            elif isinstance(args[0], file):
                whole = []
                for line in args[0]:
                    if line.startswith('# Jump end'): break
                    whole.append(line[2:])

            # Two possibilities: output from levmarq which starts with
            # 'Parameter values' or from genlmq.py which starts with 
            # 'RMS values'
            pdict = {}
            rread = False
            pread = False
            self.par = []
            self.sig = []
            nvar = 0
            for line in whole:
                if line.startswith('RMS values'):
                    rread = True
                elif line.startswith('Parameter values'):
                    pread = True
                elif line.startswith('Correlation'):
                    break
                elif rread and line.find('=') > -1:
                    eq = line.find('=')
                    self.par.append(line[:eq].strip())
                    self.sig.append(float(line[eq+1:].strip()))
                    pdict[self.par[-1]] = nvar 
                    nvar += 1
                elif pread and line.find('=') > -1:
                    eq = line.find('=')
                    self.par.append(line[:eq].strip())
                    porm = line.find('+/-')
                    if porm < 0:
                        raise MCMCError('Jump.__init__: failed to find +/- for parameter ' + self.par[-1])
                    self.sig.append(float(line[porm+3:].strip()))
                    pdict[self.par[-1]] = nvar 
                    nvar += 1

            if not rread and not pread:
                raise MCMCError('Jump.__init__: failed to find start of RMS or Parameter')

            if len(self.par) == 0:
                raise MCMCError('Jump.__init__: failed to read any RMS values')

            self.sig  = np.array(self.sig)
            self.corr = np.diag(np.ones_like(self.sig))

            for line in whole:
                if line.startswith('r('):
                    br = line.find(')')
                    (p1,p2) = line[2:br].split(',')
                    r = float(line[line.find('=')+1:].strip())
                    indx1 = pdict[p1]
                    indx2 = pdict[p2]
                    self.corr[indx1][indx2] = r
                    self.corr[indx2][indx1] = r

        elif len(args) == 2:

            # parameter names and standard deviations
            self.par  = list(args[0])
            self.sig  = np.array(args[1])
            self.corr = np.diag(np.ones_like(self.sig))

        elif len(args) == 3:

            # parameter names, standard deviations and correlations
            self.par  = list(args[0])
            self.sig  = np.array(args[1])
            if len(self.par) != len(self.sig):
                raise MCMCError('Jump.__init__: differing numbers of parameters and sigma')
            self.corr = np.array(args[2])

        else:
            raise MCMCError('Jump.__init__: unrecognised number and/or types of argumenst')
        
        if not (np.linalg.eig(self.corr)[0] > 0.).all():
            print 'WARNING: Covariance matrix appears not to be positive definite:'
            print 'Determinant             =',np.linalg.det(self.corr)
            print 'Eigenvalues             =',np.linalg.eig(self.corr)[0]
            print 'Product of eigenvalues =',np.prod(np.linalg.eig(self.corr)[0])
            print 'This will cause problems with the mcmc jump distribution.'

    def __str__(self):
        """
        Straightforward output of contents
        """
        return str(self.par) + '\n' + str(self.sig) + '\n' + str(self.corr)

    def trim(self, vars):
        """
        Trim the Jump to only cover the names listed in vars
        """
        # extract list of indices
        pdict = dict(zip(self.par, range(len(self.par))))
        cols  = [pdict[v] for v in vars]
        rows  = [[nc] for nc in cols]
        
        # now trim
        self.par  = list(vars)
        self.sig  = self.sig[cols]
        self.corr = self.corr[cols, rows]

    def trial(self, pd):
        """
        Generate new model with random jump according to the Jump parameters.

        pd  -- current model dictionary (not all of which will be altered in general)
        
        Returns new model dictionary
        """
        td = dict(pd)
        delta = self.sig*np.random.multivariate_normal(np.zeros_like(self.sig), self.corr)
        for k,v in zip(self.par, delta):
            td[k] += v
        return td

    def check_par(self, params):
        """
        Checks that a list of parameter names is consistent with
        the parameters in 'Jump'

        params -- dictionary where key names are the parameters to be checked
        """
        local    = set(self.par.keys())
        external = set(params.keys())
        return local == external

    def tofile(self, fname, means=None):
        """
        Writes out to a Jump object to a file

        fname  -- file to write to
        means  -- mean values of parameters (will produce
                  genuine levmarq format)
        """
        fp = open(fname, 'w')
        self.tofobj(fp, start='', means=means)
        fp.close()

    def tofobj(self, fo, start='', means=None):
        """
        Writes out to a file object in lmq format, with each line starting with 'start'

        fo    -- file object, opened for write access
        start -- string appended to start of each line.
        means  -- mean values of parameters (will produce
                  genuine levmarq format)
        """

        if means is None:
            fo.write(start + 'RMS values:\n' + start + '\n')
            for p,s in zip(self.par, self.sig):
                fo.write(start + p + ' = ' + str(s) + '\n')
        else:
            if len(means) != len(self.sig):
                raise MCMCError('Jump.tofobj: means and sigmas have differing lengths')

            fo.write(start + 'Parameter values:\n' + start + '\n')
            for p,m,s in zip(self.par, means, self.sig):
                fo.write(start + p + ' = ' + str(m) + ' +/- ' + str(s) + '\n')

        fo.write(start + '\n')
        fo.write(start + 'Correlation coefficients:\n' + start + '\n')
        for i,pi in enumerate(self.par[:-1]):
            for j,pj in enumerate(self.par[i+1:]):
                fo.write(start + 'r(' + pi + ',' + pj + ') = ' + str(self.corr[i,j+i+1]) + '\n')

    def cov(self):
        """
        Returns covariance
        """
        return self.sig*((self.sig*np.array(self.corr)).transpose())

    def __neq__(self, other):
        """
        Tests for inequality between two jump objects which means different parameters or covariances
        """
        if self.par != other.par or self.sig != other.sig or self.corr != other.corr:
            return True
        else:
            return False

    def __eq__(self, other):
        """
        Tests for equality between two jump objects which means same parameters and covariances
        """
        return not self.__neq__(other)

    def __imul__(self, scale):
        self.sig *= scale
        return self

    def __mul__(self, scale):
        temp  = Jump(self.par, self.sig, self.corr)
        temp *= scale
        return temp

    def __rmul__(self, scale):
        temp  = Jump(self.par, self.sig, self.corr)
        temp *= scale
        return temp

class MCMCError(Exception):
    pass
