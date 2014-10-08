#!/usr/bin/env python

"""
a module to make SIMBAD requests


Classes
=======

Query  -- defines and makes queries.

"""

import exceptions, urllib

class Query(object):

    def __init__(self, Target = None, RA = None, Dec = None, \
                     Radius = 1., url = 'http://simbad.u-strasbg.fr/simbad'):
        """
        Defines a SIMBAD query.

        Target   -- Target name. Overrides any RA and Dec
        RA, Dec  -- Position. If strings then they will be assumed to be HMS format,
                    if decimal they will be taken to be both in degrees. ICRS assumed.
        Radius   -- search radius for coordinate queries, arcminutes
        url      -- base url to query
        """

        if Target is None and (RA is None or Dec is None):
            raise SimbadError('You must specify either a target or a position')

        self.Target = Target
        self.RA     = RA
        self.Dec    = Dec
        self.Radius = Radius
        if url.endswith('/'):
            self.url = url
        else:
            self.url = url + '/'

    def query(self, maxobj=100):
        """
        Makes a query. Returns a list of object name, position and reference frame (currently
        always 'ICRS'. This will have multiple objects if there are multiple matches to your query.
        'maxobj' is the maximum number to return.
        """
        if self.Target:
            q = 'set limit ' + str(maxobj) + '\nformat object form1 "Target: %IDLIST(1) | %COO(A D;ICRS)"\nquery ' + self.Target
        else:
            if isinstance(self.RA, str):
                coord = self.RA
            else:
                coord = _d2hms(self.RA, sep=' ')
            if isinstance(self.Dec, str):
                coord += self.Dec
            else:
                coord += _d2hms(self.Dec, sep=' ', sign=True)
            q = 'set limit ' + str(maxobj) + '\nformat object form1 "Target: %IDLIST(1) | %COO(A D;ICRS)"\nquery ' + coord + ' radius=' + \
                str(self.Radius) + 'm' + ' frame=ICRS'

        uptr = urllib.urlopen(self.url + 'sim-script?' + urllib.urlencode({'submit' : 'submit script', 'script' : q}))
        data  = False
        error = False
        results = []
        for line in uptr:
            if line.startswith('::data::'):
               data = True
            if line.startswith('::error::'):
               error = True
            if data:
                if line.startswith('Target:'):
                    (name,coords) = line[7:].split(' | ')
                    results.append({'Name' : name.strip(), 'Position' : coords.strip(), 'Frame' : 'ICRS'})
        uptr.close()
        if error and len(results):
            print 'trm.simbad.Query.query: there appear to be some results but an error was unexpectedly raised.'
        return results

def _d2hms(hour, precision=2, sep=':', dp=3, sign=False):
    """hms = d2hms(hour, precision=2) -- returns HH MM SS.SS time string from a decimal hour

    hour      -- decimal hour
    precision -- 0 gives hh:mm, 1 gives hh:mm:ss, 2 gives hh:mm:ss.sss
    sep       -- separator
    dp        -- number of decimal places in the case of precision = 2
    sign      -- True if you want to make sure a '+' is appended at the start
    """
    h = int(abs(hour))
    if precision == 0:
        m = int(60.*(abs(hour)-h)+0.5)
        st = '%02.2d%s%02.2d' % (h,sep,m)
    else:
        m = int(60.*(abs(hour)-h))
    if precision == 1:
        s = int(3600.*(abs(hour)-h-m/60.)+0.5)
        st = '%02.2d%s%02.2d%s%02.2d' % (h,sep,m,sep,s)
    else:
        s = 3600.*(abs(hour)-h-m/60.)
        st = ('%02.2d%s%02.2d%s%0' + str(dp+3) + '.' + str(dp) +'f') % (h,sep,m,sep,s)
    
    if sign:
        if hour < 0:
            st = '-' + st
        else:
            st = '+' + st
    return st

# Exception class
class SimbadError(exceptions.Exception):
    pass


if __name__ == '__main__':
    q = Query(RA = '01 23 34.5', Dec = '+12 34 50.5', Radius=10)
    print q.query()


