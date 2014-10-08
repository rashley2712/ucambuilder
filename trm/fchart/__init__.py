
"""
Module to make finding charts
"""

import os
import math as m
from trm import subs, sla
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.stats as stats
import pyfits
import astLib.astWCS as astWCS
import urllib
import urllib2

# Surveys to get from DSS
SURVEYS = {'POSS2RED'  : 'poss2ukstu_red',
           'POSS2BLUE' : 'poss2ukstu_blue',
           'POSS2IR'   : 'poss2ukstu_ir',
           'POSS1RED'  : 'poss1_red',
           'POSS1BLUE' : 'poss1_blue',
           'QUICKV'    : 'quickv'}
# colours
RED  = '#aa0000'
BLUE = '#0000aa'

def get_from_dss(ra, dec, size, fname=None, prec=1, survey='POSS2BLUE', clobber=True):
    """
    Gets a FITS file from the DSS

    ra    -- ra (string hh mm ss.ss)
    dec   -- dec (string +/-dd mm ss.ss)
    fname -- name of output (including .fits.gz). If not specified it will be
             constructed from the RA and Dec.
    prec  -- used to construct fname if not specified. See precision parameter of
             trm.subs.d2hms to understand this.
    size  -- field size (arcmin)

    Returns (fname, position)

    fname    -- the name of the FITS file created.
    """
    if not survey in SURVEYS:
        raise Exception('get_from_dss: survey name = ' + survey + ' is invalid.')

    rad  = subs.hms2d(ra)
    decd = subs.hms2d(dec)
    if rad < 0. or rad >= 24 or decd < -90. or decd > +90.:
        raise Exception('get_from_dss: position out of range ra,dec = ' + ra + ' ' + dec)

    if fname is None:
        ras   = ra[:ra.find('.')].replace(' ','')
        decs  = dec[:dec.find('.')].replace(' ','')
        fname = ras + decs + '.fits.gz'
    elif not fname.endswith('.fits.gz'):
        raise Exception('get_from_dss: filename = ' + fname + ' does not end with .fits.gz')

    if clobber or not os.path.exists(fname):

        # URL at STScI
        url = 'http://stdatu.stsci.edu/cgi-bin/dss_search'

        # encode the request string.
        data = urllib.urlencode({'v' : SURVEYS[survey], 'r': ra, 'd' : dec, 'e' : 'J2000', \
                                     'h' : size, 'w' : size, 'f' : 'fits', \
                                     'c' : 'gz', 's' : 'on'})

        # get FITS data
        results = urllib2.urlopen(url, data)
        with open(fname, 'w') as f:
            f.write(results.read())

    return fname

def make_eso_chart(fname, cname, target, info, pid, pi, pos, field, scale, stype='sec', \
                       source=None, angle=None, mark=True, s1=0.04, s2=0.1, plo=50.0, phi=99.8, \
                       sloc=-0.43, pwidth=16, aspect=0.6, pm=None, obsdate=None):
    """
    Make a chart suitable for VLT. Arguments:

    fname   -- fits file from DSS, e.g. 'PG1018-047.fits'
    cname   -- chart name (e.g. "name.pdf". The type of file is deduced from the suffix)
    target  -- target name, e.g. 'PG1018-047'
    info    -- extra information, e.g. 'Slit angle: parallactic'. A list of strings will
               be printed line by line
    pid     -- programme ID, e.g. '088.D-0041(A)'
    pi      -- PI name, e.g. 'Marsh'
    pos     -- position, e.g. '10 21 10.6 -04 56 19.6' (assumed 2000)
    field   -- field width to display, arcmin, e.g. 1.0
    scale   -- size of scale to indicate in arcsec, e.g. 30
    stype   -- 'sec' or 'min' to indicate how to mark the scale
    source  -- source of fits file, e.g. 'DSS blue'. Ignore for automatic version.
    angle   -- indicator / slit angle to display, degrees. None to ignore.
    mark    -- True to indicate object
    s1      -- fraction of field for start of slit and object markers (displaced from object)
    s2      -- fraction of field for end of slit and object markers (displaced from object)
    plo     -- Low image display level, percentile 
    phi     -- High image display level, percentile
    sloc    -- vertical location of scale bar in terms of 'field'
    pwidth  -- (approx) plot width in inches
    aspect  -- aspect (vertical/horizontal) [should be < 1]
    pm      -- tuple of proper motion in RA, and dec, arcsec/year in both coords. This is used to correct
               the predicted position of the target in the image by finding the date on which the image 
               was taken.
    obsdate -- if you set pm, then if you also set date, an arrow will be drawn from target to its 
               expected position at obsdate. This should be is YYYY-MM-DD format. 
    Returns (ilo,ihi) display levels used
    """

    if cname.endswith('.pdf'):
        form = 'pdf'
    elif cname.endswith('.jpg'):
        form = 'jpg'
    else:
        raise Exception('make_eso_chart: do not recognize the format corresponding to ' + cname)

    if stype != 'sec' and stype != 'min':
        raise Exception('make_eso_chart: stype must = "sec" or "min"')

    # ra, dec in degrees
    ra,dec,sys = subs.str2radec(pos)
    ra *= 15

    # Read data
    hdulist   = pyfits.open(fname)
    data      = hdulist[0].data
    head      = hdulist[0].header
    hdulist.close()

    arrow = False
    if pm:
        # try to correct position
        if 'DATE-OBS' in head:
            dobs    = head['DATE-OBS']
            year    = int(dobs[:4])
            month   = int(dobs[5:7])
            day     = int(dobs[8:10])
            deltat  = (sla.cldj(year,month,day) - sla.cldj(2000,1,1))/365.25
            dra     = pm[0]*deltat/np.cos(np.radians(dec))/3600.
            ddec    = pm[1]*deltat/3600.
            ra  += dra
            dec += ddec
            if obsdate:
                yearp,monthp,dayp = obsdate.split('-')
                deltat  = (sla.cldj(int(yearp),int(monthp),int(dayp)) - sla.cldj(year,month,day))/365.25
                arrow   = True
                dra     = pm[0]*deltat/np.cos(np.radians(dec))/3600.
                ddec    = pm[1]*deltat/3600.
        else:
            print 'Could not find DATE-OBS in header'

    # Read WCS info
    wcs = astWCS.WCS(fname)
    dx  = 60.*wcs.getXPixelSizeDeg()
    dy  = 60.*wcs.getYPixelSizeDeg()
    rot = wcs.getRotationDeg()
    if rot > 180.:
        rot -= 360.
    if not wcs.coordsAreInImage(ra,dec):
        print 'Warning: coordinates not inside image'

    x,y = wcs.wcs2pix(ra, dec)
    if arrow:
        if not wcs.coordsAreInImage(ra+dra,dec+ddec):
            print 'Warning: end of proper motion arrow is not inside image'
        xa,ya = wcs.wcs2pix(ra+dra, dec+ddec)

    ny, nx = data.shape

    # plot limits
    limits = (dx*(0.5-x),dx*(nx+0.5-x),dy*(0.5-y),dy*(ny+0.5-y))
 
    if source is None:
        if 'SURVEY' in head:
            source = head['SURVEY']
            if source == 'POSSII-F':
                source = 'POSS-II red'
            elif source == 'POSSII-J':
                source = 'POSS-II blue'
            elif source == 'POSSII-N':
                source = 'POSS-II ir'
            elif source == 'POSSI-E':
                source = 'POSS-I red'
            elif source == 'POSSI-O':
                source = 'POSS-I blue'
            elif source == 'SERC-I':
                source = 'SERC ir'
            elif source == 'SERC-J':
                source = 'SERC blue'
            elif source == 'AAO-SES':
                source = 'AAO red'
            elif source == 'AAO-GR':
                source = 'AAO red'
        else:
            source = 'UNKNOWN'

    # Start plotting
    fig  = plt.figure(figsize=(pwidth,pwidth))
    axes = fig.add_subplot(111,aspect='equal')
    mpl.rc('font', size=24)

    # Derive the plot range
    ilo  = stats.scoreatpercentile(data.flat, plo)
    ihi  = stats.scoreatpercentile(data.flat, phi)
    
    # Plot
    plt.imshow(data, vmin=ilo, vmax=ihi, cmap=cm.binary, extent=limits, interpolation='nearest')
    axes.autoscale(False)

    # draw slit angle
    if angle is not None:
        ca = m.cos(m.radians(angle+rot))
        sa = m.sin(m.radians(angle+rot))
        t1 = s1*field
        t2 = s2*field
        plt.plot([-sa*t1,-sa*t2],[+ca*t1,+ca*t2],BLUE,lw=3)
        plt.plot([+sa*t1,+sa*t2],[-ca*t1,-ca*t2],BLUE,lw=3)

    # Mark object
    if mark:
        t1 = s1*field
        t2 = s2*field
        plt.plot([t1,t2],[0,0],RED,lw=3)
        plt.plot([0,0],[t1,t2],RED,lw=3)

    # draw arrow
    if arrow:
        plt.arrow( 0, 0, dx*(xa-x), dy*(ya-y))

    # Draw scale bar
    if stype == 'sec':
        plt.text(0, (sloc-0.05)*field, str(scale) + ' arcsec', horizontalalignment='center',color=RED)
    else:
        plt.text(0, (sloc-0.05)*field, str(scale/60) + ' arcmin', horizontalalignment='center',color=RED)

    scale /= 60.
    plt.plot([-0.5*scale,+0.5*scale],[sloc*field,sloc*field],RED,lw=3)
    plt.plot([+0.5*scale,+0.5*scale],[(sloc+0.02)*field,(sloc-0.02)*field],RED,lw=3)
    plt.plot([-0.5*scale,-0.5*scale],[(sloc+0.02)*field,(sloc-0.02)*field],RED,lw=3)

    # North-East indicator. ev, nv = East / North vectors
    xc, yc = 0.4*field, -0.4*field
    rot   = m.radians(rot)
    rmat  = np.matrix(((m.cos(rot),-m.sin(rot)),(m.sin(rot),m.cos(rot))))
    
    nv = np.array((0.,0.2*field))
    nv = nv.reshape((2,1))
    nv = rmat*nv
    nv = np.array(nv).flatten()
    
    ev = np.array((-0.2*field,0.))
    ev = ev.reshape((2,1))
    ev = rmat*ev
    ev = np.array(ev).flatten()

    plt.plot([xc,xc+nv[0]],[yc,yc+nv[1]],RED,lw=3)
    plt.text(xc+1.1*nv[0], yc+1.1*nv[1], 'N', horizontalalignment='center', color=RED)
    plt.plot([xc,xc+ev[0]],[yc,yc+ev[1]],RED,lw=3)
    plt.text(xc+1.15*ev[0], yc+1.1*ev[1], 'E', verticalalignment='center', color=RED)

    # finally the textual info with a helper function
    def ptext(x, y, limits, tstr):
        """
        Converts relative (0-1) x,y limits into data coords and plots
        a string.
        """
        xd = limits[0]+x*(limits[1]-limits[0])
        yd = limits[2]+y*(limits[3]-limits[2])
        plt.text(xd, yd, tstr, horizontalalignment='left',color=BLUE)

    xoff = 0.02
    dely = 0.035
    yoff = 0.96
    
    ptext(xoff, yoff, limits, pid)

    yoff -= dely
    ptext(xoff, yoff, limits, 'PI: ' + pi)

    yoff -= 1.5*dely
    ptext(xoff, yoff, limits, 'Target: ' + target)

    yoff -= 1.5*dely
    rah,ram,ras,decd,decm,decs = pos.split()
    ptext(xoff, yoff, limits, 'RA (2000) : ' + rah + ' ' + ram + ' ' + ras)

    yoff -= dely
    ptext(xoff, yoff, limits, 'Dec (2000): ' + decd + ' ' + decm + ' ' + decs)

    yoff -= 1.5*dely
    ptext(xoff, yoff, limits, "Field: " + str(field) + "' x " + str(field) + "'")

    yoff -= dely
    ptext(xoff, yoff, limits, 'Date: ' + head['DATE-OBS'][:10])

    yoff -= dely
    ptext(xoff, yoff, limits, 'Survey: ' + source)    

    yoff -= 2.*dely
    if info is not None:
        if isinstance(info, list):
            for line in info:
                ptext(xoff, yoff, limits, line)
                yoff -= dely
        else:
            ptext(xoff, yoff, limits, info)

    plt.savefig(cname,format=form,bbox_inches='tight')
    return (ilo,ihi)

