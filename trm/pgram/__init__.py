#!/usr/bin/env python  
"""spectral analysis package


This package implements various functions to calculate periodograms. One set
of routines implements the Lomb-Scargle periodogram which is essentially the
result of fitting a sinusoid to some data. Some are implemented using Press &
Rybicki's method for converting from irregularly spaced to uniform grids that
can then be treated with FFTs. This can be fast for full frequency searches, but 
it can be also rather expensive in terms of memory.

Functions
========

fls   : Press & Rybicki's fast FFT-based Lomb-Scargle routine. 
fwls  : Weighted version of fls
fwmls : Weighted, floating mean version of fls
genf  : generates frequencies in same way as fls, fwls and fwmls
pdm   : phase-dispersion minimisation.
wmls  : Weighted, floating mean, multi-harmonic L-S, direct calculation

Example 
=======

The following produces a spectrum with a sharp peak at 2 cycles/unit x:

import numpy as np
import trm.pgram as pgram

# Fake some data
x = np.linspace(0.,10.,100)
f = 2.
y = np.sin(2.*np.pi*f*x)
e = 0.1*np.ones_like(x)

# Compute periodogram
ofac  = 2.
hifac = 1.
f,p = pgram.fwmls(x,y,e,ofac,hifac)

"""  

import numpy as np
import numpy.fft as fft
  
def __spread__(y, yy, n, x, m):  
  """ 
  Given an array yy(0:n-1), extirpolate (spread) a value y into 
  m actual array elements that best approximate the "fictional" 
  (i.e., possible noninteger) array element number x. The weights 
  used are coefficients of the Lagrange interpolating polynomial 
  Arguments: 
    y :  
    yy :  
    n :  
    x :  
    m :  
  Returns: 
     
  """  
  nfac=[0,1,1,2,6,24,120,720,5040,40320,362880]  
  if m > 10. :  
    print 'factorial table too small in spread'  
    return  
  
  ix=long(x)  
  if x == float(ix):   
    yy[ix]=yy[ix]+y  
  else:  
    ilo = long(x-0.5*float(m)+1.0)  
    ilo = min( max( ilo , 1 ), n-m+1 )   
    ihi = ilo+m-1  
    nden = nfac[m]  
    fac=x-ilo  
    for j in range(ilo+1,ihi+1): fac = fac*(x-j)  
    yy[ihi] = yy[ihi] + y*fac/(nden*(x-ihi))  
    for j in range(ihi-1,ilo-1,-1):  
      nden=(nden/(j+1-ilo))*(j-ihi)  
      yy[j] = yy[j] + y*fac/(nden*(x-j))  
  
def fls(x,y,ofac,hifac, MACC=4):  
  """ function fls
    Given abscissas x (which need not be equally spaced) and ordinates 
    y, and given a desired oversampling factor ofac (a typical value 
    being 4 or larger), this routine creates two arrays, the first with 
    a sequence of increasing frequencies up to hifac times the "average" 
    Nyquist frequency, the second with values of the Lomb-normalized 
    periodogram at those frequencies. The arrays x and y are not altered. 
    The mean is automatically subtracted, but this is not a "floating
    mean" periodogram routine.
   
    Arguments: 
       x     : abscissas array, (e.g. an array of times). 
       y     : ordinates array, (e.g. corresponding counts). 
       ofac  : oversampling factor. 
       hifac : hifac * "average" Nyquist frequency = highest frequency 
               for which values of the Lomb normalized periodogram will 
               be calculated. 
       MACC  : number of interpolation points per 1/4 cycle 
               of highest frequency        

    Returns: 
      f : An array of Lomb periodogram frequencies. 
      p : An array of corresponding values of the Lomb periodogram. 
 
    Reference:  
      Press, W. H. & Rybicki, G. B. 1989, ApJ vol. 338, p. 277-280. 
      Fast algorithm for spectral analysis of unevenly sampled data 
      (1989ApJ...338..277P) 

    History: 
      23/12/2011, v1.1, TRM (modified from version posted on astropy)
      23/02/2009, v1.0, MF 
        Translation of IDL code (orig. Numerical recipies) 
  """  

  # Check dimensions of input arrays  
  n = long(len(x))  
  if n != len(y):
    raise Exception('fls: incompatible arrays lengths')

  nout   = int(np.ceil(0.5*ofac*hifac*n))  

  nfreqt = long(ofac*hifac*n*MACC)   # Size the FFT as next power  
  nfreq  = 64L                       # of 2 above nfreqt.    
  while nfreq < nfreqt:   
    nfreq = 2*nfreq    
  ndim = long(2*nfreq)  
    
  # Compute the mean, variance  
  ave = y.mean() 
 
  # sample variance because the divisor is N-1  
  var = ((y-y.mean())**2).sum()/(len(y)-1)
   
  # and range of the data.  
  xmin = x.min()  
  xmax = x.max()  
  xdif = xmax-xmin  
  
  # extirpolate the data into the workspaces  
  wk1 = np.zeros(ndim, dtype='complex')  
  wk2 = np.zeros(ndim, dtype='complex')  
  
  fac   = ndim/(xdif*ofac)  
  fndim = ndim  
  ck    = ((x-xmin)*fac) % fndim  
  ckk   = (2.0*ck) % fndim  
  
  for j in range(0L, n):  
    __spread__(y[j]-ave,wk1,ndim,ck[j],MACC)  
    __spread__(1.0,wk2,ndim,ckk[j],MACC)  
  
  # Take the Fast Fourier Transforms  
  wk1 = fft.ifft( wk1 )*len(wk1)  
  wk2 = fft.ifft( wk2 )*len(wk1)  
  
  wk1  = wk1[1:nout+1]  
  wk2  = wk2[1:nout+1]  
  rwk1 = wk1.real  
  iwk1 = wk1.imag  
  rwk2 = wk2.real  
  iwk2 = wk2.imag  
    
  df  = 1.0/(xdif*ofac)  
    
  # Compute the Lomb value for each frequency  
  hypo2 = 2.0 * abs( wk2 )  
  hc2wt = rwk2/hypo2  
  hs2wt = iwk2/hypo2  
  
  cwt  = np.sqrt(0.5+hc2wt)  
  swt  = np.sign(hs2wt)*(np.sqrt(0.5-hc2wt))  
  den  = 0.5*n+hc2wt*rwk2+hs2wt*iwk2  
  cterm = (cwt*rwk1+swt*iwk1)**2./den  
  sterm = (cwt*iwk1-swt*rwk1)**2./(n-den)  
  
  wk1 = df*(np.arange(nout, dtype='float')+1.)  
  wk2 = (cterm+sterm)/(2.0*var)  
  
  return wk1,wk2

def fwls( x, y, e, ofac, hifac, MACC=4):  
  """ function fwls
    
  Weighted Lomb-Scargle periodogram, evaluated using essentially
  the same method as fls. See that for more details of the methods.

  The difference here is that it uses inverse variance weights. Note that
  a weighted average is subtracted and so you do not need to make any such
  adjustment before using the routine. This is not however a "floating mean"
  periodogram because the mean is subtracted from the start.
   
  Arguments: 
     x     : abscissas array, (e.g. an array of times). 
     y     : ordinates array, (e.g. corresponding counts). 
     e     : uncertainties in the y values
     ofac  : oversampling factor. 
     hifac : hifac * "average" Nyquist frequency = highest frequency 
             for which values of the Lomb normalized periodogram will 
             be calculated. 
     MACC  : number of interpolation points per 1/4 cycle 
             of highest frequency        

    Returns: 
      f : An array of Lomb periodogram frequencies. 
      p : An array of corresponding values of the Lomb periodogram. 
 
    Reference:  
      Press, W. H. & Rybicki, G. B. 1989, ApJ vol. 338, p. 277-280. 
      Fast algorithm for spectral analysis of unevenly sampled data 
      (1989ApJ...338..277P) 

    History: 
      23/12/2011, v1.0, TRM, adapted from fls. 
  """  

  # Check dimensions of input arrays  
  n = long(len(x))  
  if n != len(y):
    raise Exception('fwls: incompatible arrays lengths')

  nout   = int(np.ceil(0.5*ofac*hifac*n))  

  nfreqt = long(ofac*hifac*n*MACC)   # Size the FFT as next power  
  nfreq  = 64L                       # of 2 above nfreqt.  
  while nfreq < nfreqt:   
    nfreq = 2*nfreq  
  ndim = long(2*nfreq)  
    
  # Compute weighted mean
  w     = 1/e**2
  wsum  = w.sum()
  wave  = (w*y).sum() / wsum
 
  # and range of the data.  
  xmin = x.min()  
  xmax = x.max()  
  xdif = xmax-xmin  
  
  # extirpolate the data into the workspaces  
  wk1 = np.zeros(ndim, dtype='complex')  
  wk2 = np.zeros(ndim, dtype='complex')  
  
  fac   = ndim/(xdif*ofac)  
  fndim = ndim  
  ck    = ((x-xmin)*fac) % fndim  
  ckk   = (2.0*ck) % fndim 

  for j in range(0L, n):
    __spread__(w[j]*(y[j]-wave),wk1,ndim,ck[j],MACC)  
    __spread__(w[j],wk2,ndim,ckk[j],MACC)  
  
  # Take the Fast Fourier Transforms  
  wk1 = fft.ifft( wk1 )*len(wk1)  
  wk2 = fft.ifft( wk2 )*len(wk1)  
  
  wk1  = wk1[1:nout+1]  
  wk2  = wk2[1:nout+1]  
  rwk1 = wk1.real  
  iwk1 = wk1.imag  
  rwk2 = wk2.real  
  iwk2 = wk2.imag  
    
  df  = 1.0/(xdif*ofac)  
    
  # Compute the Lomb value for each frequency  
  hypo2 = 2.0 * abs( wk2 )  
  hc2wt = rwk2/hypo2  
  hs2wt = iwk2/hypo2  
  
  cwt  = np.sqrt(0.5+hc2wt)  
  swt  = np.sign(hs2wt)*(np.sqrt(0.5-hc2wt))  
  den  = 0.5*wsum+hc2wt*rwk2+hs2wt*iwk2  
  cterm = (cwt*rwk1+swt*iwk1)**2./den  
  sterm = (cwt*iwk1-swt*rwk1)**2./(wsum-den)  
  
  wk1 = df*(np.arange(nout, dtype='float')+1.)  
  wk2 = (cterm+sterm)/2.0  
  
  return wk1,wk2

def fwmls( x, y, e, ofac, hifac, amp=False, MACC=4):  
  """ function fwmls
    
  Weighted, modified Lomb-Scargle periodogram, evaluated using essentially
  the same method as fls. See that for more details of the methods.

  The difference here is that it uses inverse variance weights and allows
  for a "floating mean". The final result = 1/2 the difference in chi**2
  between fitting a constant plus a sinusoid versus just a constant.

  For L-S type work, this routine should be the default one to choose.

  Arguments: 
     x     : abscissas array, (e.g. an array of times). 
     y     : ordinates array, (e.g. corresponding counts). 
     e     : uncertainties in the y values
     ofac  : oversampling factor. 
     hifac : hifac * "average" Nyquist frequency = highest frequency 
             for which values of the Lomb normalized periodogram will 
             be calculated. 
     amp   : If True, will return the amplitude rather than the usual
             half change in chi**2.
     MACC  : number of interpolation points per 1/4 cycle 
             of highest frequency        

    Returns: 
      f : An array of Lomb periodogram frequencies. 
      p : An array of corresponding values of the Lomb periodogram (or amplitudes) 
 
    Reference:  
      Press, W. H. & Rybicki, G. B. 1989, ApJ vol. 338, p. 277-280. 
      Fast algorithm for spectral analysis of unevenly sampled data 
      (1989ApJ...338..277P) 

    History: 
      23/12/2011, v1.0, TRM, adapted from fls. 
  """  

  # Check dimensions of input arrays  
  n = long(len(x))  
  if n != len(y):
    raise Exception('fwls: incompatible arrays lengths')

  nout   = int(np.ceil(0.5*ofac*hifac*n))  

  nfreqt = long(ofac*hifac*n*MACC)   # Size the FFT as next power  
  nfreq  = 64L                       # of 2 above nfreqt.  
  while nfreq < nfreqt:   
    nfreq = 2*nfreq  
  ndim = long(2*nfreq)  
    
  # Compute weighted mean
  w     = 1/e**2
  wsum  = w.sum()
  wave  = (w*y).sum() / wsum
 
  # and range of the data.  
  xmin = x.min()  
  xmax = x.max()  
  xdif = xmax-xmin  
  
  # "extirpolate" the data into the workspaces  
  # a trick to carry out double frequency FFTs
  # is used.
  wk1 = np.zeros(ndim, dtype='complex')  
  wk2 = np.zeros(ndim, dtype='complex')  
  wk3 = np.zeros(ndim, dtype='complex')  
  
  fac   = ndim/(xdif*ofac)  
  fndim = ndim  
  ck    = ((x-xmin)*fac) % fndim  
  ckk   = (2.0*ck) % fndim 

  for j in range(0L, n):
    __spread__(w[j]*(y[j]-wave),wk1,ndim,ck[j],MACC)  
    __spread__(w[j],wk2,ndim,ck[j],MACC)  
    __spread__(w[j],wk3,ndim,ckk[j],MACC)  

  # Take the Fast Fourier Transforms  
  wk1 = fft.ifft( wk1 )*len(wk1)  
  wk2 = fft.ifft( wk2 )*len(wk1)  
  wk3 = fft.ifft( wk3 )*len(wk1)  
  
  # Some reduction of memory usage occurs here
  wk1  = wk1[1:nout+1]  
  wk2  = wk2[1:nout+1]  
  wk3  = wk3[1:nout+1]
  
  # now have sums of weight*cos()*y etc.
  # wc2 = sum weight*cos(2*k*x), etc, allowing
  # computation of sum weight*cos()**2 etc.
  wcy = wk1.real  
  wsy = wk1.imag  
  wc  = wk2.real  
  ws  = wk2.imag  
  wc2 = wk3.real  
  ws2 = wk3.imag  
    
  df  = 1.0/(xdif*ofac)  

  # "Normal equation" approach followed here where we think
  # of linear least-squares fit of y = c + a*cos() + b*sin()
  # End up with 3x3 matrix A of form:
  #
  # A = w   wc  ws
  #     wc  wcc wcs
  #     ws  wcs wss
  # 
  # and a vector b of the form
  #
  # b = wy
  #     wcy
  #     wsy
  #
  # (all cases summed over the data points). In fact we have one of these
  # for each frequency. The FFTs above have given us the values of the sums
  # evaluated at all frequencies. The solution vector is Inv(A)*b while the
  # difference in chi**2 we are after is the quadratic form Trans(b)*Inv(A)*b.
  # The need to invert A gives rise to the references to cofactors and 
  # determinants in what follows. The one other feature to realise is that the
  # subtraction of wave = wy / w above, ensures that the wy (first element in the
  # b vector) at this point = 0, saving some computation.

  # Compute (twice) the determinant for each frequency,
  # cf13 = co-factor row 1, col 3, etc.
  cf13 = (wc*ws2-ws*(wsum+wc2))/2.
  cf23 = wc*ws-(wsum/2.)*ws2
  cf33 = (wsum/2.)*(wsum+wc2)-wc**2

  det2 = 2.*ws*cf13 + ws2*cf23 + (wsum-wc2)*cf33

  if amp:
    cf22 = (wsum/2.)*(wsum-wc2)-ws**2
    wk2  = 2*np.sqrt((cf22*wcy + cf23*wsy)**2 + (cf23*wcy + cf33*wsy)**2) / det2
  else:
    wk2 = (((wsum/2.)*(wsum-wc2)-ws**2)*wcy**2 + 2*cf23*wcy*wsy + cf33*wsy**2) / det2

  # Generate frequency grid and return.
  wk1 = df*(np.arange(nout, dtype='float')+1.)  

  return wk1,wk2

def genf( x, ofac, hifac):
  """ function genf
    
  Generates frequency grid in the same way as fls, fwls, and fwmls etc, but without any 
  computation of periodograms. Useful for feeding frequencies into wmls

  Arguments: 
     x     : abscissas array, (e.g. an array of times). 
     ofac  : oversampling factor. 
     hifac : hifac * "average" Nyquist frequency = highest frequency 
             for which values of the Lomb normalized periodogram will 
             be calculated. 

    Returns: 
      f : An array of frequencies. 
  """  

  nout = int(np.ceil(0.5*ofac*hifac*len(x)))  
  xdif = x.max()-x.min()
  df   = 1.0/(xdif*ofac)  
  return df*(np.arange(nout, dtype='float')+1.)  

def wmls( x, y, e, f, amp=False, nh=1):  
  """ function wmls
    
  Weighted, modified, Lomb-Scargle periodogram, evaluated directly, with
  relatively little memory usage, but potentially slowly for full frequency
  evaluation, although this does depend upon the exact distribution of x values.
  Extremely clumpy data might actually do better with this routine than the
  P&R based ones which "extirpolate" the data onto a fine grid before FFT-ing.
  It should return essentially the same values as fwmls for the same problem.

  Like fwmls this uses inverse variance weights and allows for a "floating 
  mean". However, in addition it can allow for harmonics. The final result 
  = 1/2 the difference in chi**2 between fitting a constant plus a sinusoid 
  (or a harmonic series of sinusoids) versus just a constant.
   
  The harmonic series should in principle allow more sensitive searches for 
  non-sinusoidal signals.

  Arguments: 
     x     : abscissas array, (e.g. an array of times). 
     y     : ordinates array, (e.g. corresponding counts). 
     e     : uncertainties in the y values
     f     : frequency array to compute periodogram upon.
     amp   : True to return the amplitude rather than the delta chi**2 (only if nh=1)
     nh    : number of harmonics (1 = fundamental only)

    Returns: 
      p : An array of corresponding values of the Lomb periodogram. 
 
  """  

  if amp and nh > 1:
    raise Exception('wmls: amp=True not supported for nh > 1')

  # Check dimensions of input arrays  
  n = long(len(x))  
  if n != len(y):
    raise Exception('wmls: incompatible arrays lengths')
 
  wk = np.empty((1+2*nh,n),dtype='float')
  wk[0,:] = 1.

  # Prepare x and y
  xoff = (x.min()+x.max())/2.
  xd   = x - xoff

  w     = 1/e**2
  wsum  = w.sum()
  wave  = (w*y).sum() / wsum
  yd    = y - wave

  p   = np.empty_like(f)
  A   = np.matrix(np.empty((2*nh+1,2*nh+1)))
  b   = np.matrix(np.empty((2*nh+1))).transpose()
  for n,freq in enumerate(f):

    for m in range(nh):
      wk[2*m+1,:] = np.cos((2.*np.pi*(m+1)*freq)*xd)
      wk[2*m+2,:] = np.sin((2.*np.pi*(m+1)*freq)*xd)

    for m in range(2*nh+1):
      b[m] = (w*yd*wk[m,:]).sum()
      for k in range(m+1):
        A[m,k] = (w*wk[m,:]*wk[k,:]).sum()
        A[k,m] = A[m,k]
    
    AI = A.getI()
    bt = b.transpose() 
    if amp:
      p[n] = np.sqrt((AI[1,1]*b[1]+AI[1,2]*b[2])**2 + (AI[2,1]*b[1]+AI[2,2]*b[2])**2)
    else:
      p[n] = float(bt*AI*b)/2.
  return p

def pdm( x, y, e, f, nbin, xshift=0.):  
  """ function pdm
    
  Phase dispersion minimisation periodogram. Bins data into nbin bins,
  returns dispersion calculated with separate offsets for each bin versus
  overall dispersion. Advantages: good for non-sinusoidal data. 
  Disadvantages: does not vary smoothly with f, O(N**2) calculation.
   
  Arguments: 
     x      : abscissas array, (e.g. an array of times). 
     y      : ordinates array, (e.g. corresponding counts). 
     e      : uncertainties in the y values
     f      : frequency array to compute periodogram upon
     nbin   : number of bins
     xshift : amount to shift in phase

    Returns: 
      p : An array of phase dispersion values
 
  """  

  # Check dimensions of input arrays  
  n = long(len(x))  
  if n != len(y):
    raise Exception('wmls: incompatible arrays lengths')
 
  # Prepare x and y
  xoff = (x.min()+x.max())/2.
  xd   = x - xoff

  w     = 1/e**2
  wsum  = w.sum()
  wave  = (w*y).sum() / wsum
  yd    = y - wave
  vtot  = (w*(yd-wave)**2).sum()

  p = np.empty_like(f)

  for n, freq in enumerate(f):
    bin = (nbin*freq*xd % nbin).astype('int')
    vbin = 0.
    for m in xrange(nbin):
      ok = bin == m
      if len(bin[ok]) > 1:
        wsbin = w[ok].sum()
        wabin = (w[ok]*yd[ok]).sum() / wsbin
        vbin += (w[ok]*(yd[ok]-wabin)**2).sum()
    p[n] = vbin / vtot

  return p
  
