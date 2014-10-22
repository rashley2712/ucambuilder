"""
A module of binary star related routines.
a      : calculates orbital separation
period : calculates period
chirp  : calculates chirp mass 
jdotgr : computes GR-induced angular momentum loss rate.
j      : calculates angular momentum
"""

import math as m
from trm import subs

def chirp(m1, m2):
    """
    Returns the chirp mass of a binary
    """
    M  = m1+m2
    mu = m1*m2/M
    return mu**0.6*M**0.4

def jdotgr(m1, m2, a):
    """
    Computes the rate of change of orbital angular momentum due to GRW losses
    divided by the orbital angular momentum (Jdot/J) for two point masses in a 
    circular orbit. Returned units are sec^-1
 
    m1 -- mass of star 1 (solar masses)
    m2 -- mass of star 2 (solar masses)
    a  -- orbital separation (solar radii)
    """
    
    return -32./5.*(subs.G*subs.MSUN/subs.C/subs.RSUN)**3/subs.C**2/subs.RSUN*m1*m2*(m1+m2)/a**4

def j(m1, m2, a):
    """
    Computes orbital angular momentum of two point masses in a 
    circular orbit. SI units returned.
 
    m1 -- mass of star 1 (solar masses)
    m2 -- mass of star 2 (solar masses)
    a  -- orbital separation (solar radii)
    """
    return subs.MSUN*m1*m2*m.sqrt(subs.G*subs.MSUN*subs.RSUN*a/(m1+m2))

def a(m1, m2, period):
    """
    Computes orbital separation in solar radii of two point masses in a circular
    ornit.
 
    m1      -- mass of star 1 (solar masses)
    m2      -- mass of star 2 (solar masses)
    period  -- orbital period (seconds)
    """
    return (subs.G*subs.MSUN*(m1+m2)*(period/2./m.pi)**2)**(1./3.)/subs.RSUN

def period(m1, m2, a):
    """
    Computes orbital period in seconds
    ornit.
 
    m1      -- mass of star 1 (solar masses)
    m2      -- mass of star 2 (solar masses)
    a       -- solar radii (solar radii)
    """
    return 2.*m.pi*((subs.RSUN*a)**3/(subs.G*subs.MSUN*(m1+m2)))**0.5
