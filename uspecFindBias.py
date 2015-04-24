#!/usr/bin/env python
import sys
try: 
	sys.path.remove('/home/astro/phsgan/Python64/lib/python/site-packages/astropy-0.3.2-py2.6-linux-x86_64.egg')
except (ValueError):
	print "No need to fix sys.path"
import os
import ultracamutils
import matplotlib.pyplot
import argparse
import numpy, math
import classes
import ultraspecClasses
from trm import ultracam
from trm.ultracam.UErrors import PowerOnOffError, UendError, UltracamError
import ultracam_shift
import time, datetime
import json
from scipy import ndimage
import Image
import ucamObjectClass
from photutils import datasets
from photutils import daofind
from photutils import aperture_photometry, CircularAperture, psf_photometry, GaussianPSF
import astropy.table, astropy.io
from astropy.stats import median_absolute_deviation as mad
import astropy.stats.sigma_clipping
from astropy.stats import sigma_clipped_stats
from astropy.convolution import Gaussian2DKernel
from photutils.detection import detect_sources
from scipy.ndimage import binary_dilation
from photutils.background import Background
import re

def checkWindowBinMatch(windows1, windows2):
	if len(windows1) != len(windows2): return False
	for index, w1 in enumerate(windows1):
		w2 = windows2[index]
		if w1['xll']  != w2['xll']:  return False
		if w1['yll']  != w2['yll']:  return False
		if w1['nx']   != w2['nx']:   return False
		if w1['ny']   != w2['ny']:   return False
		if w1['xbin'] != w2['xbin']: return False
		if w1['ybin'] != w2['ybin']: return False
	
	return True
	
def compareWindowLists(windows1, windows2):
	if len(windows1) != len(windows2): return False
	
	return True

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Looks in the ULTRASPEC archive for a run that might be used as a bias for this run.')
	parser.add_argument('runname', type=str, help='ULTRASPEC  [eg 2013-07-21/run010]')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf')
	parser.add_argument('-d', '--debuglevel', type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	
	arg = parser.parse_args()
	print arg
	config = ultracamutils.readConfigFile(arg.configfile)

	if arg.debuglevel==None:
		debug = classes.debugObject(config.DEBUG)
	else:
		debug = classes.debugObject(arg.debuglevel)
	
	debug.toggleTimeLog()
	
	runInfo = ultraspecClasses.runInfo(arg.runname)
	found = runInfo.loadFromJSON(config.RUNINFO)
	if not found:
		debug.write("Could not get info for this run from the ultra.json file.", level = 1)
		xmlRead = runInfo.loadFromXML(config.ULTRASPECRAW)
	debug.write(runInfo, 2)
		
	runFilename = ultracamutils.addPaths(config.ULTRASPECRAW, arg.runname)

	debug.write("Opening the ULTRASPEC raw file at: " + runFilename, level = 2)
	
	runDate, runID = ultracamutils.separateRunNameAndDate(arg.runname)
	
	""" Check that the working folders and the output folders are there
	"""
	(runDate, runNumber) = ultracamutils.separateRunNameAndDate(arg.runname)
	
	runInfo = ultraspecClasses.runInfo(arg.runname)
	runInfo.loadFromXML(config.ULTRASPECRAW)
	thisRunWindows = runInfo.windowInfo
	print "RunInfo", runInfo
	
	workingFolder = ultracamutils.addPaths(config.WORKINGDIR, runDate)
	ultracamutils.createFolder(workingFolder)
	
	maxDays = 7
	runDay = datetime.datetime.strptime(runDate,  "%Y-%M-%d")
	oneDay = datetime.timedelta(days=1)
	biasRun = ""
	for dayNumber in range(maxDays):
		dateString = runDay.strftime("%Y-%M-%d")
		debug.write("Looking for biases on this date: %s"%dateString, 2)
		runList = ultracamutils.getRunsByDay(config.ULTRASPECRAW, dateString)
		for run in runList:
			runPath = config.ULTRASPECRAW + "/" + run
			runInfo = ultraspecClasses.runInfo(run)
			runInfo.loadFromXML(config.ULTRASPECRAW)
			if runInfo.flags=='bias':
				if (checkWindowBinMatch(thisRunWindows, runInfo.windowInfo)):
					biasRun = run
					break
		runDay = runDay + oneDay
		if biasRun!="": break
	if biasRun=="":
		print "No biases found in the next %d days."%(maxDays), " Searching backwards for %d days."%(maxDays)
		runDay = datetime.datetime.strptime(runDate,  "%Y-%M-%d")
		runDay-= oneDay
		for dayNumber in range(maxDays-1):
			dateString = runDay.strftime("%Y-%M-%d")
			debug.write("Looking for biases on this date: %s"%dateString, 2)
			runList = ultracamutils.getRunsByDay(config.ULTRASPECRAW, dateString)
			for run in runList:
				runPath = config.ULTRASPECRAW + "/" + run
				runInfo = ultraspecClasses.runInfo(run)
				runInfo.loadFromXML(config.ULTRASPECRAW)
				if runInfo.flags=='bias':
					if (checkWindowBinMatch(thisRunWindows, runInfo.windowInfo)):
						biasRun = run
						break
			runDay = runDay - oneDay
			if biasRun!="": break
			
	if biasRun=="":
		print "Sorry no biases found after looking forward and backward %d days."%maxDays
		sys.exit()
	
	print "The bias for %s is %s."%(arg.runname, biasRun)
	
	
	
	sys.exit()
		
	
	"""debug.write("Examining the first few frames of the run...", level = 2)
	originalStackedFrame = numpy.zeros((1057, 1040))
	
	ccdFrame = rdat()
	window = ccdFrame[0]
	allWindows = []	
	for windowIndex, w in enumerate(window):
		# Set up some info about the window sizes and extents
		window = ultraspecClasses.window()
		window.setExtents(w.llx, w.lly, w.nx, w.ny)
		window.setBinning(w.xbin, w.ybin)
		
		image = w._data
		allWindows.append(window)
		
	(xmin, ymin, xmax, ymax) = ultracamutils.determineFullFrameSize(allWindows)
	fullFramexsize = xmax - xmin
	fullFrameysize = ymax - ymin

	print "(xmin, xmax, ymin, ymax) (%d, %d, %d, %d)"%(xmin, xmax, ymin, ymax) 
	
	for w in allWindows:
		print "Window:", w
	"""
