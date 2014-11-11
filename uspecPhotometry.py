#!/usr/bin/env python
import sys
sys.path.remove('/home/astro/phsgan/Python64/lib/python/site-packages/astropy-0.3.2-py2.6-linux-x86_64.egg')

import ultracamutils
import matplotlib.pyplot
import argparse
import numpy, math
import classes
import ultraspecClasses
import ultraspecutils
#import rashley_utils as utils
from trm import ultracam
from trm.ultracam.UErrors import PowerOnOffError, UendError, UltracamError
import ultracam_shift
import time, datetime
import json
import Image,ImageDraw
import ucamObjectClass
from photutils import datasets
from photutils import daofind
from photutils import aperture_photometry, CircularAperture, psf_photometry, GaussianPSF
import astropy.table
import astropy.time
from astropy.stats import median_absolute_deviation as mad

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Reads the Ultraspec [dd-mm-yyyy/runxxx.dat] files produces previews of the images')
	parser.add_argument('runname', type=str, help='Ultracam run name  [eg 2013-07-21/run010]')
	parser.add_argument('-p', '--preview', action='store_true', help='Show image previews with Matplotlib')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf')
	parser.add_argument('-d', '--debuglevel', type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	parser.add_argument('--startframe', default=1, type=int, help='Start frame. \'1\' is the default')
	parser.add_argument('-n', '--numframes', type=int, help='Number of frames. No parameter means all frames, or from startframe to the end of the run')
	parser.add_argument('-t', '--sleep', default=0, type=int, help='Sleep time (in seconds) between frames. \'0\' is the default')
	parser.add_argument('-o', default = 0, type=int, help = "Only plot photometry for the top [o] objects.")
	
	arg = parser.parse_args()

	config = ultracamutils.readConfigFile(arg.configfile)
	
	debug = classes.debugObject(config.DEBUG)
	debug.toggleTimeLog()
	if (arg.debuglevel!=None): debug.setLevel(arg.debuglevel);
	
	
	runFilename = ultracamutils.addPaths(config.ULTRASPECRAW, arg.runname)

	debug.write("Opening the Ultraspec raw file at: " + runFilename, level = 3)
	
	runDate, runID = ultracamutils.separateRunNameAndDate(arg.runname)
	
	""" Check that the working folders and the output folders are there
	"""
	(runDate, runNumber) = ultracamutils.separateRunNameAndDate(arg.runname)
	
	workingFolder = ultracamutils.addPaths(config.WORKINGDIR, runDate)
	outputFolder = ultracamutils.addPaths(config.SITE_PATH, runDate)
	ultracamutils.createFolder(workingFolder)
	ultracamutils.createFolder(outputFolder)
	
	""" Load the apertures from an existing aperture file
	"""
	apertureFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_apertures.json"
	debug.write("Loading apertures from file: %s"%apertureFilename, level = 2)
	apertures = ultraspecutils.loadApertures(apertureFilename)
	debug.write("Loaded %d apertures"%(len(apertures)), level = 2)
	# Set a radius and sky aperture for all of the apertures
	for a in apertures:
		a.setRadius(5)
		a.setSkyAperture(8, 12)
		
	x = []
	y = []
	for a in apertures:
		x.append(a.position[0])
		y.append(a.position[1])

	positions = zip(x, y)
	
	startFrame = arg.startframe
	if startFrame<1:
		debug.error("startframe cannot be less than 1")
		sys.exit()

	rdat  = ultracam.Rdata(runFilename, startFrame, server=False)

	maximumFrames = rdat.ntotal()
	debug.write("Total number of frames in the run is %d"%maximumFrames, level = 2 )
	if startFrame>maximumFrames:
		debug.error("startframe " + str(startFrame) + ", is beyond the end of the run, which has only " + str(maximumFrames) + " frames in it.")
		sys.exit()
		
	frameRange = maximumFrames - startFrame + 1
	
	if arg.numframes!=None:
		requestedNumFrames = arg.numframes
		if requestedNumFrames<(frameRange):
			frameRange = requestedNumFrames
	
	startTime = datetime.datetime.now()
	timeLeftString = "??:??"
	""" Run through all the frames in the .dat file.
	"""
	if arg.preview:
		matplotlib.pyplot.figure(figsize=(8, 8))
		matplotlib.pyplot.ion()
		imageFigure = matplotlib.pyplot.gcf()
		matplotlib.pyplot.figure(1)
		matplotlib.pyplot.title("Frame image")
			
	#fullFrame = numpy.zeros((1057, 1040))
			
	allWindows = []

	ccdFrame = rdat()
	#ccdFrame.rback()
	window = ccdFrame[0]
	
	for windowIndex, w in enumerate(window):
		# Set up some info about the window sizes and extents
		window = ultraspecClasses.window()
		window.setExtents(w.llx, w.lly, w.nx, w.ny)
		window.setBinning(w.xbin, w.ybin)
		image = w._data
		window.setData(image)
		allWindows.append(window)
		
	(xmin, ymin, xmax, ymax) = ultraspecutils.determineFullFrameSize(allWindows)
	fullFramexsize = xmax - xmin
	fullFrameysize = ymax - ymin
	
	allPhotometry = []
			
	for frameIndex in range(2, frameRange + 1):
		framesToGo = frameRange - frameIndex
		currentTime = datetime.datetime.now()
		trueFrameNumber = startFrame + frameIndex - 1
		completionPercent = (float(frameIndex) / float(frameRange) * 100.)
		timePassed = ultracamutils.timedeltaTotalSeconds(currentTime - startTime)
		totalTime = timePassed * 100. / completionPercent
		etaTime = startTime + datetime.timedelta(seconds = totalTime)
		timeLeft = etaTime - currentTime
		(hours, mins, secs) = ultracamutils.timedeltaHoursMinsSeconds(timeLeft)
		timeLeftString = str(hours).zfill(2) + ":" + str(mins).zfill(2) + ":" + str(secs).zfill(2)
		
		""" Read the frame and the timestamp
		"""
		ccdFrame = rdat()
		windows = ccdFrame[0]
		MJD = windows.time.mjd
		
		""" Display the status of the reduction so far
		"""
		statusString = "\r%s MJD: %5.7f Frame: [%d/%d]"%(timeLeftString, MJD, trueFrameNumber, frameRange)
		sys.stdout.write(statusString)
		sys.stdout.flush()
			
		for windowIndex, w in enumerate(windows):
			image = w._data
			allWindows[windowIndex].setData(image)
			
		if arg.preview: 
			fullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
			for w in allWindows:
				boostedImage = ultracamutils.percentiles(w.data, 20, 99)
				xll = w.xll/w.xbin - xmin
				xsize = w.nx
				yll = w.yll/w.ybin - ymin
				ysize = w.ny
				fullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + boostedImage
					
			
			matplotlib.pyplot.imshow(fullFrame, cmap='gray_r')
			
			for a in apertures:
				(x, y) = a.position
				matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), a.radius, color='green', fill=False, linewidth=1.0))
				matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), a.skyInnerRadius, color='blue', fill=False, linewidth=1.0))
				matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), a.skyOuterRadius, color='blue', fill=False, linewidth=1.0))
			
			
			matplotlib.pyplot.figure(1)
			matplotlib.pyplot.title("Frame image [%d/%d]"%(trueFrameNumber, frameRange))
			matplotlib.pyplot.gca().invert_yaxis()
			matplotlib.pyplot.draw()
			matplotlib.pyplot.show()
			matplotlib.pyplot.clf()    # This clears the figure in matplotlib and fixes the 'memory leak'
		
		""" Now perform the photometry.
		"""
		fullDataFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
		for w in allWindows:
			#data = ultracamutils.percentiles(w.data, 20, 99)
			data = w.data
			xll = w.xll/w.xbin - xmin
			xsize = w.nx
			yll = w.yll/w.ybin - ymin
			ysize = w.ny
			fullDataFrame[yll:yll+ysize, xll:xll+xsize] = fullDataFrame[yll:yll+ysize, xll:xll+xsize] + data
		
		# Try annulus subtraction of the background
		from photutils import CircularAnnulus
		radius = apertures[0].radius
		skyInner = apertures[0].skyInnerRadius
		skyOuter = apertures[0].skyOuterRadius
		photutilsApertures = CircularAperture(positions, r=radius)
		annulusApertures = CircularAnnulus(positions, r_in= skyInner, r_out=skyOuter)
		rawflux_table = aperture_photometry(fullDataFrame, photutilsApertures)
		bkgflux_table = aperture_photometry(fullDataFrame, annulusApertures)
		phot_table = astropy.table.hstack([rawflux_table, bkgflux_table], table_names=['raw', 'bkg'])
		aperture_area = numpy.pi * radius ** 2
		annulus_area = numpy.pi * (skyOuter ** 2 - skyInner ** 2)
		bkg_sum = phot_table['aperture_sum_bkg'] * aperture_area / annulus_area
		final_sum = phot_table['aperture_sum_raw'] - bkg_sum
		phot_table['residual_aperture_sum'] = final_sum
		#print phot_table['residual_aperture_sum']   
		
		framePhotometry = {}
		framePhotometry['MJD'] = MJD
		results = phot_table['residual_aperture_sum']
		for index, r in enumerate(results):
			framePhotometry[index] = r
			
		allPhotometry.append(framePhotometry)
			
		if arg.sleep!=0:
			time.sleep(arg.sleep)
			
		""" End of the main loop 
		"""
	sys.stdout.write("\rProcessed %d frames                                \n"%frameRange)
	sys.stdout.flush()
	
	matplotlib.pyplot.figure(figsize=(12, 12))
	xValues = []
	for p in allPhotometry:
		xValues.append(p['MJD'])
		
	MJDoffset = int(min(xValues))
	print "MJD offset:",MJDoffset
	
	
	xValues = [x - MJDoffset for x in xValues]
	
	
	totalObjects = len(apertures)
	for index in range(totalObjects):
		yValues = []
		for p in allPhotometry:
			yValues.append(p[index])
		matplotlib.pyplot.plot(xValues, yValues, '.', label = 'i')
		
	matplotlib.pyplot.xlabel('MJD' + ' +' + str(MJDoffset))
	matplotlib.pyplot.show(block=True)
	
	
