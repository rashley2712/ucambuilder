#!/usr/bin/env python

import ultracamutils
import matplotlib.pyplot
import argparse
import numpy, math
import classes
import ultraspecClasses
#import rashley_utils as utils
from trm import ultracam
from trm.ultracam.UErrors import PowerOnOffError, UendError, UltracamError
import sys
import ultracam_shift
import time, datetime
import json
import Image,ImageDraw
import ucamObjectClass
from photutils import datasets
from photutils import daofind
from astropy.stats import median_absolute_deviation as mad
from   scipy.ndimage.filters import gaussian_filter


def determineFullFrameSize(windows):
	leftestPixel = 1057
	rightestPixel = 0
	topestPixel = 0 
	bottomestPixel = 1040
	for w in windows:
		if w.xll/w.xbin < leftestPixel: leftestPixel = w.xll/w.xbin
		if w.yll/w.ybin < bottomestPixel: bottomestPixel = w.yll/w.ybin
		if (w.xll/w.xbin + w.nx) > rightestPixel: rightestPixel = w.xll/w.xbin + w.nx
		if (w.yll/w.ybin + w.ny) > topestPixel: topestPixel = w.yll/w.ybin + w.ny
			
	return leftestPixel, bottomestPixel, rightestPixel, topestPixel
		

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Reads the Ultraspec [dd-mm-yyyy/runxxx.dat] files produces previews of the images')
	parser.add_argument('runname', type=str, help='Ultracam run name  [eg 2013-07-21/run010]')
	parser.add_argument('-p', '--preview', action='store_true', help='Show image previews with Matplotlib')
	parser.add_argument('-s', '--stack', action='store_true', help='Stack the images in the preview window')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf')
	parser.add_argument('-d', '--debuglevel', type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	parser.add_argument('--startframe', default=1, type=int, help='Start frame. \'1\' is the default')
	parser.add_argument('-n', '--numframes', type=int, help='Number of frames. No parameter means all frames, or from startframe to the end of the run')
	parser.add_argument('-t', '--sleep', default=0, type=int, help='Sleep time (in seconds) between frames. \'0\' is the default')
	
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
		fig = matplotlib.pyplot.gcf()
		matplotlib.pyplot.title("Frame image")
		if arg.stack:
			matplotlib.pyplot.title("Stacked image")
			
	fullFrame = numpy.zeros((1057, 1040))
			
	allWindows = []

	ccdFrame = rdat()
	ccdFrame.rback()
	window = ccdFrame[0]
	
	for windowIndex, w in enumerate(window):
		# Set up some info about the window sizes and extents
		window = ultraspecClasses.window()
		window.setExtents(w.llx, w.lly, w.nx, w.ny)
		window.setBinning(w.xbin, w.ybin)
		
		image = w._data
		image -= numpy.median(image)
		window.setData(image)
		
		bkg_sigma = 1.48 * mad(image)
		print "bkg_sigma", bkg_sigma   
		sources = daofind(image, fwhm=4.0, threshold=3*bkg_sigma)   
		window.setSources(sources)	
		
		allWindows.append(window)
		
	(xmin, ymin, xmax, ymax) = determineFullFrameSize(allWindows)
	fullFramexsize = xmax - xmin
	fullFrameysize = ymax - ymin
	sourceMap = ultraspecClasses.sourceMap((fullFrameysize, fullFramexsize))
			

	for frameIndex in range(2, frameRange + 1):
		framesToGo = frameRange - frameIndex
		currentTime = datetime.datetime.now()
		trueFrameNumber = startFrame + frameIndex - 1
		ccdFrame = rdat()
		
		print currentTime, "Frame: [%d/%d]"%(trueFrameNumber, frameRange)
		
		ccdFrame.rback()
		window = ccdFrame[0]
		
		for windowIndex, w in enumerate(window):
			image = w._data
			image -= numpy.median(image)
			allWindows[windowIndex].addData(image)
			bkg_sigma = 1.48 * mad(image)
			sources = daofind(image, fwhm=4.0, threshold=3*bkg_sigma)   
			allWindows[windowIndex].setSources(sources)	
		
			
		# Combine the sources in all of the windows
		allSources = []
		for index, w in enumerate(allWindows):
			xll = w.xll/w.xbin - xmin
			yll = w.yll/w.ybin - ymin
			sources = w.getSources()
			for s in sources:
				(x, y) = s['xcentroid'], s['ycentroid']
				xAbs = x + xll
				yAbs = y + yll
				#print xll, yll, "(%d, %d) -> (%d, %d)"%(x, y, xAbs, yAbs)
				allSources.append((xAbs, yAbs))
		
		sourceMap.updateMap(allSources)
			
			
		if arg.preview: 
			fullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
			for w in allWindows:
				if (arg.stack):
					boostedImage = ultracamutils.percentiles(w.stackedData, 20, 99)
				else:
					boostedImage = ultracamutils.percentiles(w.data, 20, 99)
				xll = w.xll/w.xbin - xmin
				xsize = w.nx
				yll = w.yll/w.ybin - ymin
				ysize = w.ny
				fullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + boostedImage
					
			
			matplotlib.pyplot.imshow(fullFrame, cmap='gray')
			
			for s in allSources:
				(x, y) = s
				matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), 15, color='green', fill=False, linewidth=2.0))
			
			
			matplotlib.pyplot.title("Frame image [%d/%d]"%(trueFrameNumber, frameRange))
			if arg.stack:
				matplotlib.pyplot.title("Stacked image [%d/%d]"%(trueFrameNumber, frameRange))
			matplotlib.pyplot.gca().invert_yaxis()
			matplotlib.pyplot.draw()
			matplotlib.pyplot.show()
			matplotlib.pyplot.clf()    # This clears the figure in matplotlib and fixes the 'memory leak'
			
		if arg.sleep!=0:
			time.sleep(arg.sleep)

	# Generate the source map
	sourceMapImage = matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Source map")
	fwhm = 3
	psize = 1
	boostedImage = gaussian_filter(sourceMap.getSourceMap(),fwhm/psize/2.3548,mode='constant')
	#boostedImage = ultracamutils.percentiles(boostedImage, 1, 99)
	matplotlib.pyplot.imshow(boostedImage, cmap='hot')
	matplotlib.pyplot.gca().invert_yaxis()			
	matplotlib.pyplot.show(block=False)
	outputFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_sourcemap.png"
	matplotlib.pyplot.savefig(outputFilename)
	

	# Generate the stacked image for writing to disc
	stackedFigure = matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Stacked image")
	fullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
	for w in allWindows:
		boostedImage = ultracamutils.percentiles(w.stackedData, 40, 99)
		xll = w.xll/w.xbin - xmin
		xsize = w.nx
		yll = w.yll/w.ybin - ymin
		ysize = w.ny
		fullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + boostedImage
	
	image = matplotlib.pyplot.imshow(fullFrame, cmap='gray')
	matplotlib.pyplot.gca().invert_yaxis()			
	matplotlib.pyplot.show(block=True)
	
	outputFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + ".png"
	matplotlib.pyplot.savefig(outputFilename)
