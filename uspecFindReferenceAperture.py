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
from photutils.morphology import (centroid_com, centroid_1dg, centroid_2dg)
from photutils import CircularAperture
from photutils import CircularAnnulus
import photutils
import ppgplot

def shift_func(output_coords, xoffset, yoffset):
	return (output_coords[0] - yoffset, output_coords[1] - xoffset)


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
	
	parser = argparse.ArgumentParser(description='Performs aperture photometry for Ultraspec runs. [dd-mm-yyyy/runxxx.dat]')
	parser.add_argument('runname', type=str, help='Ultraspec run name  [eg 2013-07-21/run010]')
	parser.add_argument('-p', '--preview', action='store_true', help='Show image previews with Matplotlib')
	parser.add_argument('-s', '--stack', action='store_true', help='Stack the images in the preview window')
	parser.add_argument('--noshift', action='store_true', help='Don''t apply a linear shift to the stacked images to correct for drift.')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf')
	parser.add_argument('-d', '--debuglevel', type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	parser.add_argument('--startframe', default=1, type=int, help='Start frame. \'1\' is the default')
	parser.add_argument('-n', '--numframes', type=int, help='Number of frames. No parameter means all frames, or from startframe to the end of the run')
	parser.add_argument('-t', '--sleep', default=0, type=int, help='Sleep time (in seconds) between frames. \'0\' is the default')
	parser.add_argument('--xyls', action='store_true', help='Write an XYLS (FITS) file output catalog that can be used as input to Astronomy.net')
	parser.add_argument('--usefirstframe', action='store_true', help='Use the first frame of the run. Usually the first frame will be discarded.')
	parser.add_argument('-i', '--keyimages', action='store_true', help='Show some key images during the processing of this run.')
	parser.add_argument('--createconfig', action='store_true', help='Use this option to create a default configuration file.')
	
	arg = parser.parse_args()
	
	if arg.createconfig:
		ultracamutils.createConfigFile()
		sys.exit()
	
	config = ultracamutils.readConfigFile(arg.configfile)
	
	innerSkyRadius = float(config.INNER_SKY)
	outerSkyRadius = float(config.OUTER_SKY)
	apertureRadius = float(config.APERTURE_RADIUS)
	applyShift = True
	if (arg.noshift):
		applyShift = False
	debug = classes.debugObject(config.DEBUG)
	debug.toggleTimeLog()
	if (arg.debuglevel!=None): debug.setLevel(arg.debuglevel);
	debug.write(arg, level = 2)
	debug.write("Astropy version %s"%(astropy.__version__), level = 3)
	
	sourcesFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_sources.csv"
	debug.write("Loading source list from: " + sourcesFilename, 2)	
	sourceList = ultraspecClasses.sourceList()
	success = sourceList.loadFromCSV(sourcesFilename)
	if (not success):
		debug.error("Unable to open the list of sources. Have you run uspecCreateSourceMap yet?")
		sys.exit()
	else:
		debug.write("Loaded %d sources from the CSV file."%sourceList.getNumSources(), 2)
	
	runInfo = ultraspecClasses.runInfo(arg.runname)
	found = runInfo.loadFromJSON(config.RUNINFO)
	if not found:
		debug.write("Could not get info for this run from the ultra.json file.", level = 1)
		xmlRead = runInfo.loadFromXML(config.ULTRASPECRAW)
		
	debug.write(runInfo, 2)
		
	runFilename = ultracamutils.addPaths(config.ULTRASPECRAW, arg.runname)

	debug.write("Opening the Ultraspec raw file at: " + runFilename, level = 2)
	
	runDate, runID = ultracamutils.separateRunNameAndDate(arg.runname)
	
	""" Check that the working folders and the output folders are there
	"""
	(runDate, runNumber) = ultracamutils.separateRunNameAndDate(arg.runname)
	
	workingFolder = ultracamutils.addPaths(config.WORKINGDIR, runDate)
	ultracamutils.createFolder(workingFolder)
	
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
			
	fullFrame = numpy.zeros((1057, 1040))
			
	allWindows = []

	ccdFrame = rdat()
	frameWindows = ccdFrame[0]
	
	for windowIndex, w in enumerate(frameWindows):
		# Set up some info about the window sizes and extents
		window = ultraspecClasses.window()
		window.setExtents(w.llx, w.lly, w.nx, w.ny)
		window.setBinning(w.xbin, w.ybin)
		
		image = w._data
		if (arg.usefirstframe):
			window.setData(image)
			bkg_sigma = 1.48 * mad(image)
			sources = daofind(image, fwhm=4.0, threshold=3*bkg_sigma)   
			window.setSourcesAvoidBorders(sources)	
		else: 
			window.setBlankData(image)
		
		
		allWindows.append(window)
		
	(xmin, ymin, xmax, ymax) = determineFullFrameSize(allWindows)
	fullFramexsize = xmax - xmin
	fullFrameysize = ymax - ymin
	
	lightcurveView = ppgplot.pgopen('/xs')
	ppgplot.pgenv(startFrame, startFrame + frameRange, 0, 100, 0, 0)
	ppgplot.pgask(False)
	
	if (arg.preview):		
		bitmapView = ppgplot.pgopen('/xs')
		ppgplot.pgenv(0.,fullFramexsize,0.,fullFrameysize, 1, 0)
		pgPlotTransform = [0, 1, 0, 0, 0, 1]
		ppgplot.pgsfs(2)
	
					
	xValues = []
	yValues = []	
	yAxisMax= 100	
	referenceApertures = ultraspecClasses.referenceApertures()
	referenceApertures.initFromSourceList(sourceList)
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
		
		ccdFrame = rdat()
		
		statusString = "\r%s Frame: [%d/%d]"%(timeLeftString, trueFrameNumber, frameRange)
		sys.stdout.write(statusString)
		sys.stdout.flush()
		
		windows = ccdFrame[0]
		for windowIndex, w in enumerate(windows):
			image = w._data		
			allWindows[windowIndex].setData(image)
		
		
		if arg.preview:
			ppgplot.pgslct(bitmapView) 
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
			
			dimensions = numpy.shape(fullFrame)
			rows = dimensions[0]
			cols = dimensions[1]

			# Draw the grayscale bitmap
			ppgplot.pggray(fullFrame, 0, cols-1 , 0, rows-1 , 0, 255, pgPlotTransform)
	
			# Draw the full reference aperture list
			ppgplot.pgsci(3)
			for s in sourceList.getSources():
				(x, y) = s.abs_position
				ppgplot.pgcirc(x, y, 10)
		
		margins = 10
		
		if arg.preview: 
			ppgplot.pgslct(bitmapView)
			ppgplot.pgsci(2)
			
		plotColour = [1, 2, 3, 4, 5, 6]
		for index, s in enumerate(referenceApertures.getSources()):
			window = allWindows[s.windowIndex]
			center = s.latestPosition
			xcenterInt = int(center[0])
			xcenterOffset = center[0] - margins
			#xcenterOffset = xcenterInt - margins
			
			ycenterInt = int(center[1])
			ycenterOffset = center[1] - margins
			#ycenterOffset = ycenterInt - margins
			
			if (ycenterOffset<0) or (xcenterOffset<0): continue
			zoomImageData = window.data[ycenterInt-margins:ycenterInt+margins, xcenterInt-margins:xcenterInt+margins]
			(xcen, ycen) = photutils.morphology.centroid_2dg(zoomImageData, error=None, mask=None)
			xCollapsed = numpy.sum(zoomImageData, 0)
			yCollapsed = numpy.sum(zoomImageData, 1)
			
			xPeak = numpy.argmax(xCollapsed)
			yPeak = numpy.argmax(yCollapsed)
			print "xPeak, yPeak:", xPeak, yPeak
			
			xMat = [ [ (xPeak-1)**2, (xPeak-1) , 1 ] , \
			         [ (xPeak)**2,   (xPeak)   , 1 ] , \
			         [ (xPeak+1)**2, (xPeak+1) , 1 ] ]
			yMat = [ xCollapsed[xPeak-1], xCollapsed[xPeak], xCollapsed[xPeak+1] ]
			(a, b, c) = numpy.linalg.solve(xMat, yMat)
			newxPeak = -1.0 * b / (2.0 * a)
			xStartPlot = newxPeak - 5
			xEndPlot = newxPeak + 5
			xRange = xEndPlot - xStartPlot
			xpPoints = []
			ypPoints = []
			numPoints = 300
			xScale = xRange/numPoints
			for point in range(numPoints):
				xp = point*(xScale) + xStartPlot
				yp = a*xp*xp + b*xp + c
				xpPoints.append(xp)
				ypPoints.append(yp)

			xMat = [ [ (yPeak-1)**2, (yPeak-1) , 1 ] , \
			         [ (yPeak)**2,   (yPeak)   , 1 ] , \
			         [ (yPeak+1)**2, (yPeak+1) , 1 ] ]
			yMat = [ yCollapsed[yPeak-1], yCollapsed[yPeak], yCollapsed[yPeak+1] ]
			(a, b, c) = numpy.linalg.solve(xMat, yMat)
			newyPeak = -1.0 * b / (2.0 * a)
			xStartPlot = newyPeak - 5
			xEndPlot = newyPeak + 5
			xRange = xEndPlot - xStartPlot
			xyPoints = []
			yyPoints = []
			numPoints = 300
			xScale = xRange/numPoints
			for point in range(numPoints):
				xp = point*(xScale) + xStartPlot
				yp = a*xp*xp + b*xp + c
				xyPoints.append(xp)
				yyPoints.append(yp)

			print "Centroid method: (%f, %f)   vs   Quadratic fit: (%f, %f)"%(xcen, ycen, newxPeak, newyPeak) 
			
			centroidView = ppgplot.pgopen('/xs')
			ppgplot.pgenv(0, len(xCollapsed), 0, max(xCollapsed), 0, 0)
			pgPlotTransform = [0, 1, 0, 0, 0, 1]
			ppgplot.pglab("Pixels", "Counts position", "Source number: %d"%s.id)
			ppgplot.pgsci(2)
			ppgplot.pgpt(range(len(xCollapsed)), xCollapsed, 3)
			ppgplot.pgline(xpPoints, ypPoints)
			ppgplot.pgsci(3)
			ppgplot.pgpt(range(len(xCollapsed)), yCollapsed, 2)
			ppgplot.pgline(xyPoints, yyPoints)
			ppgplot.pgclos()
			#time.sleep(2)
			
			xcen+= xcenterOffset
			ycen+= ycenterOffset
			s.setLatestPosition(trueFrameNumber, (xcen, ycen))
			apertures = CircularAperture((xcen, ycen), r=apertureRadius)
			annulus_apertures = CircularAnnulus((xcen, ycen), r_in=innerSkyRadius, r_out=outerSkyRadius)
			
			# Draw the re-positioned apertures
			xll = window.xll/window.xbin - xmin
			yll = window.yll/window.ybin - ymin
			if arg.preview:
				ppgplot.pgslct(bitmapView)
				plotx= xcen + xll
				ploty= ycen + yll
				#print xll, yll, center, xcen, ycen
				ppgplot.pgcirc(plotx, ploty, apertureRadius)
				ppgplot.pgcirc(plotx, ploty, innerSkyRadius)
				ppgplot.pgcirc(plotx, ploty, outerSkyRadius)
				ppgplot.pgptxt(plotx-10, ploty-10, 0, 0, str(index)) 
				
			
			data = window.data
			
			innerFluxes = aperture_photometry(data, apertures)
			outerFluxes = aperture_photometry(data, annulus_apertures)
			
			innerFlux = innerFluxes[0]['aperture_sum']
			outerFlux = outerFluxes[0]['aperture_sum']
			
			bkg_mean = outerFlux / annulus_apertures.area()
			bkg_sum = bkg_mean * apertures.area()
			final_sum = innerFlux - bkg_sum
			#print "Sky subtracted flux:", final_sum
			s.addFluxMeasurement(trueFrameNumber, final_sum)
			
			xValues.append(trueFrameNumber)
			yValues.append(final_sum)
			yMin = 0
			yMax = numpy.max(yValues)
			shortXArray = [trueFrameNumber]
			shortYArray = [final_sum]
			if arg.preview: ppgplot.pgslct(lightcurveView)
			if yMax > yAxisMax:
				yAxisMax = yMax * 1.1
				ppgplot.pgenv(startFrame, startFrame + frameRange, yMin, yAxisMax, 0, 0)
				ppgplot.pgpt(xValues, yValues, 1)
				
			ppgplot.pgsci(plotColour[index])	
			ppgplot.pgpt(shortXArray, shortYArray, 1)
		
	
			"""ppgplot.pgsci(2)
			for index, s in enumerate(sourceList.getSources()):
				if index>5: break
				center = s.latestPosition
				window = allWindows[s.windowIndex]
				xll = window.xll/window.xbin - xmin
				yll = window.yll/window.ybin - ymin
				xcen= center[0] + xll
				ycen= center[1] + yll
				#print xll, yll, center, xcen, ycen
				ppgplot.pgcirc(xcen, ycen, 5)"""
			
			
		if arg.sleep!=0:
			time.sleep(arg.sleep)
	sys.stdout.write("\rProcessed %d frames      \n"%frameRange)
	sys.stdout.flush()
	
	if arg.preview:
		ppgplot.pgslct(bitmapView)
		ppgplot.pgclos()
	ppgplot.pgslct(lightcurveView)
	ppgplot.pgclos()
	
	
	# Find the reference aperture with most coverage
	referenceApertures.calculateFrameCoverage(frameRange-1)   # The '-1' is because we don't get photometry from the first frame
	referenceApertures.sortByCoverage()
	maxCoverage = referenceApertures.getSources()[0].coverage
	print "Max coverage:", maxCoverage
	referenceApertures.limitCoverage(maxCoverage)
	referenceApertures.sortByFlux()
	referenceAperture = referenceApertures.getSources()[0]
	print "Our reference aperture is:", referenceAperture
	
	# Write the reference aperture data as a CSV file
	outputFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_reference_aperture.csv"
	outputFile = open(outputFilename, 'w')
	
	lineString = "Frame, Window, X, Y, X_ABS, Y_ABS\n"
	outputFile.write(lineString)
	apertureData = referenceAperture.getData()
	windowIndex = referenceAperture.windowIndex
	window = allWindows[windowIndex]
	xll = window.xll/window.xbin - xmin
	yll = window.yll/window.ybin - ymin
		
	for d in apertureData:
		xAbs = d[2] + xll
		yAbs = d[3] + yll
		lineString = "%d, %d, %f, %f, %f, %f\n"%(d[0], windowIndex, d[2], d[3], xAbs, yAbs)
		outputFile.write(lineString)
	outputFile.close()
	
	# Plot the x-y movement of the reference aperture...
	frameValues = []
	xValues = []
	yValues = []
	for d in apertureData:
		frameValues.append(d[0])
		xValues.append(d[2])
		yValues.append(d[3])
	xStart = xValues[0]
	xValues = [ x - xStart for x in xValues]
	yStart = yValues[0]
	yValues = [ y - yStart for y in yValues]
	xrange = numpy.max(numpy.abs(xValues))
	yrange = numpy.max(numpy.abs(yValues))
	if xrange > yrange:
		plotRange = xrange
	else:
		plotRange = yrange
		
	print "Plot range", plotRange
	
	positionGraph = ppgplot.pgopen('/xs')
	ppgplot.pgenv(numpy.min(frameValues), numpy.max(frameValues), -1.0*plotRange, plotRange,0, 0)
	pgPlotTransform = [0, 1, 0, 0, 0, 1]
	ppgplot.pglab("Frame number", "Pixel position", "Source movement in pixels")
	ppgplot.pgsfs(2)
	ppgplot.pgpt(frameValues, xValues, 2)
	ppgplot.pgpt(frameValues, yValues, 3)
				
	ppgplot.pgclos()
	sys.exit()
	
	
