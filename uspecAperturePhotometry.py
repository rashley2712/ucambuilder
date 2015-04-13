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
	
	arg = parser.parse_args()
	
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
	
	referenceApertures.getFrameCoverage(frameRange-1)   # The '-1' is because we don't get photometry from the first frame
	
	
	sys.exit()
	
	""" We have run through all of the images now. """
	
	allSources = []
	sourceList = ultraspecClasses.sourceList()
	for index, w in enumerate(allWindows):
		xll = w.xll/w.xbin - xmin
		yll = w.yll/w.ybin - ymin
		image = w.stackedData
		mean, median, std = sigma_clipped_stats(image, sigma=3.0)
		maximum = numpy.max(image)
		minimum = numpy.min(image)
		debug.write("Mean: %f,  Median: %f,  Std (clipped 3sigma):%f"%(mean, median, std) , 2)
		debug.write("Minimum: %f,  Maximum: %f"%(minimum, maximum), 2)
		threshold = median + (std * 2.)
		segm_img = detect_sources(image, threshold, npixels=5)
		mask = segm_img.astype(numpy.bool)
		mean, median, std = sigma_clipped_stats(image, sigma=3.0, mask=mask) 
		debug.write("After source masking", 2)
		debug.write("Mean: %f,  Median: %f,  Std (clipped 3sigma): %f"%(mean, median, std), 2)
		selem = numpy.ones((5, 5))    # dilate using a 5x5 box
		mask2 = binary_dilation(mask, selem)
		mean, median, std = sigma_clipped_stats(image, sigma=3.0, mask=mask2)
		debug.write("After dilation", 2)
		debug.write("Mean: %f,  Median: %f,  Std (clipped 3sigma): %f"%(mean, median, std), 2)
		
		# Check the window image for any areas that should be masked...
		lowerLimitBkg = median - std*5.
		debug.write("5 sigma below the median is the lowerLimitBkg for the mask: %f"%(lowerLimitBkg), 2)
		mask = (image < lowerLimitBkg)
		maskBitmap = numpy.zeros(numpy.shape(mask))
		maskBitmap = 255 * (mask) 
		if (arg.keyimages):
			maskImage = matplotlib.pyplot.figure(figsize=(10, 10))
			matplotlib.pyplot.title("Mask for window:%d"%(index))
			matplotlib.pyplot.imshow(maskBitmap, origin='lower', cmap='Greys_r',  interpolation = 'nearest')
			matplotlib.pyplot.show(block=False)
		bkg = Background(image, (10, 10), filter_shape=(3, 3), method='median', mask=mask)
		background = bkg.background
		if (arg.keyimages):
			bgImage = matplotlib.pyplot.figure(figsize=(10, 10))
			matplotlib.pyplot.title("Fitted background, window:%d"%(index))
			matplotlib.pyplot.imshow(background, origin='lower', cmap='Greys_r',  interpolation = 'nearest')
			matplotlib.pyplot.show(block=False)
		image = image - background
		
		# Final stage source detection
		bkg_sigma = 1.48 * mad(image)
		sigmaThreshold = float(config.SIGMA_THRESHOLD)* float(bkg_sigma)
		debug.write("Threshold for source detection is %f sigma or %f counts."%(float(config.SIGMA_THRESHOLD), sigmaThreshold), 2)
		sources = daofind(image, fwhm=4.0, threshold=sigmaThreshold)   
		
		w.setSourcesAvoidBorders(sources)
		w.BGSubtractedImage = image	
		
		sources = w.getSources()
		for s in sources:
			position = (s['xcentroid'], s['ycentroid'])
			sourceObject = ultraspecClasses.source(0, position, index)
			sourceObject.setDAOPhotData(s['sharpness'], s['roundness1'], s['roundness2'], s['npix'], s['sky'], s['peak'], s['flux'], s['mag'])
			sourceObject.setOffsets((xll, yll))
			sourceList.addSource(sourceObject)
		positions = zip(sources['xcentroid'], sources['ycentroid'], sources['flux'])
		new_positions = [(x + xll, y + yll, flux) for (x, y, flux) in positions]
		allSources+=new_positions
		
	# Get the final stacked image
	fullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
	for w in allWindows:
		imageData = w.stackedData
		boostedImageData = ultracamutils.percentiles(w.BGSubtractedImage, 40, 99.8)
		xll = w.xll/w.xbin - xmin
		xsize = w.nx
		yll = w.yll/w.ybin - ymin
		ysize = w.ny
		fullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + imageData
		boostedFullFrame[yll:yll+ysize, xll:xll+xsize] = boostedFullFrame[yll:yll+ysize, xll:xll+xsize] + boostedImageData

	
	#tempSources = [ (x, y) for (x, y, flux) in allSources]
	#allSources = tempSources	
	allSources = sorted(allSources, key=lambda object: object[2], reverse = True)
	if (arg.keyimages):
		finalFigure = matplotlib.pyplot.figure(figsize=(10, 10))
		matplotlib.pyplot.title("Final stacked image")
		matplotlib.pyplot.imshow(boostedFullFrame, cmap='gray_r')
		for s in allSources:
			(x, y) = s[0], s[1]
			matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), 5, color='blue', fill=False, linewidth=1.0))
		matplotlib.pyplot.gca().invert_yaxis()
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show()
	
	# Output the source list for debug purposes
	if (arg.debuglevel>1):
		sourceString = "Sources"
		for s in allSources:
			sourceString+= "\n(%3.2f, %3.2f) %.2f"%(s[0], s[1], s[2])
		debug.write(sourceString, 2)
		
	
	sourceList.sortByFlux()
	sourcesFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_sources.csv"
	debug.write("Writing source list to CSV file: " + sourcesFilename, 2)
	sourceList.writeToCSV(sourcesFilename)
	
	# Write the XYLS FITS file
	if (arg.xyls):
		IDs = []
		x_values = []
		y_values = []
		fluxes = []
		
		for num, s in enumerate(allSources):
			IDs.append(num)
			x_values.append(s[0])
			y_values.append(s[1])
			fluxes.append(s[2])
			

		FITSFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_sources.xyls"
		debug.write("Writing FITS file: " + FITSFilename, level=2)
		col1 = astropy.io.fits.Column(name='ID', format='I', array=IDs)
		col2 = astropy.io.fits.Column(name='X', format='E', array=x_values)
		col3 = astropy.io.fits.Column(name='Y', format='E', array=y_values)
		col4 = astropy.io.fits.Column(name='FLUX', format='E', array=fluxes)
		cols = astropy.io.fits.ColDefs([col1, col2, col3, col4])	
		#tbhdu =astropy.io.fits.new_table(cols)
		tbhdu =astropy.io.fits.TableHDU.from_columns(cols)
		
		prihdr = astropy.io.fits.Header()
		prihdr['TARGET'] = runInfo.target
		prihdr['RA'] = runInfo.ra
		prihdr['DEC'] = runInfo.dec
		prihdr['COMMENT'] = "This file created by uspecCreateSourceMap.py from the Ultracam pipeline."
		prihdr['RUNIDENT'] = arg.runname
		prihdr['WIDTH'] = fullFramexsize
		prihdr['HEIGHT'] = fullFrameysize
		
		prihdu = astropy.io.fits.PrimaryHDU(header=prihdr)
		thdulist = astropy.io.fits.HDUList([prihdu, tbhdu])
		thdulist.writeto(FITSFilename, clobber=True)
	
	# Generate the stacked image for writing to disc
	stackedFigure = matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Stacked image")
	fullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
	for w in allWindows:
		boostedImage = ultracamutils.percentiles(w.BGSubtractedImage, 10, 99.5)
		xll = w.xll/w.xbin - xmin
		xsize = w.nx
		yll = w.yll/w.ybin - ymin
		ysize = w.ny
		fullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + boostedImage
	
	outputFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + ".png"
	# Write the image data with PIL library, rather than matplotlib
	imgData = numpy.rot90(fullFrame, 3)
	imgSize = numpy.shape(imgData)
	imgLength = imgSize[0] * imgSize[1]
	testData = numpy.reshape(imgData, imgLength, order="F")
	img = Image.new("L", imgSize)
	palette = []
	for i in range(256):
		palette.extend((i, i, i)) # grey scale
		img.putpalette(palette)
	img.putdata(testData)
	debug.write("Writing PNG file: " + outputFilename, level = 2) 
	img.save(outputFilename, "PNG", clobber=True)
	
	palette = []
	for i in range(256):
		palette.extend((255-i, 255-i, 255-i)) # inverse grey scale
		img.putpalette(palette)
	
	outputFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_inverted.png"
	debug.write("Writing PNG file: " + outputFilename, level = 2) 
	img.save(outputFilename, "PNG", clobber=True)
	
	# Write out the stacked image as a non-normalised FITS image
	FITSFilename =  ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_stacked.fits"
	fullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
	for w in allWindows:
		imageData = w.BGSubtractedImage
		xll = w.xll/w.xbin - xmin
		xsize = w.nx
		yll = w.yll/w.ybin - ymin
		ysize = w.ny
		fullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + imageData

	ra = runInfo.ra  # Convert RA to degrees
	dec = runInfo.dec
	fieldScaleX = -8.3E-05
	fieldScaleY = 8.3E-05
	
	prihdr = astropy.io.fits.Header()
	prihdr['COMMENT'] = "This file created by the Ultracam pipeline."
	prihdr['TARGET'] = runInfo.target
	prihdr['COMMENT'] = runInfo.comment
	prihdr['EQUINOX'] = 2000
	prihdr['RADECSYS'] = "FK5"
	prihdr['CTYPE1'] = "RA---TAN"
	prihdr['CTYPE2'] = "DEC--TAN"
	prihdr['CRPIX1'] = fullFramexsize/2
	prihdr['CRPIX2'] = fullFrameysize/2
	prihdr['CRVAL1'] = ra
	prihdr['CRVAL2'] = dec
	prihdr['CDELT1'] = fieldScaleX
	prihdr['CDELT2'] = fieldScaleY
	
	hdu = astropy.io.fits.PrimaryHDU(fullFrame, header=prihdr)
	
	hdulist = astropy.io.fits.HDUList([hdu])
	hdulist.writeto(FITSFilename, clobber=True)
			
	sys.exit()
	
