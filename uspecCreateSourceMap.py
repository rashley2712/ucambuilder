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
	
	print "Astropy version", astropy.__version__
	
	parser = argparse.ArgumentParser(description='Reads the Ultraspec [dd-mm-yyyy/runxxx.dat] files produces previews of the images')
	parser.add_argument('runname', type=str, help='Ultracam run name  [eg 2013-07-21/run010]')
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
	
	arg = parser.parse_args()
	
	
	config = ultracamutils.readConfigFile(arg.configfile)
	
	applyShift = True
	if (arg.noshift):
		applyShift = False
	debug = classes.debugObject(config.DEBUG)
	debug.toggleTimeLog()
	if (arg.debuglevel!=None): debug.setLevel(arg.debuglevel);
	debug.write(arg)
	
	runInfo = ultraspecClasses.runInfo(arg.runname)
	found = runInfo.loadFromJSON(config.RUNINFO)
	if not found:
		debug.write("Could not get info for this run from the ultra.json file.", level = 1)
		xmlRead = runInfo.loadFromXML(config.ULTRASPECRAW)
		
	runFilename = ultracamutils.addPaths(config.ULTRASPECRAW, arg.runname)

	debug.write("Opening the Ultraspec raw file at: " + runFilename, level = 3)
	
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
		
	if maximumFrames<10:
		print "The total number of frames in this run is less than 10. We need more to create a source map. Exiting."
		sys.exit()
	
	""" We are going to use the first 10 frames to build the original source map that is merely used as a guide for creating the main stacked image.
	"""	
	
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
		if (arg.usefirstframe):
			window.setData(image)
			bkg_sigma = 1.48 * mad(image)
			print "bkg_sigma", bkg_sigma   
			sources = daofind(image, fwhm=4.0, threshold=3*bkg_sigma)   
			window.setSources(sources)	
		else: 
			window.setBlankData(image)
		
		
		allWindows.append(window)
		
	(xmin, ymin, xmax, ymax) = determineFullFrameSize(allWindows)
	fullFramexsize = xmax - xmin
	fullFrameysize = ymax - ymin
	sourceMap = ultraspecClasses.sourceMap((fullFrameysize, fullFramexsize))
	
	for frame in range(10):
		ccdFrame = rdat()
		windows = ccdFrame[0]
		
		for windowIndex, w in enumerate(windows):
			image = w._data
			allWindows[windowIndex].addData(image)
			
		
	# Reconstruct the full frame from the windows	
	stackedFigure = matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Initial 10 frame stacked image")
	boostedFullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
	fullFrame = numpy.zeros((fullFrameysize, fullFramexsize))	
	for w in allWindows:
		boostedImage = ultracamutils.percentiles(w.stackedData, 10, 99.8)
		image = w.stackedData
		xll = w.xll/w.xbin - xmin
		xsize = w.nx
		yll = w.yll/w.ybin - ymin
		ysize = w.ny
		boostedFullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + boostedImage
		fullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + image
		
		bkg_sigma = 1.48 * mad(image)
		print "bkg_sigma", bkg_sigma   
		sources = daofind(image, fwhm=4.0, threshold=3*bkg_sigma) 
		print "Num sources in this window:", len(sources)  
		w.setSourcesAvoidBorders(sources)	
		
		
	# Get the source list from this image
	# Combine the sources in all of the windows
	allSources = []
	for index, w in enumerate(allWindows):
		xll = w.xll/w.xbin - xmin
		yll = w.yll/w.ybin - ymin
		sources = w.getSources()
		positions = zip(sources['xcentroid'], sources['ycentroid'], sources['flux'])
		new_positions = [(x + xll, y + yll, flux) for (x, y, flux) in positions]
		allSources+=new_positions
		
		
	allSources = sorted(allSources, key=lambda object: object[2], reverse = True)
	numSources = len(allSources)
	maxSources = int(round((numSources)*0.4))
	topSources = allSources[0:maxSources]
	print "Number of sources: %d, number of top sources: %d"%(numSources, maxSources)
	# Display the image on the user's screen
	image = matplotlib.pyplot.imshow(boostedFullFrame, cmap='gray_r')
	for s in allSources:
		x, y = s[0], s[1]
		matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), 10, color='green', fill=False, linewidth=1.0))
	for s in topSources:
		x, y = s[0], s[1]
		matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), 10, color='blue', fill=False, linewidth=1.0))
	
	matplotlib.pyplot.gca().invert_yaxis()			
	matplotlib.pyplot.show(block=False)

	masterApertureList = [ (x, y) for (x, y, flux) in topSources]
		
	""" End of the prework """

	rdat.set(1)		# Reset back to the first frame
	
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
		zoomedImage = matplotlib.pyplot.figure(figsize=(5, 5))
		
			
	fullFrame = numpy.zeros((1057, 1040))
			
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
		#image -= numpy.median(image)
		if (arg.usefirstframe):
			window.setData(image)
			bkg_sigma = 1.48 * mad(image)
			print "bkg_sigma", bkg_sigma   
			sources = daofind(image, fwhm=4.0, threshold=3*bkg_sigma)   
			window.setSourcesAvoidBorders(sources)	
		else: 
			window.setBlankData(image)
		
		
		allWindows.append(window)
		
	(xmin, ymin, xmax, ymax) = determineFullFrameSize(allWindows)
	fullFramexsize = xmax - xmin
	fullFrameysize = ymax - ymin
	sourceMap = ultraspecClasses.sourceMap((fullFrameysize, fullFramexsize))
			
	debug.write("Building a map of sources in order to define the apertures...", level = 2)
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
			#image -= numpy.median(image)
			#allWindows[windowIndex].addData(image)
			bkg_sigma = 1.48 * mad(image)
			sources = daofind(image, fwhm=4.0, threshold=3*bkg_sigma)   
			#print sources
			allWindows[windowIndex].setSourcesAvoidBorders(sources)	
			
		# Combine the sources in all of the windows
		allSources = []
		for index, w in enumerate(allWindows):
			xll = w.xll/w.xbin - xmin
			yll = w.yll/w.ybin - ymin
			sources = w.getSources()
			positions = zip(sources['xcentroid'], sources['ycentroid'], sources['flux'])
			new_positions = [(x + xll, y + yll, flux) for (x, y, flux) in positions]
			allSources+=new_positions

		allSources = sorted(allSources, key=lambda object: object[2], reverse = True)
		#print allSources
		tempSources = [ (x, y) for (x, y, flux) in allSources]
		#print tempSources
		allSources = tempSources
		sourceMap.updateMap(allSources)
		
		if (applyShift):
		
			oldCatalog = numpy.array(masterApertureList)
			newCatalog = numpy.array(allSources)

			psize  = 0.1
			fwhm   = 4.
			dmax   = 10.
			mmax   = 10.

			(gaussImage, xp, yp, xr, yr) = ultracam_shift.vimage(oldCatalog, newCatalog, dmax, psize, fwhm)
			# (nmatch, inds) = ultracam_shift.match(oldCatalog, newCatalog, xp, yp, mmax)
			debug.write("Calculated offset: (%2.2f, %2.2f)"%(xr, yr), level = 2)

		for windowIndex, w in enumerate(windows):
			image = w._data
			allWindows[windowIndex].setData(image)
			if (applyShift): image = ndimage.interpolation.shift(image, (-1.0*yr, -1.0*xr), order = 3 )
			allWindows[windowIndex].addToStack(image)
			
			
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
					
			matplotlib.pyplot.figure(fig.number)
			matplotlib.pyplot.imshow(fullFrame, cmap='gray_r')
			
			for s in allSources:
				(x, y) = s
				matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), 15, color='green', fill=False, linewidth=1.0))
			
			
			matplotlib.pyplot.title("Frame image [%d/%d]"%(trueFrameNumber, frameRange))
			if arg.stack:
				matplotlib.pyplot.title("Stacked image [%d/%d]"%(trueFrameNumber, frameRange))
			matplotlib.pyplot.gca().invert_yaxis()
			matplotlib.pyplot.draw()
			matplotlib.pyplot.show()
			matplotlib.pyplot.clf()    # This clears the figure in matplotlib and fixes the 'memory leak'
			
			# Now also draw the zoomed in region around the first aperture
			(x, y) = masterApertureList[0]
			croppedFrame = fullFrame[y-10:y+10, x-9:x+10]
			
			matplotlib.pyplot.figure(zoomedImage.number)
			matplotlib.pyplot.imshow(croppedFrame, cmap='gray_r')
			matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((10,10), 1, color='green', fill=False, linewidth=1.0))
			matplotlib.pyplot.plot([10, 10+xr], [ 10, 10+yr], lw=1, color='green')
			matplotlib.pyplot.title("Zoom on aperture number 1: Frame [%d/%d]"%(trueFrameNumber, frameRange))
			matplotlib.pyplot.xlim([0, 20])
			matplotlib.pyplot.ylim([0, 20])
			matplotlib.pyplot.gca().invert_yaxis()
			matplotlib.pyplot.draw()
			matplotlib.pyplot.show()
			
			matplotlib.pyplot.clf()    # This clears the figure in matplotlib and fixes the 'memory leak'
			
		if arg.sleep!=0:
			time.sleep(arg.sleep)
	sys.stdout.write("\rProcessed %d frames      \n"%frameRange)
	sys.stdout.flush()
	
	
	
	
	allSources = []
	for index, w in enumerate(allWindows):
		xll = w.xll/w.xbin - xmin
		yll = w.yll/w.ybin - ymin
		image = w.stackedData
		mean, median, std = sigma_clipped_stats(image, sigma=3.0)
		maximum = numpy.max(image)
		minimum = numpy.min(image)
		print "Mean:", mean, " Median:", median, " Std (clipped 3sigma):", std
		print "Minimum:", minimum, " Maximum:", maximum
		threshold = median + (std * 2.)
		segm_img = detect_sources(image, threshold, npixels=5)
		mask = segm_img.astype(numpy.bool)
		mean, median, std = sigma_clipped_stats(image, sigma=3.0, mask=mask) 
		print "After source masking"
		print "Mean:", mean, " Median:", median, " Std (clipped 3sigma):", std
		selem = numpy.ones((5, 5))    # dilate using a 5x5 box
		mask2 = binary_dilation(mask, selem)
		mean, median, std = sigma_clipped_stats(image, sigma=3.0, mask=mask2)
		print "After dilation"
		print "Mean:", mean, " Median:", median, " Std (clipped 3sigma):", std
		bkg = Background(image, (25, 25), filter_shape=(3, 3), method='median')
		bgImage = matplotlib.pyplot.figure(figsize=(10, 10))
		matplotlib.pyplot.imshow(bkg.background, origin='lower', cmap='Greys_r')
		matplotlib.pyplot.show(block=False)
		image = image - bkg.background
		bkg_sigma = 1.48 * mad(image)
		sources = daofind(image, fwhm=4.0, threshold=3*bkg_sigma)   
		print sources
		w.setSourcesAvoidBorders(sources)
		w.BGSubtractedImage = image	
		
		sources = w.getSources()
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

	
	tempSources = [ (x, y) for (x, y, flux) in allSources]
	allSources = tempSources	
	finalFigure = matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Final stacked image")
	matplotlib.pyplot.imshow(boostedFullFrame, cmap='gray_r')
	for s in allSources:
		(x, y) = s
		matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), 5, color='blue', fill=False, linewidth=1.0))
	matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.draw()
	matplotlib.pyplot.show()
	
	# Get the source map
	smoothedSourceMap = sourceMap.getSmoothMap()
	
	# Now use this source map to generate a set of apertures
	bkg_sigma = 1.48 * mad(smoothedSourceMap)
	print "sourceMap median:", numpy.median(smoothedSourceMap)
	print "sourceMap mean:", numpy.mean(smoothedSourceMap)
	print "sourceMap max:", numpy.max(smoothedSourceMap)
	threshold = frameRange/100.
	print "threshold:", threshold
	apertureSources = daofind(smoothedSourceMap, fwhm=4.0, threshold=threshold)   
	
	# Draw the source map 
	sourceMapImage = matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Source map")
	matplotlib.pyplot.imshow(smoothedSourceMap, cmap='hot')
	for s in apertureSources:
		x, y = s['xcentroid'], s['ycentroid']
		matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), 10, color='green', fill=False, linewidth=1.0))
	matplotlib.pyplot.gca().invert_yaxis()			
	#matplotlib.pyplot.show(block=False)
	outputFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_sourcemap.png"
	matplotlib.pyplot.savefig(outputFilename)

	# Draw the source map with no apertures
	sourceMapImage = matplotlib.pyplot.figure(figsize=(10, 10))
	matplotlib.pyplot.title("Source map")
	matplotlib.pyplot.imshow(smoothedSourceMap, cmap='hot')
	matplotlib.pyplot.gca().invert_yaxis()
	outputFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_sourcemap_clean.png"
	matplotlib.pyplot.savefig(outputFilename)


	# Write the XYLS FITS file
	if (arg.xyls):
		IDs = []
		x_values = []
		y_values = []
		fluxes = []
		sortedObjects = sorted(apertureSources, key= lambda p:p['flux'], reverse=True)
		
		for num, s in enumerate(sortedObjects):
			IDs.append(num)
			x_values.append(s['xcentroid'])
			y_values.append(s['ycentroid'])
			fluxes.append(s['flux'])
			

		FITSFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_sources.xyls"
		debug.write("Writing FITS file: " + FITSFilename, level=2)
		col1 = astropy.io.fits.Column(name='ID', format='I', array=IDs)
		col2 = astropy.io.fits.Column(name='X', format='E', array=x_values)
		col3 = astropy.io.fits.Column(name='Y', format='E', array=y_values)
		col4 = astropy.io.fits.Column(name='FLUX', format='E', array=fluxes)
		cols = astropy.io.fits.ColDefs([col1, col2, col3, col4])	
		tbhdu =astropy.io.fits.new_table(cols)
		
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
		boostedImage = ultracamutils.percentiles(w.stackedData, 40, 99.8)
		xll = w.xll/w.xbin - xmin
		xsize = w.nx
		yll = w.yll/w.ybin - ymin
		ysize = w.ny
		fullFrame[yll:yll+ysize, xll:xll+xsize] = fullFrame[yll:yll+ysize, xll:xll+xsize] + boostedImage
	
	outputFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + ".png"
	
	# Write the image with PIL library, rather than matplotlib
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
		imageData = w.stackedData
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
	
	# Now write out the aperture data
	
	outputFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_apertures.json"
	outputFile = open(outputFilename, "w")
	apertureList = []
	for i, s in enumerate(apertureSources):
		aperture = {}
		aperture['id'] = i
		aperture['x'] = s['xcentroid']
		aperture['y'] = s['ycentroid']
		aperture['sharpness'] = s['sharpness']
		aperture['roundness1'] = s['roundness1']
		aperture['roundness2'] = s['roundness2']
		aperture['flux'] = s['flux']
		apertureList.append(aperture)
	json.dump(apertureList, outputFile)
	outputFile.close()
		

