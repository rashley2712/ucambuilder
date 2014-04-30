#!/usr/bin/env python

import ultracamutils
import matplotlib.pyplot
import argparse
import numpy, math
import classes
import rashley_utils as utils
from trm import ultracam
from trm.ultracam.UErrors import PowerOnOffError, UendError, UltracamError
import sys
import ultracam_shift
import time, datetime
import json
import Image,ImageDraw

def getNextFrame():
	""" Reads the next frame from the CCD using trm routines. Returns all three channels in a dict object (an array of windows)
	    If we are at the first frame, also creates a 'frameInfo' object with information about window dimensions
	"""
	global tempCatalogs
	try:
		ccdFrame = rdat()
	except UltracamError:
		print "There was an error reading the next frame."
		print UltracamError
		return

	frameR = ccdFrame[0]
	frameG = ccdFrame[1]
	frameB = ccdFrame[2]
	
	if (frameIndex == 1):
		# Work out the window structure for this run...
		nxmax = frameR.nxmax
		nymax = frameR.nymax
		numWindows = len(frameR)
		debug.write("Frame dimensions: (%d, %d), Num windows %d"%(nxmax, nymax, numWindows), level = 3)
		frameInfo.nxmax = nxmax
		frameInfo.nymax = nymax
		for nwin,win in enumerate(frameR._data):
			frameInfo.addWindow(win.llx, win.lly, win.nx, win.ny)
			debug.write("Window: %d  llx: %d, lly: %d, nx: %d, ny: %d"%(nwin, win.llx, win.lly, win.nx, win.ny), level = 3)
		""" For the object tracking piece we are going to keep a temporary store of catalogues (of objects) for each window independently.
		    We will use the class FrameCatalogObject for this.
		"""
		#tempCatalogs = classes.FrameCatalogObject(numWindows)
		frameInfo.calcMaxExtents()
		for c in channelNames:
			for w in range(frameInfo.numWindows):
				windowInfo = frameInfo.getWindow(w)
				stackedImageObject = stackedImages[c]
				winnx, winny = (windowInfo.xsize, windowInfo.ysize) 
				blankImage = numpy.zeros((winnx, winny))
				stackedImageObject.addWindow(windowInfo, blankImage)
		

			
   	frameMJD = frameR.time.mjd
	goodTime = frameR.time.good 

	fullFrame = {}
	fullFrame['MJD'] = frameMJD
	frameWindows=[]
	for i in frameR._data:
		frameWindows.append(i.data.T)
	fullFrame["r"] = frameWindows

	frameWindows=[]
	for i in frameG._data:
		frameWindows.append(i.data.T)
	fullFrame["g"] = frameWindows

	frameWindows=[]
	for i in frameB._data:
		frameWindows.append(i.data.T)
	fullFrame["b"] = frameWindows

	return fullFrame

def updateMasterCatalog(newObjects, offSet):
	offSetX , offSetY = offSet
	recognisedCount = 0
	for o in newObjects:
		testX, testY = o['absX'] + offSetX, o['absY'] + offSetY
		
		# First test if this object is close to another....
		objectRecognised = False

		for eo in masterObjectList:
			if (eo.isDistanceMatch(o)!=-1):
				eo.addExposureByObject(o, wholeFrame['MJD'])
				objectRecognised = True
				recognisedCount+= 1
				break
		if (not objectRecognised):
			# Add this to the new object list
			newID = ultracamutils.getUniqueID(masterObjectList)
			newObject = classes.ObservedObject(newID)
			newObject.addExposureByObject(o, wholeFrame['MJD'])
			masterObjectList.append(newObject)
	totalObjects = len(masterObjectList)
	debug.write("%d objects being tracked. %d%% matches in this frame."%(totalObjects, (float(recognisedCount)/float(totalObjects))*100.0))

	
def updateCatalog(MJD, frameNumber, newObjects):
	global prevCatalog
	debug.write("Frame number: %d"%(frameNumber))
	debug.write("Number of objects in this frame: %d"%(len(newObjects)))

	if len(prevCatalog)==0:
		#This is probably the first frame and there are currently no objects to compare to ... add them all to the catalog
		x = []
		y = []
		for o in newObjects:
			x.append(o['absX'])   # Use the x, y values that include the Window offset
			y.append(o['absY'])
		prevCatalog = numpy.array(zip(x, y))
		updateMasterCatalog(newObjects, (0, 0))
		return
		
	# Create a new catalog
	x = []
	y = []
	for o in newObjects:
		x.append(o['absX'])
		y.append(o['absY'])
	newCatalog = numpy.array(zip(x, y))
		
	psize  = 0.5
	fwhm   = 4.
	dmax   = 30.
	mmax   = 30.

	(gaussImage, xp, yp, xr, yr) = ultracam_shift.vimage(prevCatalog, newCatalog, dmax, psize, fwhm)
	(nmatch, inds) = ultracam_shift.match(prevCatalog, newCatalog, xp, yp, mmax)

	prevCatalog = newCatalog

	offsetMag = numpy.sqrt(xr*xr + yr*yr)
	debug.write("Channel: " + channel + " -> Matched objects: %d   Offset distance: %f"%(nmatch, offsetMag))
	
	updateMasterCatalog(newObjects, (xr, yr))
	
	#print inds
	
	if arg.preview:
		matplotlib.pyplot.figure(channel+"_gauss",  figsize=(4,4))
		gaussPlot = matplotlib.pyplot.imshow(gaussImage, cmap=colourMaps[channel], interpolation='nearest')

if __name__ == "__main__":
	
	
	parser = argparse.ArgumentParser(description='Reads the Ultracam [dd-mm-yyyy/runxxx.dat] files and identifies and tracks the objects')
	parser.add_argument('runname', type=str, help='Ultracam run name  [eg 2013-07-21/run010]')
	parser.add_argument('-p', '--preview', action='store_true', help='Show image previews with Matplotlib')
	parser.add_argument('-r', '--crop', action='store_true', help='Crop the images in the preview windows to show only areas with exposed pixels')
	parser.add_argument('-g', '--png', action='store_true', help='Write PNG images of the preview to local disk')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf')
	parser.add_argument('-d', '--debuglevel', type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	parser.add_argument('-n', '--numframes', type=int, help='Number of frames. No parameter means all frames, or from startframe to the end of the run')
	parser.add_argument('-s', '--startframe', default=1, type=int, help='Start frame. \'1\' is the default')
	parser.add_argument('-t', '--sleep', default=0, type=int, help='Sleep time (in seconds) between frames. \'0\' is the default')
	parser.add_argument('-C', '--channels', default='rgb', type=str, help='Which channels to run the extraction on \'rgb\'')
	arg = parser.parse_args()

	channelNames = []
	for channelRequested in arg.channels:
		channelNames.append(channelRequested)
		
	colourMaps = {'r': 'Reds', 'g':'Greens', 'b':'Blues'}
	allObjects = { 'r': [], 'g': [], 'b':[]}

	config = ultracamutils.readConfigFile(arg.configfile)
	debug = classes.debugObject(config.DEBUG)
	debug.toggleTimeLog()
	if (arg.debuglevel!=None): debug.setLevel(arg.debuglevel);
	
	runFilename = utils.addPaths(config.ULTRACAMRAW, arg.runname)

	debug.write("Opening the Ultracam raw file at: " + runFilename, level = 3)
	
	debug.write("Getting run info from the file:" + config.RUNINFO, level = 1)
	runInfo = ultracamutils.getRunInfo(config.RUNINFO, arg.runname)
	
	debug.write("Run Info:\n----------------------", level = 2)
	debug.write(runInfo, level = 2)
	debug.write("----------------------", level = 2)
	
	if arg.preview: 
		matplotlib.pyplot.ion()
		fig = matplotlib.pyplot.gcf()
		
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
	if startFrame>maximumFrames:
		debug.error("startframe " + str(startFrame) + ", is beyond the end of the run, which has only " + str(maximumFrames) + " frames in it.")
		sys.exit()
		
	frameRange = maximumFrames - startFrame + 1
	
	if arg.numframes!=None:
		requestedNumFrames = arg.numframes
		if requestedNumFrames<(frameRange):
			frameRange = requestedNumFrames
	
	frameInfo = classes.FrameObject()
	channelTempCatalogs = {'r':[], 'g':[], 'b':[]}
	stackedImages = {'r':[], 'g':[], 'b':[]}
	""" Created a storage array for keeping the stacked images...
	"""
	for c in channelNames:
		stackedImageStore = classes.stackedImage()
		stackedImages[c] = stackedImageStore
		
	startTime = datetime.datetime.now()
	timeLeftString = "??:??"
	""" Run through all the frames in the .dat file.
	"""
	for frameIndex in range(1, frameRange + 1):
		framesToGo = frameRange - frameIndex
		currentTime = datetime.datetime.now()
		trueFrameNumber = startFrame + frameIndex - 1
		wholeFrame = getNextFrame()
		completionPercent = (float(frameIndex) / float(frameRange) * 100.)
		if frameIndex>2:
			timePassed = ultracamutils.timedeltaTotalSeconds(currentTime - startTime)
			totalTime = timePassed * 100. / completionPercent
			etaTime = startTime + datetime.timedelta(seconds = totalTime)
			timeLeft = etaTime - currentTime
			(hours, mins, secs) = ultracamutils.timedeltaHoursMinsSeconds(timeLeft)
			timeLeftString = str(hours).zfill(2) + ":" + str(mins).zfill(2) + ":" + str(secs).zfill(2)
		
		
		trackingObjectString = ""
		for c in channelNames:
			trackingObjectString+= " " + c + ":" + str(len(allObjects[c]))
		trackingObjectString+= " "
		debug.write(timeLeftString + " Frame: [" + str(frameIndex) + "," + str(trueFrameNumber) +", " + str(framesToGo) + " %d"%(int(completionPercent)) + "%] MJD:" + "%5.7f"%(wholeFrame['MJD'] ) + trackingObjectString + '\r', level = 2)
		
		for channel in channelNames:
			if (channel == 'b') & (trueFrameNumber % rdat.nblue != 0):      # This is an empty blue frame so skip it
				continue   

			prevCatalog = channelTempCatalogs[channel]
			masterObjectList = allObjects[channel]
		
			singleChannelFrame = wholeFrame[channel]
			assembledChannelFrame = numpy.zeros((frameInfo.nxmax, frameInfo.nymax))
			#stackedFrame = numpy.zeros((frameInfo.nxmax, frameInfo.nymax))
			stackedImageObject = stackedImages[channel]
		
			newObjects = []
			for j in range(frameInfo.numWindows): 
				windowImage = singleChannelFrame[j]
				
				tmpFilename = ultracamutils.createFITS(trueFrameNumber, j, channel, windowImage, arg.runname)
				catFilename = ultracamutils.runSex(tmpFilename)
				newObjectsinWindow = ultracamutils.readSexObjects(catFilename)
			
				newObjectsinWindow = ultracamutils.rejectBadObjects(newObjectsinWindow)
			
				if config.KEEP_TMP_FILES!="1":
					ultracamutils.removeTMPFile(tmpFilename)
					ultracamutils.removeTMPFile(catFilename)

				xll = frameInfo.getWindow(j).xll 
				yll = frameInfo.getWindow(j).yll 
				xsize = frameInfo.getWindow(j).xsize 
				ysize = frameInfo.getWindow(j).ysize 
			
				assembledChannelFrame[xll:xll+xsize, yll:yll+ysize] = assembledChannelFrame[xll:xll+xsize, yll:yll+ysize] + ultracamutils.percentiles(windowImage, 20, 98)
				#stackedFrame[xll:xll+xsize, yll:yll+ysize] = stackedFrame[xll:xll+xsize, yll:yll+ysize] + windowImage

				for o in newObjectsinWindow:
					(windowX, windowY) = ( o['y'], o['x'] )
					(absoluteX, absoluteY) = (windowX + xll - 1, windowY + yll - 1)
					o['absX'] = absoluteX
					o['absY'] = absoluteY
					newObjects.append(o)

			if len(newObjects)>0:
				updateCatalog(wholeFrame['MJD'], frameIndex, newObjects)
		
			allObjects[channel] = masterObjectList
			channelTempCatalogs[channel] = prevCatalog
			#stackedImages[channel] = numpy.add(stackedImages[channel], stackedFrame)
			
			for j in range(frameInfo.numWindows):
				windowImage = singleChannelFrame[j]
				stackedImageObject.addNewData(windowImage, j)
			
			if arg.crop:
				xmin, xmax, ymin, ymax = frameInfo.getMaxExtents()
				print "Cropping to:", xmin, xmax, ymin, ymax
				assembledChannelFrame = assembledChannelFrame[xmin:xmax, ymin:ymax]

			if arg.preview:
				# Rotate the image 90 degrees just to make it appear in Matplotlib in the right orientation
				mplFrame = numpy.rot90(assembledChannelFrame)
				#mplFrame = numpy.rot90(stackedImages[channel])
				mplFrame = numpy.flipud(mplFrame)
				fig = matplotlib.pyplot.figure(channel + "_main", figsize=(10,10))
				windowTitle =  "[" + str(trueFrameNumber) + "] " + str(wholeFrame['MJD'])
				fig.canvas.set_window_title(windowTitle)
				imgplot = matplotlib.pyplot.imshow(mplFrame, cmap=colourMaps[channel], interpolation='nearest')
				#matplotlib.pyplot.gca().invert_yaxis()
				for i in newObjects:
					x = i['absX'] if (not arg.crop) else i['absX'] - xmin
					y = i['absY'] if (not arg.crop) else i['absY'] - ymin
					fwhm = i['radius']
					matplotlib.pyplot.gca().add_artist(matplotlib.pyplot.Circle((x,y), fwhm*3, color='green', fill=False, linewidth=2.0))

				if arg.png: 
					matplotlib.pyplot.savefig("r_" + str(frameIndex).zfill(5) +  ".png")
			
			if arg.sleep!=0:
				time.sleep(arg.sleep)

		if arg.preview:
			matplotlib.pyplot.draw()
			matplotlib.pyplot.clf()    # This clears the figure in matplotlib and fixes the 'memory leak'

	""" Write the stacked images to a PNG file
	"""
	for channel in channelNames:
		runIdent = arg.runname
		imageFilename = utils.addPaths(config.SITE_PATH,runIdent) + "_" + channel + ".png"
		stackedImage = stackedImages[channel]
		fullImage =  numpy.zeros((frameInfo.nxmax, frameInfo.nymax))
		for j in range(frameInfo.numWindows):
			windowImage = stackedImage.getWindow(j)
			xll = frameInfo.getWindow(j).xll 
			yll = frameInfo.getWindow(j).yll 
			xsize = frameInfo.getWindow(j).xsize 
			ysize = frameInfo.getWindow(j).ysize 
			fullImage[xll:xll+xsize, yll:yll+ysize] = fullImage[xll:xll+xsize, yll:yll+ysize] + ultracamutils.percentiles(windowImage, 20, 98)
			
		fullImage = numpy.fliplr(fullImage)
		imgData = fullImage
		imgSize = numpy.shape(imgData)
		imgLength = imgSize[0] * imgSize[1]
		testData = numpy.reshape(imgData, imgLength, order="F")
		img = Image.new("L", imgSize)
		img.putdata(testData)
		debug.write("Writing PNG file: " + imageFilename, level = 2) 
		img.save(imageFilename, "PNG")

		
	if (int(config.WRITE_JSON)==1):

		""" First write the raw json files for each channel
		"""
		for channel in channelNames:
			masterObjectList = allObjects[channel]
			allChannelObjects = []
		
			for m in masterObjectList:
				allChannelObjects.append(m.toJSON())
				debug.write(m)
		
			runIdent = arg.runname
		
			outputFilename = utils.addPaths(config.WORKINGDIR,runIdent) 
			outputFilename+= "_" + channel + "_raw.json"
			debug.write("Writing JSON file: " + outputFilename, level = 2)
	
			outputfile = open( outputFilename, "w" )
			json.dump(allChannelObjects, outputfile)
			outputfile.close()

		""" Now write out some meta-data on the run, such as window dimensions, etc
		"""
		outputFilename = utils.addPaths(config.SITE_PATH,runIdent) 
		outputFilename+= "_info.json"
		outputData = {}
		outputData['Target'] = runInfo.target
		xmin, xmax, ymin, ymax = frameInfo.getMaxExtents()
		outputData['numWindows'] = frameInfo.numWindows
		outputData['maxExtents'] = [xmin, xmax, ymin, ymax]
		
		outputfile = open( outputFilename, "w" )
		json.dump(outputData, outputfile)
		outputfile.close()
