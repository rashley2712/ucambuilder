#!/usr/bin/env python

import ultracamutils
import matplotlib.pyplot
import argparse
import numpy
import classes
import rashley_utils as utils
from trm import ultracam
from trm.ultracam.UErrors import PowerOnOffError, UendError, UltracamError
import sys

def getNextFrame():
	""" Reads the next frame from the CCD using trm routines. Returns all three channels in a dict object (an array of windows)
	    If we are at the first frame, also creates a 'frameInfo' object with information about window dimensions
	"""
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
		    We will create a global variable, dictionary called windowCat to do this.
		"""
		for i in range(frameInfo.numWindows):
			window = frameInfo.getWindow(i)
			 

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

	
def updateCatalog(MJD, windowNumber, frameNumber, newObjects):
	debug.write("Updating the catalog...")
	debug.write("Frame number: %d, Window number: %d"%(frameNumber, windowNumber))
	debug.write("Number of objects in this window: %d"%(len(newObjects)))
	debug.write("Number of objects in the same window of the previous frame: %d"%(len(oldcat)))
	#for o in newObjects:
	#	print o
		


if __name__ == "__main__":
	
	
	parser = argparse.ArgumentParser(description='Reads the Ultracam [dd-mm-yyyy/runxxx.dat] files and identifies and tracks the objects')
	parser.add_argument('runname', type=str, help='Ultracam run name  [eg 2013-07-21/run010]')
	parser.add_argument('-p', '--preview', action='store_true', help='Show image previews with Matplotlib')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf')
	parser.add_argument('-d', '--debuglevel', type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	parser.add_argument('-n', '--numframes', type=int, help='Number of frames. No parameter means all frames, or from startframe to the end of the run')
	parser.add_argument('-s', '--startframe', default=1, type=int, help='Start frame. \'1\' is the default')
	arg = parser.parse_args()

	config = ultracamutils.readConfigFile(arg.configfile)
	debug = classes.debugObject(config.DEBUG)
	debug.toggleTimeLog()
	if (arg.debuglevel!=None): debug.setLevel(arg.debuglevel);
	
	runFilename = utils.addPaths(config.ULTRACAMRAW, arg.runname)

	debug.write("Opening the Ultracam raw file at: " + runFilename, level = 3)
	
	if arg.preview: 
		matplotlib.pyplot.ion()
	
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
	masterObjectList = []
	""" Run through all the frames in the .dat file.
	"""
	for frameIndex in range(1, frameRange + 1):
		trueFrameNumber = startFrame + frameIndex - 1
		wholeFrame = getNextFrame()
		debug.write("Frame: [" + str(frameIndex) + "," + str(trueFrameNumber) + "] MJD:" + str(wholeFrame['MJD']), level = 2)
		
		redFrame = wholeFrame['r']
		assembledRedFrame = numpy.zeros((frameInfo.nxmax, frameInfo.nymax))
		
		for j in range(frameInfo.numWindows): 
			windowImage = redFrame[j]
			tmpFilename = ultracamutils.createFITS(trueFrameNumber, j, 'r', windowImage)
			catFilename = ultracamutils.runSex(tmpFilename)
			newObjects = ultracamutils.readSexObjects(catFilename)
			
			updateCatalog(wholeFrame['MJD'], j, frameIndex, newObjects)

			if config.KEEP_TMP_FILES!="1":
				ultracamutils.removeTMPFile(tmpFilename)
				ultracamutils.removeTMPFile(catFilename)

			xll = frameInfo.getWindow(j).xll 
			yll = frameInfo.getWindow(j).yll 
			xsize = frameInfo.getWindow(j).xsize 
			ysize = frameInfo.getWindow(j).ysize 
			
			assembledRedFrame[xll:xll+xsize, yll:yll+ysize] = assembledRedFrame[xll:xll+xsize, yll:yll+ysize] + ultracamutils.percentiles(windowImage, 20, 98)

		if arg.preview:
			# Rotate the image 90 degrees just to make it appear in Matplotlib in the right orientation
			mplFrame = numpy.rot90(assembledRedFrame)
			imgplot = matplotlib.pyplot.imshow(mplFrame, cmap='gray', interpolation='nearest')
			matplotlib.pyplot.draw()
			matplotlib.pyplot.clf()    # This clears the figure in matplotlib and fixes the 'memory leak'
	