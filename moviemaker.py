#!/usr/bin/env python
import astropy, sys, time, os, subprocess, tempfile
import astropy.io.fits, scipy.interpolate
import numpy
from trm import ultracam
import rashley_utils as utils
import classes
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import Image,ImageDraw,ImageFont

def getNextFrame():
	global rdat, frameCounter, frameMJD, nxmax, nymax, numWindows, frameInfo
	try:
		ccdFrame = rdat()
	except:
		print "There was an error reading the next frame."
		return
	
	frameR = ccdFrame[0]
	frameG = ccdFrame[1]
	frameB = ccdFrame[2]

	if (frameCounter == 1):
		# Work out the window structure for this run...
		nxmax = frameR.nxmax
		nymax = frameR.nymax
		numWindows = len(frameR)
		debug.write("(%d, %d), %d"%(nxmax, nymax, numWindows))
		frameInfo.channel = 'r'
		frameInfo.nxmax = nxmax
		frameInfo.nymax = nymax
		for nwin,win in enumerate(frameR._data):
			frameInfo.addWindow(win.llx, win.lly, win.nx, win.ny)
		debug.write(str(frameInfo))

   	timeInfo = rdat.time(frameCounter)
	frameMJD = frameR.time.mjd

	debug.write("Frame:" + str(frameCounter) + "  MJD:" + str(frameMJD), level = 1)
	
	frameWindows = []
	for i in frameR._data:
		frameWindows.append(i.data)
	for i in frameG._data:
		frameWindows.append(i.data)
	for i in frameB._data:
		frameWindows.append(i.data)
	
	return frameWindows

config = utils.readConfigFile()
debug = classes.debugObject(1)
debug.toggleTimeLog()

if (len(sys.argv) < 2):
	print "Please give me a run name."
	sys.exit()

runName = sys.argv[1]
runFilename = utils.addPaths(config.ULTRACAMRAW, runName)

startFrame = 1
requestedNumFrames = 20
CCDside = 0
keepTmpFiles = False
frameInfo = classes.FrameObject()
tmpMoviePath = config.MOVIE_TMP_PATH
movieFilename = utils.addPaths(config.SITE_PATH, runName)
movieFilename+= "movie.mp4"

for i in sys.argv[2:]:
	if i[:2]=="-n": 
		debug.write("command line: number of frames:" + i[2:]) 
		framesStr = i[2:]
		requestedNumFrames = int(framesStr)
	if i[:2]=="-d": 
		debug.write("command line: debugLevel:" + i[2:]) 
		debugLevel = int(i[2:])
		debug.setLevel(debugLevel)
	if i[:2]=="-s": 
		print "starting at frame:", i[2:]
		framesStr = i[2:]
		startFrame = int(framesStr)
	if i[:2]=="-o": 
		print "output file:",  i[2:] 
		outputFilename = i[2:]
	if i[:2]=="-c": 
		if (i[2:]=='r'): 
			CCDside = 1
			print "CCD right side" 
	if i[:2]=="-k":
		keepTmpFiles = True


date = utils.separateRunNameAndDate(runName)[0]
newFolder = utils.addPaths(config.MOVIE_TMP_PATH, date)
print "Creating folder:", newFolder
subprocess.call(["mkdir", newFolder])

rdat  = ultracam.Rdata(runFilename, startFrame, server=False)
numFrames = rdat.ntotal()
frameMJD = 0

debug.write("Total frames in %s: %d"%(runName,numFrames))
debug.write("nBlue: " + str(rdat.nblue))

frameCounter = 0
colours = ['r','g','b']

for i in range(requestedNumFrames):
	if (frameCounter<numFrames):
		frameCounter+= 1
	else: 
		break

	frameWindows = getNextFrame()
	rgbFrames = []
	
	for i, colour in enumerate(colours):
	
		# Construct full frame from Windows
		fullFrame = numpy.zeros((frameInfo.nxmax, frameInfo.nymax))

		for j in range(frameInfo.numWindows): 
	
			frameImage = frameWindows[j + frameInfo.numWindows*i]
			frameImage = numpy.rot90(frameImage)
					
			# All windows need to be rotated by 180 degrees before using the Image library to write to png
			xll = frameInfo.getWindow(j).xll 
			yll = frameInfo.getWindow(j).yll 
			xsize = frameInfo.getWindow(j).xsize 
			ysize = frameInfo.getWindow(j).ysize 
			rotatedImage = numpy.flipud(frameImage)
			normalisedImage = utils.percentiles(rotatedImage, 25, 99)
			fullFrame[xll:xll+xsize, yll:yll+ysize] = fullFrame[xll:xll+xsize, yll:yll+ysize] + normalisedImage
	
		fullFrame = numpy.fliplr(fullFrame)
	
		imgData = fullFrame
		imgSize = numpy.shape(imgData)
		imgLength = imgSize[0] * imgSize[1]
		testData = numpy.reshape(imgData, imgLength, order="F")
		img = Image.new("L", imgSize)
		img.putdata(testData)
			
		if (colour=="b"):
			if ((frameCounter % rdat.nblue)!= 0):    
				debug.write("Skipping an empty blue frame...")
				if (frameCounter==1): 
					# This is the first frame, so it is ok to just use the first blue frame, which is readout noise
					blueChannelHolder = img.copy()
				else:
					# Use last good blue image
					img = blueChannelHolder
			else:
				# This was a 'real' blue image so use it and keep it
				blueChannelHolder = img.copy()

		if (keepTmpFiles==True):
			outputFilename = utils.addPaths(config.MOVIE_TMP_PATH,runName)  + str(frameCounter).zfill(5) + colour + "."
			caption = str(frameCounter).zfill(5) + " : " + str(frameMJD) + " [" + colour + "]"
			utils.writePNG(img, outputFilename + "png", caption)
	
		rgbFrames.append(img)

	img = Image.merge("RGB", (rgbFrames[0],rgbFrames[1],rgbFrames[2]))			
	
	caption = str(frameCounter).zfill(5) + " : " + str(frameMJD) + " [rgb]"
	
	outputFilename = utils.addPaths(config.MOVIE_TMP_PATH,runName) + str(frameCounter).zfill(5) + "rgb"
	utils.writePNG(img, outputFilename, caption)
		

print "Written %d png files to %s folder"%(frameCounter, config.MOVIE_TMP_PATH)

pngPath = utils.addPaths(config.MOVIE_TMP_PATH, date)
pngFiles = pngPath + "/*.png"
subprocess.call(["mencoder", "mf://" + pngFiles, "-o", movieFilename, "-fps", "10", "-ovc", "lavc"])

subprocess.call(["rm", "-f", pngPath + "/*"])



