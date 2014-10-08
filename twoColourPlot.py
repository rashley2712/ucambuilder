#!/usr/bin/env python

import astropy.io.fits
import argparse
import matplotlib.pyplot, numpy, math
import matplotlib.image
import Image, ImageDraw
import ultracamutils
import classes, wcsclasses
import os, subprocess, sys, json

def getArrayFromObjects(objects, propertyName):
	values = []
	for o in objects:
		value = o[propertyName]
		values.append(value)
	return numpy.array(values)
	
def getMeanFlux(photometry, colour):
	totalFlux = 0
	for p in photometry[colour]:
		flux = p[1]
		totalFlux+=flux
	meanFlux = totalFlux/len(photometry[colour])
	
	return meanFlux
	
def filter3colours(objects):
	newObjectList = []
	for o in objects:
		redID = o['colourID']['r']
		greenID = o['colourID']['g']
		blueID = o['colourID']['b']
		if (redID!=-1) & (greenID!=-1) & (blueID!=-1): 
			newObjectList.append(o)
		
	return newObjectList
	
	
if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Plots the 2-colour (u-g, g-r) plot for a run')
	parser.add_argument('runname', type=str, help='Ultracam run name  [eg 2013-07-21/run010]')
	parser.add_argument('-d', '--debuglevel', default = 1, type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf')
	
	arg = parser.parse_args()

	config = ultracamutils.readConfigFile(arg.configfile)
	debug = classes.debugObject(config.DEBUG)
	debug.toggleTimeLog()
	if (arg.debuglevel!=None): debug.setLevel(arg.debuglevel);

	runDate, runNumber = ultracamutils.separateRunNameAndDate(arg.runname)
	
	channels = ['r', 'g', 'b']
	channelDescriptions = {'r':'Red', 'g':'Green', 'b':'Blue'}
	
		
	""" Read info about the run to get the starting RA and DEC coordinates
	"""
	debug.write("Getting run info from the file:" + config.RUNINFO, level = 2)
	runInfo = ultracamutils.getRunInfo(config.RUNINFO, arg.runname)
	
	debug.write("Run Info:\n----------------------", level = 2)
	debug.write(runInfo, level = 2)
	debug.write("----------------------", level = 2)

	
	jsonFilename = ultracamutils.addPaths(config.SITE_PATH, arg.runname) + "_objects.json"
	
	JSONfile = open(jsonFilename, "r")

	wholeFileString = JSONfile.read()

	allObjectsJSON = json.loads(wholeFileString)

	objects = []

	print "Loading JSON file with object data...", jsonFilename
	for i in allObjectsJSON:
		ob = json.loads(i)
		object = {}
		object['id'] = ob['id']
		object['colourID'] = ob['colourID']
		object['meanPosition'] = ob['meanPosition']
		object['photometry'] = ob['photometry']
		objects.append(object)
		print object['id']
		
	
	""" Filter out only those objects with colours in all three channels
	"""
	filteredObjects = filter3colours(objects)
	
	""" Calculate the mean fluxes in each channel
	"""
	for o in filteredObjects:
		meanFlux = {}
		for c in channels:
			meanFlux[c] = getMeanFlux(o['photometry'], c)
		
		o['meanFlux'] = meanFlux
		print o['id'], meanFlux
		uminusg = -2.5 * math.log10(meanFlux['b']/meanFlux['g'])
		gminusr = -2.5 * math.log10(meanFlux['g']/meanFlux['r'])
		
		print uminusg, gminusr
		o['ug'] = uminusg
		o['gr'] = gminusr
	
	figure1 = matplotlib.pyplot.figure(figsize=(12, 12))
	
	x_values = []
	y_values = []
	for o in filteredObjects:
		x_values.append(o['gr'])
		y_values.append(o['ug'])	
	
	print "Number of objects plotted:", len(x_values)
	
	matplotlib.pyplot.plot(x_values, y_values, 'k.')
	matplotlib.pyplot.gca().invert_yaxis()
	matplotlib.pyplot.ylabel('(u-g)')
	matplotlib.pyplot.xlabel('(g-i)')
	

	matplotlib.pyplot.show()
	
	figure1.savefig('test1.eps',dpi=100, format='eps')
	
