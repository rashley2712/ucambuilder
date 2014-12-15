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
	
def distance(p1, p2):
	xSep = p1[0] - p2[0]
	ySep = p1[1] - p2[1]
	return math.sqrt( (xSep*xSep) + (ySep*ySep) )
	
def getMeanColourPosition(photometry, colour):
	
	xMean = 0
	yMean = 0
	for p in photometry[colour]:
		xMean+= p[3]
		yMean+= p[4]
	xMean = xMean/len(photometry[colour])
	yMean = yMean/len(photometry[colour])
		
	return (xMean,yMean)
	
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
	parser.add_argument('-n', '--numobjects', default=-1, type=int, help='Limit the number of objects to load for the plot. Default is ''all''.')
	parser.add_argument('--labels', action = 'store_true', help='Draw labels for each object plotted.')
	
	
	arg = parser.parse_args()

	if arg.numobjects!=-1:
		limitObjects = True
		numObjects = arg.numobjects
	else:
		limitObjects = False
		
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
	for index, i in enumerate(allObjectsJSON):
		ob = json.loads(i)
		object = {}
		object['id'] = ob['id']
		object['colourID'] = ob['colourID']
		object['meanPosition'] = ob['meanPosition']
		object['photometry'] = ob['photometry']
		objects.append(object)
		print object['id']
		if limitObjects:
			if index>=(numObjects-1): break;
		
	""" Filter out only those objects with colours in all three channels
	"""
	filteredObjects = filter3colours(objects)
	
	""" Calculate the mean fluxes in each channel
	"""
	gmagMean = 0
	for o in filteredObjects:
		meanFlux = {}
		meanColourPosition = {}
		for c in channels:
			meanFlux[c] = getMeanFlux(o['photometry'], c)
			meanColourPosition[c] = getMeanColourPosition(o['photometry'], c)
			
		o['meanPosition'] = meanColourPosition
		o['rg_distance'] = distance(meanColourPosition['r'], meanColourPosition['g'])	
		o['rb_distance'] = distance(meanColourPosition['r'], meanColourPosition['b'])	
		o['3D_distance'] = math.sqrt( o['rg_distance']*o['rg_distance'] + o['rb_distance']*o['rb_distance'])
		
		o['meanFlux'] = meanFlux
		print o['id'], meanFlux, o['3D_distance']
		
		uminusg = -2.5 * math.log10(meanFlux['b']/meanFlux['g'])
		gminusr = -2.5 * math.log10(meanFlux['g']/meanFlux['r'])
		gmag = -2.5 * math.log10(meanFlux['g'])
		
		#print uminusg, gminusr
		o['ug'] = uminusg
		o['gr'] = gminusr
		o['gmag'] = gmag
		gmagMean+=gmag
	
	
	# Remove the 'g' offset
	gmagMean = gmagMean/len(filteredObjects)
	for o in filteredObjects:
		o['gmag'] = o['gmag'] - gmagMean
		
	
	figure1 = matplotlib.pyplot.figure(figsize=(12, 12))
	
	x_values = []
	y_values = []
	distances = []
	labels = []
	for o in filteredObjects:
		x_values.append(o['gr'])
		y_values.append(o['ug'])	
		distances.append(o['3D_distance'])
		labels.append(o['id'])
	minDistance = numpy.min(distances)
	maxDistance = numpy.max(distances)
	rangeDistances = maxDistance - minDistance
	print minDistance, maxDistance, rangeDistances
	distances = [ (d - minDistance) / rangeDistances for d in distances]
	print distances
	
	print "Number of objects plotted:", len(x_values)
	
	matplotlib.pyplot.scatter(x_values, y_values, c=distances, cmap=matplotlib.pyplot.cm.coolwarm, s=50)
	
	matplotlib.pyplot.gca().invert_yaxis()
	plot = matplotlib.pyplot.gcf().add_subplot(111)
	plot.tick_params(axis='both', which='major', labelsize=14)
	matplotlib.pyplot.ylabel('(u-g)', size = 18)
	matplotlib.pyplot.xlabel('(g-i)', size = 18)
	
	
	#legend
	#fig, ax = matplotlib.pyplot.subplots()
	#heatmap = ax.pcolor(distances, cmap=matplotlib.pyplot.cm.coolwarm)
	cbar = matplotlib.pyplot.colorbar(shrink = 0.7)
	ticks = numpy.linspace(minDistance, maxDistance, 10)
	print ticks
	tickLabels = [str(round(t*10)/10) for t in ticks]
	print tickLabels
	cbar.ax.set_yticklabels(tickLabels)
	cbar.ax.get_yaxis().labelpad = 20
	cbar.set_label('separation distance', rotation=270, size=16)

	if (arg.labels):
		for i, txt in enumerate(labels):
			matplotlib.pyplot.annotate(txt, (x_values[i]+0.01, y_values[i]-0.01), textcoords='offset points')


	matplotlib.pyplot.show()
	
	figure1.savefig('colour-colour.eps',dpi=100, format='eps')
	
	# Now try a colour-magnitude plot
	for o in filteredObjects:
		meanFlux = o['meanFlux']
		print meanFlux
		
	figure2 = matplotlib.pyplot.figure(figsize=(12, 12))
	
	x_values = []
	y_values = []
	labels = []
	for o in filteredObjects:
		if o['id']!=782:
			x_values.append(o['gr'])
			y_values.append(o['gmag'])	
			labels.append(o['id'])
	
	print "Number of objects plotted:", len(x_values)
	matplotlib.pyplot.plot(x_values, y_values, 'ko')
	matplotlib.pyplot.gca().invert_yaxis()
	plot = matplotlib.pyplot.gcf().add_subplot(111)
	plot.tick_params(axis='both', which='major', labelsize=12)
	matplotlib.pyplot.ylabel('g', size = 18)
	matplotlib.pyplot.xlabel('g-r', size = 18)
	if (arg.labels):
		for i, txt in enumerate(labels):
			matplotlib.pyplot.annotate(txt, (x_values[i]+0.01, y_values[i]-0.01), textcoords='offset points')

	matplotlib.pyplot.show()
	
	figure2.savefig('colour-magnitude.eps',dpi=100, format='eps')
	
