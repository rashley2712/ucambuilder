#!/usr/bin/env python

import astropy.io.fits
import argparse
import matplotlib.pyplot, numpy
import Image, ImageDraw
import ultracamutils
import classes, wcsclasses
import os, subprocess

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Combines the catalog files (produced by the "postprocessor.py" program and uses "solve-field" to find WCS solutions')
	parser.add_argument('runname', type=str, help='Ultracam run name  [eg 2013-07-21/run010]')
	parser.add_argument('-d', '--debuglevel', default = 1, type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf')
	parser.add_argument('-f', '--forcesolve', action='store_true', help='Force a "solve" even if we already have a solution')
	parser.add_argument('-p', '--preview', action='store_true', help='Show a preview of the check plots.')

	arg = parser.parse_args()

	config = ultracamutils.readConfigFile(arg.configfile)
	debug = classes.debugObject(config.DEBUG)
	debug.toggleTimeLog()
	if (arg.debuglevel!=None): debug.setLevel(arg.debuglevel);

	""" First check if a directory is made for the output files in the working director folder. and create one. 
	"""
	runDate, runNumber = ultracamutils.separateRunNameAndDate(arg.runname)
	
	outputDirectory = ultracamutils.addPaths(config.WORKINGDIR, runDate)
	
	if not os.path.exists(outputDirectory):
		debug.write("Creating the directory: " + outputDirectory, level = 2)
		os.mkdir(outputDirectory)

	
	""" Produce a 'preview plot' of the catalogs and the png files 
	"""
	channels = ['r', 'g', 'b']
	channelDescriptions = {'r':'Red', 'g':'Green', 'b':'Blue'}
	for c in channels:
		pngFile = ultracamutils.addPaths(config.SITE_PATH, arg.runname) + '_' + c + ".png"
		xylsFile = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + '_' + c + ".xyls"
		checkplotFile = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + '_' + c + "_preview.png"

		plotXYLSCommand = ["plotXYLS.py"]
		plotXYLSCommand.append("-i" + xylsFile)
		plotXYLSCommand.append("-m" + pngFile)
		plotXYLSCommand.append("-o" + checkplotFile)
		if (arg.preview): plotXYLSCommand.append("-p")

		subprocess.call(plotXYLSCommand)
		
	""" Read info about the run to get the starting RA and DEC coordinates
	"""
	debug.write("Getting run info from the file:" + config.RUNINFO, level = 2)
	runInfo = ultracamutils.getRunInfo(config.RUNINFO, arg.runname)
	
	debug.write("Run Info:\n----------------------", level = 2)
	debug.write(runInfo, level = 2)
	debug.write("----------------------", level = 2)

	solved = [False, False, False]
	for n, c in enumerate(channels):
		solvedFileMarker = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + '_' + c + ".solved"
		if os.path.exists(solvedFileMarker):
			solved[n] = True

	""" solve-field -X X -Y Y test.xyls --width 1024 --height 1024 --ra 296.007 --dec 40.2954 --radius 1 --overwrite -L 5 -u amw
	"""
	for n, c in enumerate(channels):
		if (not solved[n]) or (arg.forcesolve):
			solvefieldCommand = ["solve-field"]
			solvefieldCommand.append(ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + '_' + c + ".xyls")
			solvefieldCommand.append("-XX")
			solvefieldCommand.append("-YY")
			solvefieldCommand.append("-w1024")
			solvefieldCommand.append("-e1024")
			solvefieldCommand.append("--overwrite")
			solvefieldCommand.append("-L5")
			solvefieldCommand.append("-uamw")
			solvefieldCommand.append("--ra")
			solvefieldCommand.append(str(runInfo.ra*15.))
			solvefieldCommand.append("--dec")
			solvefieldCommand.append(str(runInfo.dec))
			solvefieldCommand.append("--radius")
			solvefieldCommand.append(str(1))
			solvefieldCommand.append("-t3")
			solvefieldCommand.append("--no-plots")
	
			subprocess.call(solvefieldCommand)
		else:
			debug.write("Skipping the %s channel as it appears to be solved"%(channelDescriptions[c]), level = 2)
	
	""" Check if we have a solution for each channel
	"""
	solved = [False, False, False]
	for n, c in enumerate(channels):
		solvedFileMarker = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + '_' + c + ".solved"
		if os.path.exists(solvedFileMarker):
			solved[n] = True
	print solved
	if solved[0] & solved[1] & solved[2]:
		print "All solved!"

	""" Create new fits files with the WCS corrections in the headers
	    new-wcs -w test_blue.wcs -i run010_b.fits -o blue_new.fits -v -d
	"""
	for n, c in enumerate(channels):
		if solved[n]:
			debug.write("Merging the WCS with the FITS file", level = 2)
			inputFitsFile = ultracamutils.addPaths(config.SITE_PATH, arg.runname) + '_' + c + "_n.fits"
			wcsFile = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + '_' + c + ".wcs"
			newFitsFile = ultracamutils.addPaths(config.SITE_PATH, arg.runname) + '_' + c + '_wcs.fits'
			mergeWCSCommand = ["new-wcs"]
			mergeWCSCommand.append("-w")
			mergeWCSCommand.append(wcsFile)
			mergeWCSCommand.append("-i")
			mergeWCSCommand.append(inputFitsFile)
			mergeWCSCommand.append("-o")
			mergeWCSCommand.append(newFitsFile)
			mergeWCSCommand.append("-d")
			subprocess.call(mergeWCSCommand)
			
	""" Write a .json object for the WCS solution to use in the web pages
	"""
	for n, c in enumerate(channels):
		if solved[n]:
			wcsFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + '_' + c + ".wcs"
			wcsFile = astropy.io.fits.open(wcsFilename)
			header = wcsFile[0].header
			wcs = wcsclasses.wcsSolution()
			
			equinox = float(header['EQUINOX'])
			referenceCoord = float(header['CRVAL1']), float(header['CRVAL2'])
			referencePixel = float(header['CRPIX1']), float(header['CRPIX2'])
			
			CD_array = [ [header['CD1_1'], header['CD1_2']], [ header['CD2_1'], header['CD2_2'] ] ]
			
			wcs.setSolution(equinox, referenceCoord, referencePixel, CD_array)
			
			print wcs
			
			print wcs.getWorldCoord( (512, 512) )
			#data = sexCatalog["LDAC_OBJECTS"].data
			#columns = sexCatalog["LDAC_OBJECTS"].columns
			#objects = []
	
			#for item in data:
			#	object = {}
			#	object['id']     = item[columns.names.index("NUMBER")]
			#	objects.append(object)
	
			wcsFile.close()
	
	
