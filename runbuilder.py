#!/usr/bin/env python
import ultracamutils
import sys, argparse, subprocess, re, os
import classes

if __name__ == "__main__":
	
		
	parser = argparse.ArgumentParser(description='Chains together the necessary steps to run the pipeline for one particular run')
	parser.add_argument('runname', type=str, help='Ultracam run name  [eg 2013-07-21/run010]')
	parser.add_argument('-d', '--debuglevel', default = 2, type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	parser.add_argument('-n', '--numframes', type=int, help='Number of frames (default = all frames)')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf')
	arg = parser.parse_args()

	config = ultracamutils.readConfigFile(arg.configfile)
	debug = classes.debugObject(config.DEBUG)
	debug.toggleTimeLog()
	if (arg.debuglevel!=None): debug.setLevel(arg.debuglevel);

	if arg.numframes!=None:
		numFrames = arg.numframes

	""" First check if a directory is made for the output files in the ucamsite folder. and create one. 
	"""
	runDate, runNumber = ultracamutils.separateRunNameAndDate(arg.runname)
	
	outputDirectory = ultracamutils.addPaths(config.SITE_PATH, runDate)
	
	if not os.path.exists(outputDirectory):
		debug.write("Creating the directory: " + outputDirectory)
		os.mkdir(outputDirectory)

	print arg
	dString = "-d" + str(arg.debuglevel)
	if arg.numframes!=None:
		nString = "-n" + str(arg.numframes)
		subprocess.call(["objectdbcreator.py", arg.runname, nString, dString])
	else:
		subprocess.call(["objectdbcreator.py", arg.runname, dString])

	xylsString = "--xyls"
	subprocess.call(["postprocessor.py", arg.runname, dString, xylsString])
	
	subprocess.call(["wcssolver.py", arg.runname, dString])

	subprocess.call(["mergecolours.py", arg.runname, dString])

	subprocess.call(["create_html.py", arg.runname])
	
	outputURL = ultracamutils.addPaths(config.ROOTURL, arg.runname + ".html")
	print "The output for this run is available at: %s"%outputURL
