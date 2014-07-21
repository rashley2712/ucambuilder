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
	parser.add_argument('-v', '--version', default='primary', help="Optional version string.")
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
	
	# Run OBJECTDBCREATOR
	objectdbcreatorCommand = ['objectdbcreator.py']
	objectdbcreatorCommand.append(arg.runname)
	dString = "-d" + str(arg.debuglevel)
	objectdbcreatorCommand.append(dString)
	if arg.numframes!=None:
		objectdbcreatorCommand.append("-n" + str(arg.numframes))
	if arg.version!='primary':
		objectdbcreatorCommand.append("-v" + str(arg.version))
	subprocess.call(objectdbcreatorCommand)


	# Run POSTPROCESSOR
	postprocessorCommand = ['postprocessor.py']
	postprocessorCommand.append(arg.runname)
	postprocessorCommand.append("--xyls")
	postprocessorCommand.append(dString)
	if arg.version!='primary':
		postprocessorCommand.append('-v' + str(arg.version))
	subprocess.call(postprocessorCommand)
	
	# Run WCSSOLVER
	wcssolverCommand = ['wcssolver.py']
	wcssolverCommand.append(arg.runname)
	wcssolverCommand.append(dString)
	if arg.version!='primary':
		wcssolverCommand.append('-v' + str(arg.version))
	subprocess.call(wcssolverCommand)
	

	subprocess.call(["mergeobjects.py", arg.runname, dString])

	subprocess.call(["create_html.py", arg.runname, dString])
	
	outputURL = ultracamutils.addPaths(config.ROOTURL, arg.runname + ".html")
	print "The output for this run is available at: %s"%outputURL
