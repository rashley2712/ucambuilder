#!/usr/bin/env python
import ultracamutils
import sys, argparse, subprocess, re, os
import classes

if __name__ == "__main__":
	
		
	parser = argparse.ArgumentParser(description='Chains together the necessary steps to run the pipeline for one full day')
	parser.add_argument('date', type=str, help='Ultracam date  [eg 2013-07-21]')
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

	""" First check if a directory is made for the output files in the ucamsite folder and create one. 
	"""
	
	outputDirectory = ultracamutils.addPaths(config.SITE_PATH, arg.date)
	
	if not os.path.exists(outputDirectory):
		debug.write("Creating the directory: " + outputDirectory)
		os.mkdir(outputDirectory)

	debug.write("Arguments: " + str(arg))

	""" Get a list of all of the .dat files in the folder. This is going to be our list of 'runs'
	"""
	
	path = ultracamutils.addPaths(config.ULTRACAMRAW, arg.date)

	filer = subprocess.Popen(["ls", str(path)], stdout = subprocess.PIPE)

	fileList = filer.communicate()[0].split('\n')

	dayData = classes.dayObject(arg.date)

	runs_re = re.compile("run[0-9][0-9][0-9].xml")
	for i, f in enumerate(fileList):
		m = runs_re.match(f)
		if (m): 
			runName = m.group()[:6]
			dayData.addRun(runName)

	for run in dayData.getRuns():
		print run
