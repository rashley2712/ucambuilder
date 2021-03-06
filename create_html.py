#!/usr/bin/env python
import jinja2
import ultracamutils
import classes
import sys, argparse

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Reads some meta data and builds the .html file for a particular run')
	parser.add_argument('runname', type=str, help='Ultracam run name  [eg 2013-07-21/run010]')
	parser.add_argument('-d', '--debuglevel', default = 1, type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf')
	parser.add_argument('-v', '--version', default='primary', help="Optional version string.")
	arg = parser.parse_args()

	config = ultracamutils.readConfigFile(arg.configfile)
	debug = classes.debugObject(config.DEBUG)
	debug.toggleTimeLog()
	if (arg.debuglevel!=None): debug.setLevel(arg.debuglevel);

	templateLoader = jinja2.FileSystemLoader(searchpath="/")
	templateEnv = jinja2.Environment( loader=templateLoader )
	
	date, runID = ultracamutils.separateRunNameAndDate(arg.runname)
	
	""" Get info about the run 
	"""
	runDate, runID = ultracamutils.separateRunNameAndDate(arg.runname)
	runInfo = classes.runObject(runDate, runID)
	runInfo.version = arg.version
	runInfo.loadSelf(config)
	if runInfo.comment=="":
		runInfo.checkForComments(config.ULTRACAMRAW)

	
	debug.write("Run Info:\n----------------------", level = 2)
	debug.write(runInfo, level = 2)
	debug.write("----------------------", level = 2)
	
	runDuration = ultracamutils.writeFriendlyTimeMinutes(runInfo.expose)
	debug.write("Approx run time: %s"%(runDuration), level =2)
	
	# Read the template file using the environment object.
	# This also constructs our Template object.
	template = templateEnv.get_template(config.RUNTEMPLATE)
	
	# Specify any input variables to the template as a dictionary.
	templateVars = { 	"runName" : runID, 
						"date" : date, 
						"comments" : runInfo.comment,
						"object" : runInfo.target,
						"duration" : runDuration,
						"version" : runInfo.version }
	
	outputFilename = ultracamutils.addPaths(config.SITE_PATH, arg.runname)
	if arg.version!='primary':
		outputFilename+="_" + str(arg.version)
	outputFilename+= ".html"
	
	debug.write("Writing HTML file to: %s"%(outputFilename))
	
	outFile = open(outputFilename, "w")
	outFile.write(template.render(templateVars))
	outFile.close()
	
	
	
