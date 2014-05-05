#!/usr/bin/env python
import ultracamutils
import jinja2
import sys, argparse, subprocess, re, os, json
import classes

if __name__ == "__main__":
	
		
	parser = argparse.ArgumentParser(description='Chains together the necessary steps to run the pipeline for one full day')
	parser.add_argument('date', type=str, help='Ultracam date  [eg 2013-07-21]')
	parser.add_argument('-d', '--debuglevel', default = 2, type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	parser.add_argument('-n', '--numframes', type=int, help='Number of frames (default = all frames)')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf')
	parser.add_argument('-r', '--buildruns', action='store_true', help='Build the run output for each run (if no existing output found)')
	parser.add_argument('-f', '--forcebuildruns', action='store_true', help='Force build of each run (even if existing data is found)')
	arg = parser.parse_args()

	config = ultracamutils.readConfigFile(arg.configfile)
	debug = classes.debugObject(config.DEBUG)
	debug.toggleTimeLog()
	if (arg.debuglevel!=None): debug.setLevel(arg.debuglevel);

	if arg.numframes!=None:
		numFrames = arg.numframes

	""" First check if a directory for the output files exists in the ucamsite folder and if not, create one. 
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

	runList = []
	runs_re = re.compile("run[0-9][0-9][0-9].xml")
	for i, f in enumerate(fileList):
		m = runs_re.match(f)
		if (m): 
			runName = m.group()[:6]
			runList.append(runName)
			
	debug.write("This is the list of runs found in the folder for the date %s: %s"%(arg.date, str(runList)))

	runData = []

	for runID in runList:
		newRun = classes.runObject(arg.date, runID)
		newRun.loadSelf(config.SITE_PATH)
		if newRun.ra == 0:
			# Failed to get info from local file.... go to Tom's ultra.json file for some meta-data
			runname = ultracamutils.addPaths(arg.date, runID)
			additionalRunInfo = ultracamutils.readULTRAJSON(config.RUNINFO, runname)
			print "Falling back to:", runInfo
	
		#newRun.writeSelf(config.SITE_PATH)
		runData.append(newRun)
	
	for run in runData:
		runID = run.runID
		runDate = arg.date
		runname = ultracamutils.addPaths(runDate, runID)
		runInfo = ultracamutils.getRunInfo(config.RUNINFO, runname)
		runInfo.runDuration = ultracamutils.writeFriendlyTimeMinutes(runInfo.expose)

		# Check if the run has been processed (ie if the file runxxx_info.json exists)
		runMetaDataFilename = ultracamutils.addPaths(config.SITE_PATH, runname) + "_info.json"
		print "Checking for meta-data (JSON) file: ", runMetaDataFilename



		if (arg.buildruns):
			
			runbuilderCommand = ["runbuilder.py"]
			runbuilderCommand.append(run.runDate + "/" + run.runID)
			if (arg.numframes!=None):
				runbuilderCommand.append("-n")
				runbuilderCommand.append(arg.numframes)
			
			
			debug.write("Running runbuilder with:" + str(runbuilderCommand))
			subprocess.call(runbuilderCommand)

	
		# Check if the run has been processed (ie if the file runxxx_info.json exists)
		runMetaDataFilename = ultracamutils.addPaths(config.SITE_PATH, runname) + "_info.json"
		print "Checking for meta-data (JSON) file gat: ", runMetaDataFilename
		if os.path.exists(runMetaDataFilename):
			print "Run output exists!"
			JSONfile = open(runMetaDataFilename, "r")
			wholeFileString = JSONfile.read()
			metaDataObject = json.loads(wholeFileString)
			print metaDataObject
			maxExtents = metaDataObject['maxExtents']
			x1 = maxExtents[0]
			x2 = maxExtents[1]
			y1 = maxExtents[2]
			y2 = maxExtents[3]
			print (x1, y1, x2, y2)
			
			colours = ['r', 'g', 'b']
			for c in colours:
				imageFilename = ultracamutils.addPaths(config.SITE_PATH, runname) + '_' + c + '.png'
				outputFilename = ultracamutils.addPaths(config.SITE_PATH, runname) + '_' + c + '_thumb.png'

				thumbnailCommand = ["createthumbnail.py"]
				thumbnailCommand.append(imageFilename)
				thumbnailCommand.append("--autocrop")
				thumbnailCommand.append("-w")
				thumbnailCommand.append("128")
				thumbnailCommand.append("-o")
				thumbnailCommand.append(outputFilename)
				
				debug.write("Creating thumbnail" + str(thumbnailCommand))
				subprocess.call(thumbnailCommand)
	
			runInfo.thumbnailURI = runID + '_r_thumb.png'
			runInfo.runURL =  runID + ".html" 
		else:
			runInfo.thumbnailURI = "../default_thumbnail.png"
			runInfo.runURL = ""
			
		runData.append(runInfo)

	dayData.setRuns(runData)


	# Initialise the Jinja environment
	templateLoader = jinja2.FileSystemLoader(searchpath="/")
	templateEnv = jinja2.Environment( loader=templateLoader )

	# Read the template file using the environment object.
	# This also constructs our Template object.
	template = templateEnv.get_template(config.DAYTEMPLATE)
	
	
	# Specify any input variables to the template as a dictionary.
	templateVars = { 	"date" : arg.date,
						"runs" : runData
						 }
	
	
	outputPath = ultracamutils.addPaths(config.SITE_PATH, arg.date)
	outputFilename = "index.html"
	
	fullFilename = ultracamutils.addPaths(outputPath, outputFilename)
	outFile = open(fullFilename, "w")
	outFile.write(template.render(templateVars))
	outFile.close()

	outputURL = config.ROOTURL + arg.date + "/" + outputFilename
	print "Browse the page at: ", outputURL
	

	
