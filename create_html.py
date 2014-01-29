#!/usr/bin/env python
import jinja2
import rashley_utils as utils
import classes
import sys

def getComments(date, runName):
	""" Loads the PO comments from the log file in the directory
	"""
	comments = ""
	commentsFilename = utils.addPaths(config.ULTRACAMRAW, date)
	commentsFilename+= "/" + date + ".dat"
	commentsFile = open(commentsFilename, "r")
	
	for line in commentsFile:
		run = line.split(' ')[0]
		if (run==runName): comments = line[7:]

	commentsFile.close()
	
	return comments

templateLoader = jinja2.FileSystemLoader(searchpath="/")
templateEnv = jinja2.Environment( loader=templateLoader )

config = utils.readConfigFile()

if (len(sys.argv) < 2):
	print "Please give me a run name."
	sys.exit()

runname = sys.argv[1]

rundate = utils.separateRunNameAndDate(runname)
date = rundate[0]
runName = rundate[1]
comments = getComments(date, runName)

# Read the template file using the environment object.
# This also constructs our Template object.
template = templateEnv.get_template(config.RUNTEMPLATE)

# Specify any input variables to the template as a dictionary.
templateVars = { 	"runName" : runName, 
					"date" : date, 
					"comments" : comments }

outputFilename = utils.addPaths(config.SITE_PATH, date)
outputFilename = utils.addPaths(outputFilename, runName)
outputFilename+= ".html"

outFile = open(outputFilename, "w")
outFile.write(template.render(templateVars))
outFile.close()


