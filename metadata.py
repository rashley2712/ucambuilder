#!/usr/bin/env python
import rashley_utils as utils
import sys, subprocess, re, random
import os
import classes

class runDate: 
	def __init__(self, date):
		self.date = date
		self.runs = []
		
	def addRun(self, runName, comment=""):
		run = classes.runObject(self.date, runName)
		run.setComment(comment)
		self.runs.append(run)
		
	def getRuns(self):
		return self.runs
		
	def __str__(self):
		return str(self.date)
		

def loadAllComments(date):
	""" Returns dictionary of run names and comments for a particular date
	"""	
	comments = {}

	commentsFilename = utils.addPaths(config.ULTRACAMRAW, date)
	commentsFilename+= "/" + date + ".dat"
	commentsFile = open(commentsFilename, "r")
	
	runname_re = re.compile(r"run[0-9][0-9][0-9]")

	for line in commentsFile:
		run = line.split(' ')[0]
		r = runname_re.search(run)
		if r:
			comment = line[7:]
			comments[run] = comment

	commentsFile.close()

	return comments
				

if __name__ == "__main__":

	config = utils.readConfigFile()
	debug = classes.debugObject(config.DEBUG)

	dates = []
	rawFolders = os.listdir(config.ULTRACAMRAW) 

	date_re = re.compile(r'20[0-9]{2}(-[0-9]{2}){2}')
	for dateFolder in rawFolders:
		d = date_re.search(dateFolder)
		if (d):
			date = d.group(0)
			debug.write("Found a date: " + date)
			dateObject = runDate(date)
			dates.append(dateObject)
			

	# For debugging purposes, we run over a sample set of the dates, just 5 or so (sample_size)
	sample_size = 5
	sample_start = int(random.random()*len(dates)-sample_size)
	
	sample_size = len(dates)
	sample_start = 0

	totalRuns = 0
	totalFileSize = 0
	runs_re = re.compile(r'run[0-9][0-9][0-9].xml')
	for index in range(sample_start, sample_start+sample_size):
		d = dates[index]
		debug.write("Processing: " + d.date, level = 1)
		files = os.listdir(utils.addPaths(config.ULTRACAMRAW,d.date))
		commentsList = loadAllComments(d.date)
		for f in files:
			debug.write("file: " + f)
			r = runs_re.search(f)
			if (r):
				runname = r.group(0)[:6]
				# Get the filesize (.dat) of this run
				runPath = utils.addPaths(config.ULTRACAMRAW,d.date)
				runPath = utils.addPaths(runPath, runname)
				runPath+= ".dat"
				try: 
					totalFileSize+= os.path.getsize(runPath)
				except:
					pass
				# Get the observer's comments (if they exist)
				try:
					comment = commentsList[runname]
				except KeyError:
					debug.write("No comment found, making it blank")
					comment = ""
				# Get some more meta-data about the run, by creating a 'runObject' which will try to use trm.ultracam to read the start of the file and extract info about it
				try: 
					d.addRun(runname, comment = comment)
					totalRuns+=1
				except:
					debug.write("Couldn't add the run... " + d.date + "/" + runname, level=1)

	totalFrames = 0
	for index in range(sample_start, sample_start+sample_size):
		d = dates[index]
		for r in d.getRuns():
			totalFrames+=r.numFrames*2 + r.numFrames/r.nblue
				
	print "Total number of 'nights': ", len(dates)
	print "Total number of runs: ", totalRuns  
	print "Total number of frames: ", totalFrames
	totalTeraBytes = totalFileSize / 1000. / 1000. / 1000. / 1000.  # This is the SI decimal definition of a terabyte 
	print "Total file size (.dat files only): %d (bytes) %f (tera-bytes)"%(totalFileSize,totalTeraBytes)

