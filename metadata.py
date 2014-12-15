#!/usr/bin/env python
import rashley_utils as utils
import sys, subprocess, re, random
import os, json
import classes
import numpy, matplotlib.pyplot, scipy.stats

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
		

def loadULTRAJSON(filename):
	""" Loads the ultra.json file and returns it as an array.
	"""
	JSONfile = open(filename, "r")
	allObjectsJSON = json.load(JSONfile)
	return allObjectsJSON

def getULTRAmatch(runid, rundate, ultraObjects):
	runNumberStr = runid[3:]
	runNumber = int(runNumberStr)
	for o in ultraObjects:
		date = o['night']
		num = o['num']
		if (date==rundate) & (num == runNumber):
			return o

	print "No match for:", runid, rundate
	return None
			
			


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
	sample_size = 40
	sample_start = int(random.random()*len(dates)-sample_size)
	
	sample_size = len(dates)
	sample_start = 0

	totalRuns = 0
	totalFileSize = 0
	runs_re = re.compile(r'run[0-9][0-9][0-9].xml')
	for index in range(sample_start, sample_start+sample_size):
		d = dates[index]
		debug.write("%d/%d Processing: %s"%(index-sample_start+1, sample_size, d.date), level = 1)
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

	print "Loading ultra.json"
	ultraobjects = loadULTRAJSON(config.RUNINFO)
	totalFrames = 0
	frameTimes = []
	identifiers = []
	comments = []
	runLengths = []
	for index in range(sample_start, sample_start+sample_size):
		d = dates[index]
		for r in d.getRuns():
			totalFrames+=r.numFrames*2 + r.numFrames/r.nblue
			ultra = getULTRAmatch(r.runID, r.runDate, ultraobjects)
			print r
			if (ultra!=None):
				if ultra['expose']!=0:
					r.frameTime = float(ultra['expose'])/r.numFrames * 60
					print "Frametime:",r.frameTime
					frameTimes.append(r.frameTime)
					identifiers.append(str(ultra['night'] + ' / ' + str(ultra['num'])))
					comments.append( ultra['target'] + " : " + r.comment)
					runLengths.append(ultra['expose'])
					
	print "Total number of 'nights': ", len(dates)
	print "Total number of runs: ", totalRuns  
	print "Total number of frames: ", totalFrames
	totalTeraBytes = totalFileSize / 1000. / 1000. / 1000. / 1000.  # This is the SI decimal definition of a terabyte 
	print "Total file size (.dat files only): %d (bytes) %f (tera-bytes)"%(totalFileSize,totalTeraBytes)

	#print frameTimes
	print "run lengths:", runLengths
	
	trimAt = 10
	trimmedRunLengths = []
	for r in runLengths:
		if r>trimAt: trimmedRunLengths.append(r)
	
	binwidth = 10
	bottomofbins = trimAt
	topofbins = 600
	numbins = 1 + ( (topofbins - bottomofbins) / binwidth)
	bins = numpy.linspace(bottomofbins, topofbins, numbins)
	
	hist, bins, patches = matplotlib.pyplot.hist(trimmedRunLengths, bins, histtype='bar')
	
	fig = matplotlib.pyplot.gcf()
	
	
	matplotlib.pyplot.xlabel('Run length (minutes)')
	matplotlib.pyplot.ylabel('Number of runs')
	ax = matplotlib.pyplot.gca()
	#ax.set_yscale('log')
	matplotlib.pyplot.yscale('log', nonposy='clip')
	matplotlib.pyplot.xlim([trimAt, 600])
	#matplotlib.pyplot.axis([bottomofbins, topofbins, min(n), max(n)])
	
	print hist
	frameTimes = numpy.array(frameTimes)
	print "Number of runs used for this analysis:", len(frameTimes)
	print "Longest exposure time", max(frameTimes), identifiers[frameTimes.argmax()]
	print "Shortest exposure time", min(frameTimes), identifiers[frameTimes.argmin()]
	
	values = zip(hist, bins)
	print values
	
	matplotlib.pyplot.show()
	
	DPI = fig.get_dpi()
	print "DPI:", DPI
	DefaultSize = fig.get_size_inches()
	print "Default size in Inches", DefaultSize
	print "Which should result in a %i x %i Image"%(DPI*DefaultSize[0], DPI*DefaultSize[1])

	
	fig.savefig('hist.eps',dpi=100, format='eps')

	print "20 second examples"
	for index in range(len(frameTimes)):
		t = frameTimes[index]
		if (t>19) & (t<22):
			print "exposuretime:", t 
			print identifiers[index]
			print comments[index]
	
