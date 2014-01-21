#!/usr/bin/env python
import rashley_utils as utils
import sys, subprocess, re
import classes

config = utils.readConfigFile()

if (len(sys.argv) < 2):
	print "Please give me a run name."
	sys.exit()

requestedNumFrames = -1

for i in sys.argv[2:]:
	if i[:2]=="-n": 
		print "number of frames:",  i[2:] 
		framesStr = i[2:]
		requestedNumFrames = int(framesStr)

runname = sys.argv[1]

runFilename = utils.addPaths(config.ULTRACAMRAW, runname)

print runname, runFilename

channel = ['r','g','b']

for i in channel:
	channelParam = "-c" + str(i)
	subprocess.call(["objecttracker.py", runname, "-n" + str(requestedNumFrames), channelParam])

