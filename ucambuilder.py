#!/usr/bin/env python
import rashley_utils as utils
import sys, subprocess, re
import classes

config = utils.readConfigFile()

if (len(sys.argv) < 2):
	print "Please give me a date."
	sys.exit()


date = sys.argv[1]

path = utils.addPaths(config.ULTRACAMRAW, date)

filer = subprocess.Popen(["ls", str(path)], stdout = subprocess.PIPE)

fileList = filer.communicate()[0].split('\n')

dayData = classes.dayObject(date)

runs_re = re.compile("run[0-9][0-9][0-9].xml")
for i, f in enumerate(fileList):
	m = runs_re.match(f)
	if (m): 
	      runName = m.group()[:6]
	      dayData.addRun(runName)

newFolder = utils.addPaths(config.SITE_PATH, date)
subprocess.call(["mkdir", newFolder])

for run in dayData.runs: 
	runFilename = utils.addPaths(date, str(run.runName))
	print runFilename
<<<<<<< HEAD
	subprocess.call(["threecolours.py", runFilename, "-n100"])
	subprocess.call(["postprocessor.py", runFilename])
	subprocess.call(["create_html.py", runFilename])
=======
	subprocess.call(["objecttracker.py", runFilename, "-n100"])
>>>>>>> f97547520688a895e0b7d385fea8ba94409d03e1
