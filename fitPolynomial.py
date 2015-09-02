#!/usr/bin/env python
import sys
import argparse
import ultracamutils
import ppgplot
import csv
import numpy
import scipy.optimize


def readFromCSV(filename):
	inputFile = open(filename, 'r')
	csvReader = csv.reader(inputFile)
	firstLine = csvReader.next()
	attributes = []
	for a in firstLine:
		attributes.append(a.strip())
	
	objects = []
	for line in csvReader:
		record = {}
		for i, l in enumerate(line):
			value = float(l.strip())
			record[attributes[i]] = value
		objects.append(record)
		
	return objects
	
class referenceSource():
	def init(self):
		self.xPositions = []
		self.yPositions = []
		self.startX = 0
		self.startY = 0
		self.offsetX = []
		self.offsetY = []
		
	def setStartPosition(self):
		self.startX = self.xPositions[0]
		self.startY = self.yPositions[0]
		self.calcOffsets()
		
	def calcOffsets(self):
		offsetX = [x - self.startX for x in self.xPositions]
		offsetY = [y - self.startY for y in self.yPositions]
		self.offsetX = offsetX
		self.offsetY = offsetY
	
def chiSquared(a, *args):
	global datax, datay, n
	chiSq = 0
	for x, y in zip(datax, datay):
		chiSq+= (polynomial(n, a, x) - y)**2
		
	yfit = [polynomial(n, a, x) for x in xfit]
	ppgplot.pgsci(3)
	ppgplot.pgline(xfit, yfit)
	print chiSq, a
	return chiSq	
	
def polynomial(n, a, x):
	result = 0	
	for i in range(n):
		result+= a[i] * (x**float(i))
	return result

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Loads the positions of a reference aperture and fits a polynomial to the data.')
	parser.add_argument('runname', type=str, help='Ultraspec run name  [eg 2013-07-21/run010]')
	parser.add_argument('-d', '--debuglevel', type=int, help='Debug level: 3 - verbose, 2 - normal, 1 - warnings only')
	parser.add_argument('-c', '--configfile', default='ucambuilder.conf', help='The config file, usually ucambuilder.conf')
	
	
	arg = parser.parse_args()
	
	config = ultracamutils.readConfigFile(arg.configfile)
	
	inputFilename = ultracamutils.addPaths(config.WORKINGDIR, arg.runname) + "_reference_aperture.csv"
	
	positions = readFromCSV(inputFilename)
	
	pixelPositionView = ppgplot.pgopen('/xs')
	
	source = referenceSource()
	source.xPositions = [p['X_ABS'] for p in positions]
	source.yPositions = [p['Y_ABS'] for p in positions]
	source.setStartPosition()
	
	maxOffset = numpy.max(source.offsetX)
	minOffset = numpy.min(source.offsetX)
	if abs(minOffset) > maxOffset:
		maxOffset = - minOffset
	maxOffset*= 1.1
	ppgplot.pgenv(0, len(source.offsetX), -maxOffset, maxOffset, 0, 0)
	ppgplot.pgask(False)
	
	ppgplot.pgpt(range(len(source.offsetX)), source.offsetX, 1)
	
	n = 6
	a = numpy.zeros(n) 
	
	xfit = numpy.arange(len(source.offsetX))
	yfit = [polynomial(n, a, x) for x in xfit]
	ppgplot.pgsci(2)
	ppgplot.pgline(xfit, yfit)
	datax = xfit
	datay = source.offsetX
	results = scipy.optimize.minimize(chiSquared, a, args = n, method = 'Nelder-Mead')
	print results
	a = results.x
	yfit = [polynomial(n, a, x) for x in xfit]
	
	ppgplot.pgsci(1)
	ppgplot.pgenv(0, len(source.offsetX), -maxOffset, maxOffset, 0, 0)
	ppgplot.pgpt(range(len(source.offsetX)), source.offsetX, 1)
	ppgplot.pgsci(2)
	ppgplot.pgline(xfit, yfit)
	
	ppgplot.pgclos()
			
		
		