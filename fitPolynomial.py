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
	
def polynomial_gradient(n, a, x):
	result = 0	
	for i in range(1, n):
		result+= i * a[i] * (x**float(i-1))
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
	
	# Start by fitting a quadratic through three points in the data range
	datax = numpy.arange(len(source.offsetX))
	datay = numpy.array(source.offsetX)
	if len(datax) < 3:
		print "Too few points for a fit. Exiting."
		sys.exit()
	xpoints = [datax[0], datax[int(len(datax)/2)], datax[-1]]
	ypoints = [datay[0], datay[int(len(datax)/2)], datay[-1]]
	print zip(xpoints, ypoints)
	xMat = [ [ 1, xpoints[0], xpoints[0]**2 ] , \
	         [ 1, xpoints[1], xpoints[1]**2 ] , \
	         [ 1, xpoints[2], xpoints[2]**2 ] ]
	yMat = [ ypoints[0], ypoints[1], ypoints[2] ]
	a_quad = numpy.linalg.solve(xMat, yMat)
	print a_quad


	order = 9
	n = order+1
	a = numpy.zeros(n)
	a[0] = a_quad[0]
	a[1] = a_quad[1]
	a[2] = a_quad[2]
	fixed_args = ( n ) 
	delta = 0.1
	bounds = []
	scales = []
	for order, parameter in enumerate(a):
		limit = parameter + (delta * 10**(-(2*order)))
		scale = 1*10**(-order*3)
		bound = (-limit, limit)
		print parameter, bound
		bounds.append(bound)
		scales.append(scale)
	
	print a, bounds, scales
	
	xfit = numpy.arange(len(source.offsetX))
	yfit = [polynomial(n, a, x) for x in xfit]
	ppgplot.pgsci(2)
	ppgplot.pgline(xfit, yfit)
	
	result = numpy.polyfit(datax, datay, order)
	print "Numpy result", result
	a = result[::-1]
	yfit = [polynomial(n, a, x) for x in xfit]
	print "ChiSquared of the fit:", chiSquared(a)
	ppgplot.pgsci(2)
	ppgplot.pgline(xfit, yfit)
	
	sys.exit()
	
	results = scipy.optimize.minimize(chiSquared, a, method = 'SLSQP', jac=False, bounds=bounds)
	# results = scipy.optimize.fmin(chiSquared, a)
	# results = scipy.optimize.fmin_tnc(chiSquared, a, approx_grad=True, scale=scales)
	print results
	a = results.x
	yfit = [polynomial(n, a, x) for x in xfit]
	
	
	ppgplot.pgsci(1)
	ppgplot.pgenv(0, len(source.offsetX), -maxOffset, maxOffset, 0, 0)
	ppgplot.pgpt(range(len(source.offsetX)), source.offsetX, 1)
	ppgplot.pgsci(2)
	ppgplot.pgline(xfit, yfit)
	
	ppgplot.pgclos()
			
		
		