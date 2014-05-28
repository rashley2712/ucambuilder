import math, json, numpy
import ultracamutils

class colourObject:
	""" This class is the storage holder for a merged object
	"""
	
	colours = ['r', 'g', 'b']
	
	def __init__(self, id):
		self.id = id
		self.isComparison = False;
		self.meanPosition = {'r': (0, 0), 'g': (0, 0), 'b': (0, 0)}
		self.photometry = {'r': [], 'g': [], 'b':[] }
		self.colourID = { 'r': -1, 'g': -1, 'b': -1 }
		self.deviations = { 'r': -1, 'g': -1, 'b': -1 }
		self.comparisonFlags = { 'r': False, 'g': False, 'b': False }
		
	def setMeanPosition(self, colour, meanPosition):
		self.meanPosition[colour] = meanPosition
		
	def addExposure(self, colour, exposure):
		self.photometry[colour].append(exposure)
		
	def calculateSigma(self):
		for c in colourObject.colours:
			data = self.photometry[c]
			measurement = []
			for e in data:
				measurement.append(e['magnitude'])
			standardDeviation = numpy.std(measurement)
			mean = numpy.mean(measurement)
			fractionalstd = standardDeviation/mean
			print "%s mean: %f, stdev: %f, std/mean: %f"%(c, mean, standardDeviation, fractionalstd)
			self.deviations[c] = fractionalstd
			
	def getPhotometry(self, c, frameIndex):
		data = self.photometry[c]
		for e in data:
			if e['frameIndex'] == frameIndex:
				return e['magnitude']
			
		return -1
		

	def __str__(self):
		retStr = "ID: %d\n"%self.id
		for c in colourObject.colours:
			retStr+= '%s(%d,%d)[%d][%d] '%(c, int(self.meanPosition[c][0]), int(self.meanPosition[c][1]), len(self.photometry[c]), self.colourID[c])
			
		return retStr
		
	def setComparisonFlag(self, colour):
		self.comparisonFlags[colour] = True
		
	def testComparison(self):
		if (self.comparisonFlags['r'] & self.comparisonFlags['g'] & self.comparisonFlags['b']): 
			self.isComparison = True
			return True
		else: 
			return False
			
			 
		
	def toJSON(self):
		testObject = {}
		testObject['id'] = self.id
		testObject['isComparison'] = self.isComparison
		testObject['colourID'] = self.colourID
		testObject['meanPosition'] = self.meanPosition
		shortPhotometry = {'r': [], 'g': [], 'b':[]}
		for c in colourObject.colours:
			measurements = self.photometry[c]
			for m in measurements:
				reading = [ m['frameIndex'], m['magnitude'], m['fwhm'], m['position'][0], m['position'][1] ]
				shortPhotometry[c].append(reading)
		testObject['photometry'] = shortPhotometry
		jsonString = json.dumps(testObject)
		return jsonString
		

class timeLine:
	""" This class encapsulates info about the individual frames from a run
	"""
	
	def __init(self):
		self.MJDs = []            		# Keeps track of the MJD for each frame
		self.objectsPerFrame = []    	# Keeps track of the count of identified objects for each frame
		self.numFrames = 0
		self.duration = 0

class frameObject:
	""" Encapsulates info about a particular frame
	"""
	
	def __init__(self, *args):
		if len(args)>0: self.MJD = args[0]
		else: self.MJD = 0
		if len(args)>1: self.frameNumber = args[1]
		else: self.frameNumber = 0
		if len(args)>2: self.frameIndex = args[2]
		else: self.frameIndex = 0
		self.objects = {'r':0, 'g':0, 'b':0}
		
	def setObjectCount(self, colour, number):
		self.objects[colour] = number
		
	def __str__(self):
		outStr = "Index: %d MJD: %s "%(self.frameIndex, str(self.MJD))
		outStr+= "[r: %d, g: %d, b:%d]"%(self.objects['r'], self.objects['g'], self.objects['b'])
		outStr+= " frameNumber: %d"%self.frameNumber
		return outStr
	
	def toJSON(self):
		testObject = {'frameNumber': self.frameNumber, 'frameIndex': self.frameIndex}
		testObject['objects'] = self.objects
		testObject['MJD'] = self.MJD
		
		jsonString = json.dumps(testObject)
		
		return jsonString
		
	def setFromObject(self, object):
		self.MJD = object['MJD']
		self.frameNumber = object['frameNumber']
		self.frameIndex = object['frameIndex']
		self.objects['r'] = object['objects']['r']
		self.objects['g'] = object['objects']['g']
		self.objects['b'] = object['objects']['b']
			
	
		
