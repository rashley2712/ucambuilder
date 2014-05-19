import math, json
import ultracamutils

class ucamObject:
	""" This class is the storage holder for a merged object
	"""
	
	def __init__(self, id):
		self.id = id
		photometry = { 'magnitude': 0, 'MJDIndex': 0 }
		self.photometryMeasurements = 0
		self.isComparison = False;
		

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
			
	
		
