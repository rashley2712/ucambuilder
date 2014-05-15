import math, json
import ultracamutils


class ucamObject:
	""" This class is the storage holder for a merged object
	"""
	
	def __init__(self, id):
		self.id = id
		photometry = { 'magnitude': 0, 'MJDIndex': 0, 
		self.photometryMeasurements = {'r': []}
		self.isComparison = False;
		

class timeLine:
	""" This class encapsulates info about the individual frames from a run
	"""
	
	def __init(self):
		self.MJDs = []            # Keeps track of the MJD for each frame
		self.objectsPerFrame[]    # Keeps track of the count of identified objects for each frame
		self.numFrames = 0
		self.duration = 0
