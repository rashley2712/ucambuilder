import numpy, json
from scipy.ndimage.filters import gaussian_filter
import ultracamutils

class aperture:
	""" This is an aperture instance 
	"""
	def __init__(self, id, position, radius):
		self.position = position
		self.radius = radius
		self.isAnnular = False
		self.id = id
		
	def pixelArea(self):
		return numpy.pi * self.radius * self.radius
	
	def fluxPerPixel(self, flux):
		return flux / self.pixelArea()	
	
class window:
	""" This class defines a window for the ULTRASPEC camera
	"""
	def __init__(self):
		self.xbin = 0
		self.ybin = 0
		self.xll = 0
		self.yll = 0
		self.nx = 0
		self.ny = 0
		self.data = None
		self.stackedData = []
		self.sources = None
		
	def setExtents(self, xll, yll, nx, ny):
		self.xll = xll
		self.yll = yll
		self.nx = nx
		self.ny = ny
		
	def setBinning(self, xbin, ybin):
		self.xbin = xbin
		self.ybin = ybin
		
	def setData(self, data):
		self.data = data
		if len(self.stackedData) == 0:
			self.stackedData = data

	def addData(self, data):
		self.stackedData = self.stackedData + data
		self.data = data
		
	def setSources(self, sources):
		self.sources = sources
		
	def getSources(self):
		return self.sources

class sourceMap:
	""" This class defines a histogram that is going to be used as a 'source map' or a 'heat map of sources'
	"""
	
	def __init__(self, dimensions):
		self.heatMap = numpy.zeros(dimensions)
		self.psize = 1.0
		self.fwhm = 4.0
		self.xsize = dimensions[1]
		self.ysize = dimensions[0]
		self.border = 10
		print "Intialised a source map with dimensions:", dimensions
		
	def updateMap(self, sources):
		for s in sources:
			j, i = int(s[0]), int(s[1])
			if (j<self.border) or (j>self.xsize-self.border): continue
			if (i<self.border) or (i>self.ysize-self.border): continue
			self.heatMap[i][j]+=1 
			
	def getSourceMap(self):
		return self.heatMap
		
	def getSmoothMap(self):
		return gaussian_filter(self.heatMap, self.fwhm/self.psize/2.3548, mode='constant')
		
class runInfo:
	""" This class is used to store the meta-data for the run
	"""
	def __init__(self, runPath):
		self.runPath = runPath
		self.runDate, self.runID = ultracamutils.separateRunNameAndDate(runPath)
		
	def loadFromJSON(self, JSONFilename):
		JSONfile = open(JSONFilename, "r")
		allObjectsJSON = json.load(JSONfile)
		run = {}
		runNumberStr = self.runID[3:]
		runNumber = int(runNumberStr)
	
		for object in allObjectsJSON:
			date = object['night']
			num = object['num']
			if ((date == self.runDate) & (runNumber == num)):
				self.comment = object["comment"]
				self.ra = object['ra']
				self.dec = object['dec']
				self.objectID = object['id']
				self.target = object['target']
				self.num = object['num']
				self.expose = object['expose']
				print object
			

