import numpy

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
		self.stackedData = None
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
		if self.stackedData == None:
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
		print "Intialised a source map with dimensions:", dimensions
		
	def updateMap(self, sources):
		for s in sources:
			j, i = int(s[0]), int(s[1])
			self.heatMap[i][j]+=1
			
	def getSourceMap(self):
		return self.heatMap