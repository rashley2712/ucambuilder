import numpy, json
import xml.etree.ElementTree as ElementTree
from scipy.ndimage.filters import gaussian_filter
import ultracamutils
from astropy.table import Table

class aperture:
	""" This is an aperture instance 
	"""
	def __init__(self, id, position, radius):
		self.position = position
		self.radius = radius
		self.hasSky = False
		self.skyInnerRadius = 0
		self.skyOuterRadius = 0
		self.id = id
		
	def pixelArea(self):
		return numpy.pi * self.radius * self.radius
	
	def fluxPerPixel(self, flux):
		return flux / self.pixelArea()	
		
	def setSkyAperture(self, skyInnerRadius, skyOuterRadius):
		self.hasSky = True
		self.skyInnerRadius = skyInnerRadius
		self.skyOuterRadius = skyOuterRadius
		
	def setRadius(self, radius):
		self.radius = radius
		
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
		self.BGSubtractedImage = []
		self.sources = None
		self.borderWidth = 10
		
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
	
	def setBlankData(self, data):
		dataShape = numpy.shape(data)
		blanks = numpy.zeros(dataShape)
		if len(self.stackedData) == 0:
			self.stackedData = blanks

	def addData(self, data):
		self.stackedData = self.stackedData + data
		self.data = data
		
	def addToStack(self, data):
		self.stackedData = self.stackedData + data
		
	def filterBorderSources(self):
		currentNumSources = len(self.sources)
		if currentNumSources==0:
			return 0
		sourcesToRemove = []
		for index, s in enumerate(self.sources):
			x = s['xcentroid']
			y = s['ycentroid']
			if (x < self.borderWidth) or (x > (self.nx - self.borderWidth)): 
				sourcesToRemove.append(index)
				continue
			if (y < self.borderWidth) or (y > (self.ny - self.borderWidth)): 
				sourcesToRemove.append(index)
				continue
		self.sources.remove_rows(sourcesToRemove)
		return len(self.sources)
		
	def setSources(self, sources):
		self.sources = sources
		
	def setSourcesAvoidBorders(self, sources):
		self.sources = sources
		self.filterBorderSources()
		
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
		self.comment = ""
		self.ra = 0
		self.dec = 0
		
	def loadFromJSON(self, JSONFilename):
		JSONfile = open(JSONFilename, "r")
		allObjectsJSON = json.load(JSONfile)
		run = {}
		runNumberStr = self.runID[3:]
		runNumber = int(runNumberStr)
		
		objectNotFound = True
		
		for object in allObjectsJSON:
			date = object['night']
			num = object['num']
			if ((date == self.runDate) & (runNumber == num)):
				objectNotFound = False
				self.comment = object["comment"]
				self.ra = object['ra'] * 15. # Convert the ultra.json RA value to degrees
				self.dec = object['dec']
				self.objectID = object['id']
				self.target = object['target']
				self.num = object['num']
				self.expose = object['expose']
				print object
				return True 
		return False
            

	def loadFromXML(self, rawDataDirectory):
		""" This method is called upon if the loadFromJSON method fails to return any data
		"""
		print "Trying to load run info from the XML file"
		XMLFilename = ultracamutils.addPaths(rawDataDirectory, self.runPath) + ".xml"
		print XMLFilename
		try:
			tree = ElementTree.parse(XMLFilename)
		except IOError as e:
			print "Could not get read the run's XML file."
			return False
		
		root = tree.getroot()
		user = root.find('user')
		target = user.find('target').text
		PI = user.find('PI').text
		ID = user.find('ID').text
		observers = user.find('Observers').text
		raStr = user.find('RA').text
		decStr = user.find('Dec').text
		ra, dec = ultracamutils.fromSexagesimal(raStr, decStr)
		print "Target:", target
		print "PI: %s Observers: %s  Programme: %s"%(PI, observers, ID)
		print "(ra, dec): (%f, %f)"%(ra, dec)
		
		self.ra = ra
		self.dec = dec
		self.target = target
		
		return True
        