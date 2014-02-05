import math, json, datetime
import trm.ultracam
import rashley_utils as utils

class dayObject:
    """ This class stores metadata for a particular night's observing
	"""
    def __init__(self, date):
        self.date = date
        self.runs = []
        self.totalTime = 0
        self.totalFrames = 0
        
    def addRun(self, runName):
        run = runObject(self.date, runName)
        self.runs.append(run)

class runObject:
	def __init__(self, date, runName):
		self.runName = runName
		config = utils.readConfigFile()
		runPath = utils.addPaths(config.ULTRACAMRAW, date)
		runPath = utils.addPaths(runPath, runName)
		self.totalTime = 0
		runMetaData = trm.ultracam.Rhead(runPath, server=False)
		self.mode = runMetaData.mode
		self.userData = runMetaData.user
		try: 
			self.nblue = runMetaData.nblue
		except AttributeError:
			self.nblue = 1
		if (self.mode != "PONOFF"):
			runData =  trm.ultracam.Rdata(runPath, 1, server=False)
			self.numFrames = runData.ntotal()
		else: 
			self.numFrames = 0
		self.runClass = 0
		self.comment = ""
		
		
	def __str__(self):
		outStr = "RunName: " + self.runName + "\n"
		outStr+= "Frames: " + str(self.numFrames) + "\n"
		outStr+= "Mode: " + self.mode + "\n"
		outStr+= "nBlue: " + str(self.nblue) + "\n"
		outStr+= "Comments: " + str(self.comment) + "\n"
		return outStr
		
	def determineRunClass(self):
		""" Determines the type of run returns as an integer and label  (ntuple)
		0. Unknown
		1. Full field science run (a)... Many objects (>40), large windows (~ full fields), 100-2000 frames
		2. Intermediate science run (b)... Intermediate objects (4-40), intermediate windows (~ half fields), 100-2000 frames
		3. High speed science run (c)...  Few objects (2-4), short exposures, DRIFT mode, small windows (~ quarter size), many frames (2000-100000)
		4. Acquisition run. Science run (1, 2, 3) with few frames <50 and several glitches
		5. Bias run 
		6. Flat field. Full CCD exposed. Start or end of the evening. ~100 exposures
		7. Junk (unclassified)
		"""
		runClassLabels = ["Unknown", "Science A", "Science B", "Science C", "Acquisition", "Bias", "Flat", "Junk"]
		
		if (self.numFrames>2000):
			self.runClass = 1
			return 1, runClassLabels[1]
			
		self.runClass = 0	
		return 0, runClassLabels[0] 

	def setComment(self, comment):
		self.comment = comment


class configObject:
    """ This class stores all of the configuration for the run
	"""
    def __init__(self):
        self.ULTRACAMRAW = "/storage/astro1/phsaap/ultracam/raw_data"
        self.DEBUG = 1
        self.SITE_PATH = "/storage/astro2/phrnaw/ucamsite"
        self.TMP_PATH = "/tmp"
        self.WRITE_FITS = 0
        self.KEEP_TMP_FILES = 0
        self.WRITE_JSON = 0
        self.MOVIE_TMP_PATH = "/tmp/movie"
        self.MINPIXELDISTANCE = 5
        self.FONT = "/usr/share/fonts/truetype/ubuntu-font-family/Ubuntu-B.ttf"
        self.RUNTEMPLATE = "/home/rashley/astro/ucamsite/templates/runxxx.jinja"


    def __setitem__(self, item, value):
		setattr(self, item, value)
	
    def __str__(self):
		out = ""
		for key, value in self.__dict__.items():
			out+= str(key) + " --> " + str(value) + "\n"
		return out

class ExposureObject:
	""" This is a class to encapsulate a single exposure
	"""
	def __init__(self):
		self.centroid = (0,0)     	# tuple containing the position of the object in this exposure
		self.counts = 0             	# counts measured by sextractor for this exposure
		self.exposureTime = 0       	# the exposure time for this exposure as given by Tom's ucam library
		self.MJD = 0                	# the calculated MJD for this exposure from Tom's ucam library
		self.FWHM = 0 			# the width of the image, as calculated by sextractor  
		
	def __str__(self):
		return "MJD: " + str(self.MJD) + " counts: " + str(self.counts)
	   
class FrameObject:
	""" This class contains the definitions of the windows for a specific run
	"""
	def __init__(self):
		self.channel = 'unknown'
		self.numWindows = 0
		self.windows = []
		self.nxmax = 0
		self.nymax = 0
		
	def addWindow(self, xoffset, yoffset, xsize, ysize):
		window = WindowObject(xoffset, yoffset, xsize, ysize)
		self.windows.append(window)
		self.numWindows+= 1
		
	def getWindow(self, index):
		return self.windows[index]
		
	def __str__(self):
		out = "Channel: " + str(self.channel) + " numWindows: " + str(self.numWindows) + " dimensions (" + str(self.nxmax) + ", " + str(self.nymax) + ")"
		out += "\n"
		for i in self.windows:
			out+= str(i) + "\n"
		return out

class debugObject:
	def __init__(self, debugLevel):
		self.timeStamp = 0
		self.debugText = ""
		self.debugLevel = debugLevel
		self.timeLog = False
	
	def __str__(self):
		out = ""
		if(self.timeLog): 
			self.timeStamp = datetime.datetime.now().strftime("[%H:%M:%S]")
			out+= str(self.timeStamp) + " "
		out+= str(self.debugText)
		return out
		
	def setLevel(self, newLevel):
		self.debugLevel = newLevel
		
	def toggleTimeLog(self):
		if self.timeLog: 
			self.timeLog = False
		else:
			self.timeLog = True
		
	def write(self, debugText, level = 3):
		self.debugText = debugText
		if (int(self.debugLevel)>=int(level)): print str(self)

				

class WindowObject:
	""" This class contains the definitions of the windows for a specific run
	"""
	def __init__(self, xll, yll, xsize, ysize):
		self.xll = xll
		self.yll = yll
		self.xsize = xsize
		self.ysize = ysize
		
	def __str__(self):
		out = "lower left: (" + str(self.xll) + ", " + str(self.yll) + ") size (" + str(self.xsize) + ", " + str(self.ysize) + ")"
		return out

class combined3ColourObject:
	""" This is a class to encapsulate a single observed object with all three channel data combined
	"""
	def __init__(self, id):
		self.r = None
		self.g = None
		self.b = None
		self.id = id
		
	def setRedObject(self, object):
		self.r = object

	def setGreenObject(self, object):
		self.g = object

	def setBlueObject(self, object):
		self.b = object

	def __str__(self):
		outString = "ID: " + str(self.id) + "\n"
		
		outString+= "Red: " + str(self.r) + "\n"
		outString+= "Green: " + str(self.g) + "\n"
		outString+= "Blue: " + str(self.b)
	
		return outString
		
	def toJSON(self):
		outObject = {'id':0, 'red_id': -1, 'green_id':-1, 'blue_id':-1}
		outObject['id'] = self.id
		if (self.r!=None): outObject['red_id'] = self.r.id
		if (self.g!=None): outObject['green_id'] = self.g.id
		if (self.b!=None): outObject['blue_id'] = self.b.id
		return json.dumps(outObject)

class ObservedObject:
	""" This is a class to encapsulate a single observed object in a single channel (colour), usually a star
	"""
	def __init__(self, id):
		self.id = id                # ID for this object, unique to this run
		self.CCDchannel = 0         # which channel am I on? r = 1, g = 2, b = 3, undefined = 0 
		self.CCDportion = 0         # which half of the CCD am I on? left = 1, right = 2
		self.windowID = 0           # which window am I on? ID from 0 to n 
		self.distanceThreshold = 3.0
		self.numExposures = 0       # The number of exposures in this run for this object
		self.exposures = []         # An array containing numExposures exposure objects
		self.currentPosition = (0,0) # The object's last known x, y position
		self.lastCounts = 0
		
	def addExposure(self, x, y, counts, FWHM, MJD, exposureTime):
		""" Adds a new exposure object to this object
		"""
		exposure = ExposureObject()
		exposure.centroid = (x, y)
		exposure.counts = counts
		exposure.exposureTime = exposureTime
		exposure.MJD = MJD
		exposure.FWHM = FWHM
		self.exposures.append(exposure)
		self.numExposures+= 1
		self.currentPosition = (x, y)
		
	def resetExposures(self, exposures):
		self.exposures = exposures
		self.numExposures = len(exposures)
		
	def removeExposure(self, index):
		""" Removes an expsure by index. This is used to remove frames with bad timings (like at the start of a run in DRIFT mode)
		"""
		self.exposures.pop(index)
		print "Removing exposure at frame: ", index
		self.numExposures-= 1
		
	def calculateMeanPosition(self):
		(mx, my) = (0, 0)
		for i in range(self.numExposures):
			(x, y) = self.exposures[i].centroid
			mx += x
			my += y
		mx /= self.numExposures
		my /= self.numExposures
		
		self.meanPosition = (mx, my)
		
		return self.meanPosition
		
		
	def addExposureByObject(self, newValue, MJD):
		""" Adds a new exposure object to this object
		"""
		exposure = ExposureObject()
		exposure.centroid = ( newValue['x'], newValue['y'] ) 
		exposure.counts = newValue['counts']
		exposure.exposureTime = 0
		exposure.MJD = MJD
		exposure.FWHM = 0
		self.exposures.append(exposure)
		self.numExposures+= 1
		self.currentPosition = (newValue['x'], newValue['y'] )
		
	def isDistanceMatch(self, object):
		""" Returns -1 if object is not a match, or the distance if it is closer than the distanceThreshold
		"""
		xo = object['x'] - self.currentPosition[0]
		yo = object['y'] - self.currentPosition[1]
		distance = math.sqrt(xo*xo + yo*yo)
		if (distance>self.distanceThreshold): return -1;
		return distance;

	def isCosmicRay(self):
		""" Returns true if this object thinks it is a cosmic ray, based on the simple test that it only appears in 1 frame
		"""
		if self.numExposures == 1: return True;
		return False
		
	def isInCircle(self, x, y, radius):
		""" Returns -1 if object is not in the circle, or the distance if it is inside the circle
		"""
		xo = x - self.currentPosition[0]
		yo = y - self.currentPosition[1]
		distance = math.sqrt(xo*xo + yo*yo)
		if (distance>radius): return -1;
		return distance;	
		
	def toString(self):
		""" Returns a nicely formatted description of this object (for debug purposes)
		"""
		out = ""
		out += "ID: " + str(self.id) + " (" + str(self.currentPosition[0]) + ", " + str(self.currentPosition[1]) + ") frames: " + str(self.numExposures) + "   last counts:" + str(self.exposures[-1].counts)
		return out
		
	def __repr__(self):
		self.lastCounts = self.exposures[-1].counts
		return repr((self.id, self.lastCounts))
		
	def __str__(self):
		""" Returns a nicely formatted description of this object (for debug purposes)
		"""
		out = ""
		out += "ID: " + str(self.id) + " (" + str(self.currentPosition[0]) + ", " + str(self.currentPosition[1]) + ") frames: " + str(self.numExposures) + "   last counts:" + str(self.exposures[-1].counts)
		return out
		
	def toJSON(self):
		testObject = {'id': 0, 'x':0, 'y':0, 'data':[]}
		testObject['id'] = self.id
		testObject['x'] = self.currentPosition[0]
		testObject['y'] = self.currentPosition[1]
		exposureDataArray = []
		for c in self.exposures:
			exposureData = (c.MJD, c.counts, c.centroid[0], c.centroid[1])
			exposureDataArray.append(exposureData)
		testObject['data'] = exposureDataArray
		return json.dumps(testObject)
		
	def getDeltaXY(self):
		""" Returns the delta XY between the two most recent frames for this object or (0,0) if the object only has one frame
		"""
		if (self.numExposures<2): return (0, 0);
		deltaXY = (0, 0) 
		thisXY = self.exposures[-1].centroid
		previousXY = self.exposures[-2].centroid
		print thisXY, previousXY,
		deltaXY = (thisXY[0]-previousXY[0], thisXY[1]-previousXY[1]) 
		print deltaXY
		return deltaXY
		
	def getData(self):
		""" Returns a tuple of arrays with the MJD and counts for this object. Useful for plotting
		"""
		MJDArray = []
		CountsArray = []
		for i in self.exposures:
			MJDArray.append(i.MJD)
			CountsArray.append(i.counts)
		return (MJDArray, CountsArray) 
		
	def getMJDs(self):
		""" Returns an array with the MJD for this object. Useful for plotting
		"""
		MJDArray = []
		for i in self.exposures:
			MJDArray.append(i.MJD)
		return MJDArray 

	def getLastCounts(self):
		""" Returns an integer containing the last read counts figure for this object
		"""
		return self.exposures[-1].counts
		
	def getCountsForMJD(self, MJD):
		""" Returns counts for a certain MJD, returns 0 if not found
		"""
		counts = 0
		for e in self.exposures:
			 if e.MJD==MJD: counts = e.counts
		return counts
  
	def getCounts(self):
		""" Returns an array with the counts for this object. Useful for plotting
		"""
		CountsArray = []
		for i in self.exposures:
			CountsArray.append(i.counts)
		return CountsArray 