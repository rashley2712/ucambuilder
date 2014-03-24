import numpy

class wcsSolution:

	def __init__(self):
		self.equinox = 2000.0
		self.SIPOrder = 3
		self.pixReference = (512.0, 512.0)
		self.ra = 0.0
		self.dec = 0.0 
		self.raDeg = 0.0
		
		self.linearTransform = [ [0.0, 0.0], [0.0, 0.0] ]
		
	def setSolution(self, equinox, refCoord, refPixel, CD):
		self.equinox = equinox
		self.raDeg = refCoord[0]
		self.dec = refCoord[1]
		self.ra = self.raDeg / 15.
		self.linearTransform = CD
		self.pixReference = refPixel
		
	def getSexagesimal(self):
		ra = self.ra
		hours = int(ra)
		minutes = (ra - int(ra)) * 60
		seconds = (minutes - int(minutes)) * 60
		
		dec = self.dec
		decDegrees = int(dec)
		decMinutes = (dec - int(dec)) * 60
		decSeconds = (decMinutes - int(decMinutes)) * 60
		
		outString = "RA: %2d:%2d:%2.1f"%(hours, minutes, seconds)
		outString+= " DEC: %2d:%2d:%2.3f"%(dec, decMinutes, decSeconds)
		return outString
		
	def __str__(self):
		outString = "RA: %.4f DEC: %.4f --- "%(self.ra, self.dec)
		outString+= self.getSexagesimal()
		return outString
		
	def getWorldCoord(self, position):
		u = position[0] - self.pixReference[0]
		v = position[1] - self.pixReference[1]
		CD = numpy.array(self.linearTransform)
		pixel = numpy.array((u,v))
		world = numpy.dot(CD,pixel)
		world = (world[0] + self.raDeg, world[1] + self.dec)
				
		return world
		
