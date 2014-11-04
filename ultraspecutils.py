import json
import ultraspecClasses

def loadApertures(jsonFilename):
	jsonFile = open(jsonFilename, 'r')
	wholeFileString = jsonFile.read()
	parsedApertures = json.loads(wholeFileString)
	allApertures = []
	for apertureObject in parsedApertures:
		position = (apertureObject['x'], apertureObject['y'])
		radius = 10
		apertureID = apertureObject['id']
		aperture = ultraspecClasses.aperture(apertureID, position, radius)
		allApertures.append(aperture)
		
	return allApertures
		
def determineFullFrameSize(windows):
	leftestPixel = 1057
	rightestPixel = 0
	topestPixel = 0 
	bottomestPixel = 1040
	for w in windows:
		if w.xll/w.xbin < leftestPixel: leftestPixel = w.xll/w.xbin
		if w.yll/w.ybin < bottomestPixel: bottomestPixel = w.yll/w.ybin
		if (w.xll/w.xbin + w.nx) > rightestPixel: rightestPixel = w.xll/w.xbin + w.nx
		if (w.yll/w.ybin + w.ny) > topestPixel: topestPixel = w.yll/w.ybin + w.ny

	return leftestPixel, bottomestPixel, rightestPixel, topestPixel
	