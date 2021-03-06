import numpy, json, classes
import astropy.io.fits
import os, subprocess, math, re
import Image,ImageDraw,ImageFont

def readConfigFile():
    """ Reads the config file name ucambuilder.conf and returns a configObject
    """
    configuration = classes.configObject()
    
    try:
        configfile = open("ucambuilder.conf", "r")
    except IOError:
        print "Warning: Cannot find ucambuilder.conf in the local path... will use defaults"
        return configuration

    for line in configfile:
        if(line[0]!="#"):      # Ignore lines that are comments
            if len(line.split())>=2:
                option = line.split()[0]
                value = line.split()[1]
                configuration[option] = value
	
    configfile.close()
    return configuration
    
def getRunMetaData(runName):
    """ Opens an Ultracam raw file and gets some metadata from it
    """
    rdat  = ultracam.Rdata(runFilename, startFrame, server=False)
    numFrames = rdat.ntotal()

def getUniqueID(objects):
	""" Returns a unique ID (monotonically increasing) for a list of objects that have a property called 'id'
	""" 
	newID = 0 
	for o in objects:
		if (o.id >= newID): newID = o.id + 1	
	return newID

def percentiles(data, lo, hi):
    """ Returns a normalised array where lo percent of the pixels are 0 and hi percent of the pixels are 255
    """
    max = data.max()
    dataArray = data.flatten()
    pHi = numpy.percentile(dataArray, hi)
    pLo = numpy.percentile(dataArray, lo)
    range = pHi - pLo
    scale = range/255
    data = numpy.clip(data, pLo, pHi)
    data-= pLo
    data/=scale
    return data

def measureDistance(p1, p2):
	""" Measures the 2D pythogorean distance between to (x,y) tuples
    """
	(x1, y1) = p1
	(x2, y2) = p2
	
	distance = math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) )
	
	return distance
	
    
def buildObjectsFromJSON(filename):
	""" Reads a JSON file and re-constructs the object list... returns it as an array
	"""
	config = readConfigFile()
	debug = classes.debugObject(config.DEBUG)
	JSONfile = open(filename, "r")

	wholeFileString = JSONfile.read()

	allObjectsJSON = json.loads(wholeFileString)

	objects = []

	for i in allObjectsJSON:
		ob = json.loads(i)
		# Create a new instance of ObservedObject for this object 
		newObject = classes.ObservedObject(ob['id'])
		newObject.currentPosition = (ob['x'], ob['y'])
		debug.write("ID: %d (%d, %d)"%(ob['id'], int(ob['x']), int(ob['y'])))
		dataArray = ob['data']
		for data in dataArray:
			newObject.addExposure(data[2],data[3], data[1], 0, data[0], 0)
			debug.write(data)
		objects.append(newObject)
		
	return objects
	
def createFITS(number, imageData):
    """ Writes a FITS file for the images in a frame of the CCD
    """
    filename = "/tmp/frame" + str(number) + ".fits"
    hdu = astropy.io.fits.PrimaryHDU(imageData)
    hdulist = astropy.io.fits.HDUList([hdu])
    hdulist.writeto(filename, clobber=True)
    
def saveFITSImage(imageData, filename):
    """ Writes a FITS file for the images in a frame of the CCD
    """
    hdu = astropy.io.fits.PrimaryHDU(imageData)
    hdulist = astropy.io.fits.HDUList([hdu])
    hdulist.writeto(filename, clobber=True)

def removeFITS(number):
    """ Removes a temporary FITS file once sextractor has finished with it
    """
    filename = "/tmp/frame" + str(number) + ".fits"
    os.remove(filename)

def removeCAT(number):
    """ Removes a temporary CAT file once sextractor has finished with it
    """
    filename = "/tmp/frame" + str(number) + ".cat"
    os.remove(filename)


def runSex(frameNumber):
    """ Runs a sextractor process for an image that has just been placed in a FITS file
    """
    fitsFilename = "/tmp/frame" + str(frameNumber) + ".fits"
    catFileParameter = " CATALOG_NAME /tmp/frame" + str(frameNumber) + ".cat"
    subprocess.call(["sex", fitsFilename, "-" + catFileParameter])
    
def readSexObjects(frameNumber):
	""" Reads a sextractor .cat file and returns a list of objects in that file
    """
	filename = "/tmp/frame" + str(frameNumber) + ".cat"
	sexfile = open( filename, "r" )
 
	objects = []

	for line in sexfile:
		if(line[0]!="#"):                                    # If the line is a comment, ignore it
		# sys.stdout.write(line)                          
			object = {'id': 0, 'x':0, 'y':0, 'radius':0, 'counts':0, 'flags':0}
			params = line.split()
			id = int(params[0])
			x = float(params[1])
			y = float(params[2])
			counts = float(params[3])
			mag = float(params[4])
			radius = float(params[5])
			flags = int(params[6])
			object['id'] = id
			object['x'] = x
			object['y'] = y
			object['counts'] = counts
			object['radius'] = radius
			object['flags'] = flags
			objects.append(object)
	
	sexfile.close()

	return objects
	
def filterOutCosmicRays(objects):
	""" Returns a reduced list of objects with the cosmic rays removed
	"""
	filteredList = []
	for i in objects:
		if not i.isCosmicRay():
			filteredList.append(i)
	
	return filteredList
		
def filterOutLowFrameCountObjects(objects, percentage):
	""" For a given list of objects, filters out those objects that do not appear on > percentage of frames
	"""
	maxFrames = 0
	for i in objects:
		if (i.numExposures>maxFrames): maxFrames = i.numExposures;
	
	filteredObjects = []
	threshold = maxFrames * percentage/100
	for i in objects:
		if i.numExposures>=threshold: filteredObjects.append(i)
		
	return filteredObjects

def filterOutBadTimingFrames(objects):
	""" For a given list of objects, checks for frames with weird times (like the first few frames in a DRIFT mode run, or just some random jumps in time) and removes them from the data
	This approach is very inefficient since it checks all objects individually, when timing errors are more likely to be the same for all objects. Needs to be optimised.
	"""
	
	for o in objects:
		startTime = o.exposures[0].MJD
		endTime = o.exposures[-1].MJD
		
		newExposures = []
		for f, frame in enumerate(o.exposures):
			time = frame.MJD
			difference = time - startTime
			print("MJD: " + str(frame.MJD) + " [" + str(f) +"] Time difference: " + str(difference))
			if (time<startTime): 
				print "We need to remove this exposure"
			elif (time>endTime): 
				print "We need to remove this exposure"
			else:
				print "This exposure is ok"
				newExposures.append(frame)
						
		o.resetExposures(newExposures)
				
	return objects		
		

def rejectBadObjects(objects):
    """ Returns a filtered set of objects without the sextractor bad objects
    """
    filteredObjects = []
    for ob in objects:
      if (ob['flags'] == 4): filteredObjects.append(ob)
      if (ob['flags'] == 0): filteredObjects.append(ob)
    return filteredObjects

def getTopObjects(objectList, num):
	""" Returns the brightest n objects in the field
	"""
	brightest = []
	sortedList = sorted(objectList, key=lambda object: object.lastCounts, reverse = True)
	
	brightest = sortedList[0:num]
	
	return brightest
	
def computeOffsetArray(arrayObject, dx, dy):
	""" Returns a new array with the values linearly interpolated by the offset amount
	"""
	# Start with dx
	row = arrayObject[0]
	print row
	outputRow = numpy.zeros(shape = len(row))
	
	for i in range(len(row)-1):
		newXvalue = i + dx
		print newXvalue, math.floor(newXvalue), math.ceil(newXvalue)
		y1 = row[math.floor(newXvalue)]
		y2 = row[math.ceil(newXvalue)]
		m = y2-y1
		c = y1 - m*i
		print "m,c",m,c, 
		ynew = (newXvalue) * m + c
		print "ynew", ynew
		outputRow[i] = ynew	
		
	print row
	print outputRow
	
def separateRunNameAndDate(name):
	""" Takes a string like 2013-07-21/run111 and returns a tuple with the date and the runName separated
	"""
	runs_re = re.compile(r'run[0-9]{3}')
	date_re = re.compile(r'20[0-9]{2}(-[0-9]{2}){2}')
		
	date = ""
	runName = ""
	
	d = date_re.search(name)
	if (d):
		date = d.group(0)

	r = runs_re.search(name)
	if (r):
		runName = r.group(0)
	
	return (date, runName) 
	
def addPaths(path1, path2):
	""" Adds two paths together inserting a '/' if required
	"""
	path = path1
	if (path[-1]!='/'): path+= '/';
	path+= path2
	return path
	
def writePNG(image, filename, caption = ""):
    """ Writes to a PNG file using the PIL library. Adds a caption if sent in the parameters. Also adds a .png extension if it isn't already there in 'filename'
	"""
    imageCopy = image.copy()    # We need to copy the image so we don't alter the original when adding the caption.
    config = readConfigFile()
    if (caption!=""): 
		font = ImageFont.truetype(config.FONT, 25) 
		draw = ImageDraw.Draw(imageCopy)
		if (imageCopy.mode == "L"):
			draw.text((0, 0), caption, 255, font = font)
		else: 
			draw.text((0, 0), caption, (255, 255, 255), font = font)
	
    if (filename[-4:]!=".png"): filename+= ".png"
    imageCopy.save(filename, "PNG")
	
