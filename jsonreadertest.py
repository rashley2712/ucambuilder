#!/usr/bin/env python
import json
import classes, ucamObjectClass


filename = "run114_objects.json"

""" Reads a JSON file and re-constructs the object list... 
"""
JSONfile = open(filename, "r")
wholeFileString = JSONfile.read()
allObjectsJSON = json.loads(wholeFileString)
objects = []
for i in allObjectsJSON:
	ob = json.loads(i)
	#print ob
	"""   Create a new instance of ObservedObject for this object"""
	newObject = ucamObjectClass.colourObject(ob['id'])
	
	meanPosition = (ob['meanPosition'])
	
	print meanPosition['r']
	"""     newObject.currentPosition = (ob['x'], ob['y'])
         dataArray = ob['data']
         for data in dataArray:
                 newObject.addExposure(data[2],data[3], data[1], data[4], data[0], 0)
         objects.append(newObject)
	"""

JSONfile.close()






