#!/usr/bin/env python
import rashley_utils as utils
import sys, subprocess 
import classes
import numpy as np
import astropy.io.fits
import matplotlib.pyplot as plt

if (__name__ == "__main__"):
	
	config = utils.readConfigFile()
	debug = classes.debugObject(config.DEBUG)
	
	headFilename = "None"
	fitsFilename = "None"
	outputFilename = "None"
	imagePreview = False
	
	for i in sys.argv[1:]:
		if i[:2]=="-i": 
			debug.write("input .head file:" + str(i[2:]), level = 1)
			headFilename = str(i[2:])
		if i[:2]=="-f": 
			debug.write("input .fits file:" + str(i[2:]), level = 1)
			fitsFilename = str(i[2:])
		if i[:2]=="-p": 
			imagePreview = True
		if i[:2]=="-o": 
			debug.write("output file:" + str(i[2:]), level = 1)
			outputFilename = str(i[2:])
		
			
	if (headFilename=="None"):
		debug.write("Please give me an input .head file:  eg: -i[file.head]", level = 1)
		sys.exit()
	if (fitsFilename=="None"):
		debug.write("Please give me an input .fits file:  eg: -f[file.fits]", level = 1)
		sys.exit()
	if (outputFilename=="None"):
		outputFilename = fitsFilename[:-5] + "-modified-headers.fits"

		
	headFile = open(headFilename, "r")
	
	FITSHeaders = []
	for line in headFile:
		FITSHeader = {'keyword': "", 'value': "", 'comment': ""} 
		if len(line)>8:
			if line[8]=='=':
				FITSHeader['keyword'] =  line[0:8]
				#print "we have values", line[9:30]
				FITSHeader['value'] = line[9:30]
				if len(line)>33:
					#print "with comment", line[33:]
					FITSHeader['comment'] = line[33:]
					
				FITSHeaders.append(FITSHeader)
		
				print "keyword: ", FITSHeader['keyword'], " value: ", FITSHeader['value'], " comment: ", FITSHeader['comment']
				
	headFile.close()
	
	FITSInputhdulist = astropy.io.fits.open(fitsFilename)
	
	print FITSInputhdulist.info()
	
	print FITSInputhdulist[0].data	
	imageData = FITSInputhdulist[0].data
	
	if imagePreview:
		imgplot = plt.imshow(utils.percentiles(imageData, 20, 98), cmap='gray', interpolation='nearest')
		plt.show()

	print "Doing some more stuff"
	FITSInputhdulist.close()
	
	debug.write("Writing FITS file: " + outputFilename, level=1)
	
	prihdr = astropy.io.fits.Header()
	prihdr['COMMENT'] = "This file created by makefits.py from the Ultracam pipeline."
	prihdr['OBJECT'] = "KOI-823"
	for h in FITSHeaders:
		print "Adding: ", h["keyword"], ":", h["value"]
		keyword = str(h["keyword"])
		try: 
			prihdr[keyword] = float(h["value"])
		except ValueError:
			value = h["value"].split('\'')
			if len(value)>1: 
				print value[1]
				prihdr[keyword] = value[1]
	
	hdu = astropy.io.fits.PrimaryHDU(imageData, header=prihdr)
	hdulist = astropy.io.fits.HDUList([hdu])
	hdulist.writeto(outputFilename, clobber=True)
