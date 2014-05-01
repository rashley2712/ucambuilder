#!/usr/bin/env python
import Image
import argparse

if __name__ == "__main__":
	
		
	parser = argparse.ArgumentParser(description='Creates a thumbnail of the input image with the new dimensions as specified by the input parameters')
	parser.add_argument('input', type=str, help='Input image filename')
	parser.add_argument('-c', '--crop', action='store_true', help='Toggle to crop to maxExtents')
	parser.add_argument('-w', '--width', default=128, type=int, help='New x dimension... default is 128')
	parser.add_argument('-e', '--height', type=int, help='New y dimension... default is the width')
	parser.add_argument('-xl', default=1, type=int, help='lower x extent (for crop, if requested) - default is 1')
	parser.add_argument('-xu', default=1080, type=int, help='upper x extent (for crop, if requested) - default is 1080')
	parser.add_argument('-yl', default=1, type=int, help='lower y extent (for crop, if requested) - default is 1')
	parser.add_argument('-yu', default=1040, type=int, help='upper y extent (for crop, if requested) - default is 1040')
	arg = parser.parse_args()

	if (arg.height==None): arg.height = arg.width
	size = (arg.width, arg.height)
	print "Opening file: ", arg.input
	outputFilename = arg.input[:len(arg.input)-4] + ".thumb.png"
	print "output: ", outputFilename
	
	image = Image.open(arg.input)
	
	# Perform the 'crop'
	if (arg.crop):
		extents = (arg.xl, arg.yl, arg.xu, arg.yu)
		print "Cropping to: ", extents
		croppedImage = image.crop(extents)
		image = croppedImage
	
	print "Finished crop"
	print "New size: ", size
	image.thumbnail(size, Image.ANTIALIAS)
	image.save(outputFilenamye, "PNG")
