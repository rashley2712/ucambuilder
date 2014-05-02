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
	parser.add_argument('-o', '--output', type=str, help='Output file name. If left out, ".thumb" will be added before the .png extension')
	parser.add_argument('--flip', action='store_true', help='Flip the y-axis for cropping.')
	parser.add_argument('--autocrop', action='store_true', help='Perform an autocrop for non-zero pixels.')
	arg = parser.parse_args()

	if (arg.height==None): arg.height = arg.width
	thumbSize = (arg.width, arg.height)
	print "Opening file: ", arg.input
	if (arg.output==None):
		outputFilename = arg.input[:len(arg.input)-4] + ".thumb.png"
		print "auto generated output filename: ", outputFilename
	else:
		outputFilename = arg.output
		
	image = Image.open(arg.input)
	
	# Perform the 'crop'
	if (arg.crop):
		extents = (arg.xl, arg.yl, arg.xu, arg.yu)
		print "Cropping to: ", extents
		croppedImage = image.crop(extents)
		image = croppedImage
	
	if (arg.autocrop):
		print "Auto cropping to bounding box: ", image.getbbox()
		croppedImage = image.crop(image.getbbox())
		image = croppedImage
	
	print "New size: ", thumbSize
	image.thumbnail(thumbSize, Image.ANTIALIAS)
	image.save(outputFilename, "PNG")
	
	print "Your thumbnail is ready at:", outputFilename
