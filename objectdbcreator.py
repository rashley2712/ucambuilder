#!/usr/bin/env python

import argparse

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description='Reads the Ultracam .dat files and identifies and tracks the objects.')
	parser.add_argument('runname', type=str, help='Ultracam run name  [eg 2013-07-21/run010]')
	parser.add_argument('-p', '--preview', action='store_true', help='Toggle Matplotlib image previews')
	parser.add_argument('-c', '--configfile', default='')

	arg = parser.parse_args()
	print arg.runname
	print "Previews:", arg.preview
	
	
