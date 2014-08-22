#!/usr/bin/env python
import ultracamutils
import sys, argparse, subprocess, re, os
import classes

if __name__ == "__main__":
	
		
	parser = argparse.ArgumentParser(description='Puts a daybuilder.py into the CoWS scheduler')
	parser.add_argument('date', type=str, help='Ultracam date  [eg 2013-07-21]')
	arg = parser.parse_args()
	
	cowsFilename = arg.date + ".pbs"
	
	cowsFile = file(cowsFilename, 'w')
	
	cowsFile.write("#!/bin/bash\n")
	cowsFile.write("#PBS -l nodes=1:ppn=1,pvmem=1024mb,walltime=12:00:00\n")
	cowsFile.write("#PBS -V\n\n")
	cowsFile.write("setenv PBS_O_WORKDIR /storage/astro2/phrnaw/ucambuilder\n")
	cowsFile.write("cd $PBS_O_WORKDIR\n\n")
	cowsFile.write("/storage/astro2/phrnaw/code/ucambuilder/daybuilder.py " + arg.date + " -r\n\n")
	cowsFile.close()
	
	cowsCommand = ['qsub']
	cowsCommand.append('-qserial')
	cowsCommand.append(cowsFilename)
	subprocess.call(cowsCommand)
	
	
	
