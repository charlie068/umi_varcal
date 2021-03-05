
# THIS FUNCTION IS CALLED TO SPLIT THE INITIAL SAM FILE INTO N SAM SUBFILES, N CORRESPONDING
# TO THE NUMBER OF CORES USED TO EXECUTE THE PROGRAM. THEREFORE, EACH CORE WILL ANALYZE A 
# SUBFILE, ALLOWING THEREFORE ALL THE SUBFILES TO BE ANALYZED AT THE SAME TIME BY ALL THE 
# PRECISED CORES TO BE USED.
#
# INPUT : 
# 		  -FILENAME                     : (STR)    THE PATH OF THE INITIAL SAM FILE
# 		  -TOTALLINES                   : (FLOAT)  TOTAL LINES IN THE FILE
# 		  -CORES                        : (INT)    NUMBER OF CORES TO BE USED
#
# VALUE : -SUBFILES                     : (LIST)   A LIST CONTAINING THE NAMES OF THE CREATED SUBFILES
#		
#

import os
import sys

from func import PrintTime, PrintProgress


def PreprocessReads(fileName, totalLines, cores):

	# if only one core is to be used, no need to split the file
	# return the initial file
	if cores == 1:
		return [fileName]




	print("\n")
	PrintTime('console', "\tPreprocessing Reads...")


	# create an empty list to contain the created sub files names
	subFiles = []

	# calculate at the interval at which a new subfile is to be created
	# this allows for the subfiles to be created with approximately the
	# same size
	interval =  int(totalLines/cores)
	interval = interval if interval % 2 == 0 else interval - 1


	# set counters
	i = 0
	nLine = 0.0
	lastProgress = 0.0

	# read the file line by line
	for line in open(fileName, 'r'):

		# display progress
		lastProgress = PrintProgress(nLine+1, totalLines, lastProgress)

		# if interval value is reached
		if nLine % interval == 0:
			# if this subfile is not the last subfile
			if i != cores:
				# define the name of the subfile
				outName = fileName.replace('.sam', '_'+str(i)+'.sam')
				# create the output file
				output = open(outName, 'w')
				# add the subfile name to the subfiles list
				subFiles.append(outName)
				# increment the subfiles counter
				i += 1
		# write the line in the i subfile
		output.write(line)
		# increment the line counter
		nLine += 1


	
	print("\n")
	PrintTime('console', "\tDone")

	# return the subfiles list 
	return subFiles