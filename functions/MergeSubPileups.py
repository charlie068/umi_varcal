
# THIS FUNCTION ALLOWS TO MERGE MULTIPLE SUBPILEUPS TO CREATE ONE WHOLE PILEUP.
#
# INPUT : -PILEUP            : (DICT) THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
#       : -PILEUPS           : (LIST) A LIST CONTAINING ALL THE SUBPILEUPS TO BE MERGED
#       : -SUBFILES          : (LIST) THE LIST OF SUBFILES NAMES
#
# VALUE : -POSITION          : (INT)  THE FINAL POSITION THAT'S BEEN PARSED IN THE READ
#		  -CURSOR_POS        : (INT)  THE CURRENT POSITION IN THE READ
# 		  -CURSOR_SEQ        : (INT)  THE CURRENT POSITION IN THE SEQUENCE OF THE READ
#	


import os
from func import *

def MergeSubPileups(pileup, pileups, subFiles, OUTPUT):

	# check that the number of subfiles == the number of given cub pileups
	# the 2 number must be equal. If not, an error occured in the creation
	# of some pileups ==> fatal error ==> script exists with error
	if len(pileups) != len(subFiles):
		print('\n')
		PrintTime("error", "\tError while attempting to merge pileups : lengths differ!\n\t\t\tExiting...")

	# if the pileup lists == 1 ==> no merging has to be done ==> return the pileup directly
	if len(pileups) == 1:
		return pileups[0]
	

	### remove subpileups after loading them
	for subFile in subFiles:
		samName = subFile.split("/")[-1]
		pileupFile = OUTPUT+"/"+samName.replace('.sam', '.pileup')
		os.remove(pileupFile)



	# loop through the sub pileups list
	for p in pileups:

		# for each chromosome
		for chrom in p.keys():

			# for each position and its counters
			for position, composition in p[chrom].items():
				
				# increment the A forward counters
				pileup[chrom][position]['A'][0] += p[chrom][position]['A'][0]
				# increment the A reverse counters
				pileup[chrom][position]['A'][1] += p[chrom][position]['A'][1]
				# add the umis to the A unique umis set 
				pileup[chrom][position]['A'][2] += p[chrom][position]['A'][2]


				# increment the A forward counters
				pileup[chrom][position]['C'][0] += p[chrom][position]['C'][0]
				# increment the C reverse counters
				pileup[chrom][position]['C'][1] += p[chrom][position]['C'][1]
				# add the umis to the C unique umis set	
				pileup[chrom][position]['C'][2] += p[chrom][position]['C'][2]


				# increment the G forward counters
				pileup[chrom][position]['G'][0] += p[chrom][position]['G'][0]
				# increment the G reverse counters
				pileup[chrom][position]['G'][1] += p[chrom][position]['G'][1]
				# add the umis to the G unique umis set	
				pileup[chrom][position]['G'][2] += p[chrom][position]['G'][2]


				# increment the T forward counters
				pileup[chrom][position]['T'][0] += p[chrom][position]['T'][0]
				# increment the T reverse counters
				pileup[chrom][position]['T'][1] += p[chrom][position]['T'][1]
				# add the umis to the T unique umis set	
				pileup[chrom][position]['T'][2] += p[chrom][position]['T'][2]


				# add the insertions on the sub pileup to the big pileup
				# for each insertion in the insertion dictionnary
				for ins in p[chrom][position]['in'].keys():
					# try-except block
					try:
						# if the insertion is already present ==> increment this insertion counters
						pileup[chrom][position]['in'][ins][0] += p[chrom][position]['in'][ins][0]
						pileup[chrom][position]['in'][ins][1] += p[chrom][position]['in'][ins][1]
						pileup[chrom][position]['in'][ins][2] += p[chrom][position]['in'][ins][2]
					except:
						# if insertion is not present in the insertions dict ==> it has to be inserted
						# with its own counters and unique umis set
						pileup[chrom][position]['in'][ins] = [p[chrom][position]['in'][ins][0], p[chrom][position]['in'][ins][1], p[chrom][position]['in'][ins][2]]


				# add the deletions on the sub pileup to the big pileup
				# for each deletion in the deletion dictionnary
				for dele in p[chrom][position]['del'].keys():
					# try-except block
					try:
						# if the deletion is already present ==> increment this deletion counters
						pileup[chrom][position]['del'][dele][0] += p[chrom][position]['del'][dele][0]
						pileup[chrom][position]['del'][dele][1] += p[chrom][position]['del'][dele][1]
						pileup[chrom][position]['del'][dele][2] += p[chrom][position]['del'][dele][2]
					except:
						# if deletion is not present in the deletions dict ==> it has to be inserted
						# with its own counters and unique umis set
						pileup[chrom][position]['del'][dele] = [p[chrom][position]['del'][dele][0], p[chrom][position]['del'][dele][1], p[chrom][position]['del'][dele][2]]


				# increment the total base quality scores 
				pileup[chrom][position]['base_error_probability'] += p[chrom][position]['base_error_probability']
	
	# return the final whole pileup
	return pileup