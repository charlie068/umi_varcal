
# THIS FUNCTION WILL WRITE THE PILEUP DICTIONNARY IN A CSV FILE IN THE OUTPUT DIRECTORY 
#
# INPUT : 
#         -PILEUP            : (DICT) THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
#         -SAM               : (STR)  THE PATH OF THE SAM FILE
#         -OUTPUT            : (OUTPUT) THE DIRECTORY IN WHICH THE CSV FILE IS TO BE CREATED
#
# VALUE : NONE
#	

import os
import sys

from func import *

def PileupToCSV(pileup, SAM, OUTPUT):
	
	print("\n")
	PrintTime('console', "\tConverting PILEUP to CSV...")

	# create a csv file in the output dir with the same name as the SAM file
	csv = open(OUTPUT+"/"+SAM.replace(".sam", ".pileup.csv").split("/")[-1], "w")
	# create the header of the file
	header = "index,chr,position,A,C,G,T,in,del,total_reads,qScore,base_error_probability,nb_unique_umi,unique_umi_list"
	# write the header into the file
	csv.write(header+"\n")

	# counters to track progress
	currentPos = 1.0
	lastProgress = 0.0
	totalPos = GetPileupLength(pileup)


	# loop through the PILEUP dictionnary
	# for each chromosome
	for chrom , infos in pileup.items():
		# for each position | composition = counts
		for position, composition in infos.items():

			# function that displays progress efficiently
			lastProgress = PrintProgress(currentPos, totalPos, lastProgress)

			# create a set that will contain unique umis
			unique_umis = set()
			
			# get n insertions at this location for the + and the - strands
			# update the unique_umis set with the insertions umis
			n_ins_pos = 0
			n_ins_neg = 0
			for value in composition['in'].values():
				n_ins_pos += value[0]
				n_ins_neg += value[1]
				unique_umis.update(value[2])
			

			# get n deletions at this location for the + and the - strands
			# update the unique_umis set with the deletions umis
			n_del_pos = 0
			n_del_neg = 0
			for value in composition['del'].values():
				n_del_pos += value[0]
				n_del_neg += value[1]
				unique_umis.update(value[2])

			# get total bases count on the + strand
			total_reads_pos = composition['A'][0]+composition['C'][0]+composition['G'][0]+composition['T'][0]
			
			# get total bases count on the - strand
			total_reads_neg = composition['A'][1]+composition['C'][1]+composition['G'][1]+composition['T'][1]
			
			# get total bases count
			total_reads = total_reads_pos+total_reads_neg
			
			# update the unique_umis set with the 4 bases umis
			for base in ['A', 'C', 'G', 'T']:
				unique_umis.update(composition[base][2])

			
			# create the line to be written in the CSV file for each position  and write it to the file
			line = [chrom+"|"+str(position), chrom, str(position), str(composition['A'][0])+"/"+str(composition['A'][1]), str(composition['C'][0])+"/"+str(composition['C'][1]), str(composition['G'][0])+"/"+str(composition['G'][1]), str(composition['T'][0])+"/"+str(composition['T'][1]),str(n_ins_pos)+"/"+str(n_ins_neg), str(n_del_pos)+"/"+str(n_del_neg), str(total_reads), str(composition['qScore']), str(composition['base_error_probability']), str(len(unique_umis)), ";".join(unique_umis)]
			line = ",".join(line)
			csv.write(line+"\n")

			# increment counter to track progress			
			currentPos += 1.0


	# close csv file
	csv.close()

	# print done
	print("\n")
	PrintTime('console', "\tDone")