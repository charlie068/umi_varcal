# THIS FUNCTION WILL CALCULATE THREE STATISTICS : THE PERCENTAGE OF THE BED COVERED BY THE BAM/SAM FILE, 
# THE AVERAGE DEPTH OF THE BAM/SAM FILE ACROSS ALL POSITIONS AND THE UNIFORMITY OF THE SEQUENCING (THE 
# PERCENTAGE OF THE LOCATIONS THAT HAVE A DEPTH > 0.2*AVERGAE DEPTH)
#
# INPUT : -PILEUP            : (DICT) THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
#		  -POTENTIAL         : (INT)  NUMBER OF POTENTIAL VARIANTS FOUND
#		  -FINAL   			 : (INT)  NUMBER OF FINAL VARIANTS FOUND
#
# VALUE : NONE
#	

### import functions from the script in the same dir
from func import *

def CalculateStats(pileup, potential, final):
	print("\n")
	PrintTime('console', "\tCalculating statistics...")

	# initiate a counter for not covered positions
	not_covered = 0
	# get total number of positions in pileup
	totalPos = GetPileupLength(pileup)
	# create an empty list to stock depths
	depths = []

	# loop through the pileup items
	# for each chromosome
	for chrom , infos in pileup.items():
		# for each position | composition = counts
		for position, composition in infos.items():
			# try-except block
			# if position is covered => depth must be > 0
			try:		
				# calculate insertion counts
				n_ins_pos = 0
				n_ins_neg = 0
				for value in composition['in'].values():
					n_ins_pos += value[0]
					n_ins_neg += value[1]
				
				# calculate deletion counts
				n_del_pos = 0
				n_del_neg = 0
				for value in composition['del'].values():
					n_del_pos += value[0]
					n_del_neg += value[1]
				
				# insertions shouldn't be taken into account when calculating depth
				# n_reads = composition['A'][0]+composition['C'][0]+composition['G'][0]+composition['T'][0]+composition['A'][1]+composition['C'][1]+composition['G'][1]+composition['T'][1]+n_ins_pos+n_ins_neg+n_del_pos+n_del_neg
				n_reads = composition['A'][0]+composition['C'][0]+composition['G'][0]+composition['T'][0]+composition['A'][1]+composition['C'][1]+composition['G'][1]+composition['T'][1]+n_del_pos+n_del_neg
				# make a test to check that depth > 0
				test = 25 / n_reads
				# if test succeed => position covered => append depth to list of depths
				depths.append(n_reads)

			# if position not covered => increment not_covered counter
			except:
				not_covered += 1

	# calculate coverage
	coverage = round(float(totalPos-not_covered)/float(totalPos), 5)*100
	# calculate average depth
	avg_depth = int(round(float(sum(depths))/float(len(depths)), 0))

	# calculate uniformity
	uniform = 0
	for depth in depths:
		if depth >= 0.2*avg_depth:
			uniform += 1
	uniformity = round(float(uniform)/float(len(depths)), 2)*100


	# print out stats to console
	message = "BAM/BED coverage: "+ str(coverage)+" %"
	PrintTime('green', "\t\t"+message)
	message = "Average design depth: "+ str(avg_depth)+"x"
	PrintTime('green', "\t\t"+message)
	message = "Uniformity: "+ str(uniformity)+" %"
	PrintTime('green', "\t\t"+message)
	message = "Candidate Positions: "+ str(potential)
	PrintTime('green', "\t\t"+message)
	message = "Final Variants: "+ str(final)
	PrintTime('green', "\t\t"+message)