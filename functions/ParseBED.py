
# THIS FUNCTION ALLOWS TO PARSE A CIGAR SEGMENT CONTAINING ONLY DELETED BASES AND INCREMENT
# A SPECIFIC PILEUP DICTIONNARY ACCORDINGLY
#
# INPUT : 
# 		  -BED               : (STR)  THE PATH OF THE BED FILE
#
# VALUE : -PILEUP            : (DICT) A DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
#		



def ParseBED(bed):	

	# create an empty dictionnary for pileup counters
	pileup = {}

	# read the bed file line by line
	for line in open(bed):

		# split the line using tab as delimiter
		line = line.split('\t')

		if len(line) >= 3:
			# first element is the chromosome
			chrom = line[0]

			# second element is the start of the region
			# third element is the end of the region
			# create a list containing start and end 
			# sort to ensure that start < end
			limits = [int(line[1]), int(line[2])]
			limits.sort()

			# for each position in the [start;end] interval
			for pos in range(limits[0], limits[1]+1):
				
				# try to add the position to the pileup for the specified chromosome
				# if the chromosome was already aded to the pileup
				try:
					# add the position to the chrom with the default counters
					pileup[chrom][pos] = { 'A': [0, 0, []], 'C': [0, 0, []], 'G': [0, 0, []], 'T': [0, 0, []], 'in': {}, 'del': {}, 'base_error_probability': 0 } 
				except:
					# add the chrom to the pileup with the first position and the default counters
					pileup[chrom] = { pos: { 'A': [0, 0, []], 'C': [0, 0, []], 'G': [0, 0, []], 'T': [0, 0, []], 'in': {}, 'del': {}, 'base_error_probability': 0 } }



	# the pileup dictionnary contains :
	# for each chromosome:
	# 		for each position:
	#		{
	#			'A'   : 0 | 0 | list()
	#			'C'   : 0 | 0 | list()
	#			'T'   : 0 | 0 | list()
	#			'G'   : 0 | 0 | list()
	#			'in'  :     {}
	#			'del' :     {}
	#		}
	# for the four normal bases, this structure 0 | 0 | set() is a list containing 
	# a counter for the forward strand, a counter for the reverse strand, and a set
	# that will contain the UMIs of the reads. 
	# a set was used instead of a list to contain unique UMIs only  
	#
	# the insertion dictionnary will be updated each time an insertion is encountered
	# and will end up looking like this:
	# {
	# 	'inserted_seq_1' : N1_+ | N1_- | list()
	# 	'inserted_seq_2' : N2_+ | N2_- | list()
	# 	'inserted_seq_n' : Nn_+ | Nn_- | list()
	# }
	#
	# the deletion dictionnary will be updated each time a deletion is encountered
	# and will end up looking like this:
	# {
	# 	len of the deleted segment_1 : N1_+ | N1_- | list()
	# 	len of the deleted segment_2 : N2_+ | N2_- | list()
	# 	len of the deleted segment_n : Nn_+ | Nn_- | list()
	# }	
	# 
	# the base_error_probability by bases qscores at each position
	# at the end, this value will be divided by depth to get the mean qscore a the location
	# the the mean qscore will be used to estimate the noise at the location

	# sort pileup
	pileup = SortPileup(pileup)
	
	# return the pileup
	return pileup


