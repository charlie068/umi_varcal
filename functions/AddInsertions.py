
# THIS FUNCTION ALLOWS TO PARSE A CIGAR SEGMENT CONTAINING ONLY INSERETD BASES AND INCREMENT
# A SPECIFIC PILEUP DICTIONNARY ACCORDINGLY
#
# INPUT : -PILEUP            : (DICT) THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
# 		  -UMI               : (STR)  UMI SEQUENCE OF THE READ
# 		  -CHROM             : (STR)  THE CHROMOSOME MAPPED TO THE READ
# 		  -POSITION          : (INT)  THE POSITION OF THE INSERTION SITE IN THE READ
# 		  -SEQ               : (STR)  THE SEQUENCE OF THE READ
# 		  -STRAND            : (INT)  THE STRAND OF THE READ (0 = FORWARD | 1 = REVERSE)
# 		  -CURSOR_POS        : (INT)  THE CURRENT POSITION IN THE READ
# 		  -CURSOR_SEQ        : (INT)  THE CURRENT POSITION IN THE SEQUENCE OF THE READ
# 		  -MAXX              : (INT)  THE LENGTH OF THE CIGAR 'I' ELEMENT
# 		  -QUAL              : (STR)  THE QUALITY STRING OF THE READ
# 		  -QUALS             : (DICT) A DICTIONNARY FOR THE CONVERSION OF THE QUALITIES FROM ASCII TO INT
# 		  -MIN_BASE_QUALITY  : (INT)  MINIMUM QUALITY SCORE OF THE BASE FOR IT TO ADDED TO THE PILEUP DICTIONNARY
#		  -ALL_UMIS          : (DICT) A DICTIONNARY CONTAINING UMIS INDEXES
#
# VALUE : -POSITION          : (INT)  THE FINAL POSITION THAT'S BEEN PARSED IN THE READ
#		  -CURSOR_POS        : (INT)  THE CURRENT POSITION IN THE READ
# 		  -CURSOR_SEQ        : (INT)  THE CURRENT POSITION IN THE SEQUENCE OF THE READ

def AddInsertions(pileup, umi, chrom, position, seq, strand, cursor_pos, cursor_seq, maxx, qual, quals, MIN_BASE_QUALITY, ALL_UMIS):

	# create an empty string to contain the inserted sequence
	# create a mean qscore for the inserted sequence 
	inserted_seq = ""
	inserted_qscore = 0


	# try because position in read could not be in BED file
	# therefore, the position could not be in the PILEUP dictionnary
	# if position not in BED, skip the read
	try:
		# testing if position + 1 if in the BED file
		# position is the position of the last matched/mismatched position
		# returned by the AddMatches function, therefore the insertion 
		# occurs at position+1
		test = pileup[chrom][position+1]



		# looping through the inserted sequence
		# maxx = length of the "I" cigar element => length of the inserted_sequence
		for i in range(0, maxx):
			# get base at i
			base = seq[cursor_seq]
			# get base qual at i
			baseQual = qual[cursor_seq]
			# convert base qual to base qscore
			base_qscore = quals[baseQual]

			# add the base at i in seq to the inserted sequence
			inserted_seq += base
			# increment the inserted_seq qscore with the quality 
			# of the inserted base
			inserted_qscore += base_qscore
			# advance in the read sequence
			cursor_seq += 1



		# calculate qscore for the inserted sequence corresponding
		# to the mean of the qscores of the bases in the inserted seq
		inserted_qscore = inserted_qscore/len(inserted_seq)


			
		# check if the quality of the inserted sequence >= minimum base quality
		if inserted_qscore >= MIN_BASE_QUALITY:

			# if this insertion was already seen at this location
			try:
				# try to increment the corresponding counter in the PILEUP dictionnary
				pileup[chrom][position+1]['in'][inserted_seq][strand] += 1
				# add the umi to the corresponding set specific to the inserted sequence
				pileup[chrom][position+1]['in'][inserted_seq][2].append(umi)
				# increment the qScore of this position
				pileup[chrom][position+1]['base_error_probability'] += inserted_qscore

			# if this is the first time we see this inserted sequence at this position, 
			# an entry should be created with its appropriated structure (2 counters + set)  
			except:
				# create the entry in PILEUP dictionnary
				pileup[chrom][position+1]['in'][inserted_seq] = [0, 0, [umi]]
				# increment the corresponding counter in the PILEUP dictionnary
				pileup[chrom][position+1]['in'][inserted_seq][strand] += 1
				# increment the qScore of this position
				pileup[chrom][position+1]['base_error_probability'] += inserted_qscore


			cursor_pos += 1
			cursor_seq += 1

	
	# if position not in BED, skip the read
	except:
		cursor_seq += maxx

	# if seq == "GAATTTAGTAGCCCCCAGCACTTTAGTGATGTTTGGTAACAGCTGTAGAAGCGCTTATGCTATGGACTGAAGCTAGACCACCTGAGTTAGTATTCTGGTCCCTTCACTTACTAGCTCTGTGAACATGGGCAAATACTTAATTTCAT":
	# 	print('i'+str(cursor_seq))
	
	# return position, cursor_seq and cursor_pos to continue from where we left off
	return [position, cursor_seq, cursor_pos]


