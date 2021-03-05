
# THIS FUNCTION ALLOWS TO PARSE A CIGAR SEGMENT CONTAINING ONLY MATCHED/MISMATCHED BASES AND INCREMENT
# A SPECIFIC PILEUP DICTIONNARY ACCORDINGLY
#
# INPUT : -PILEUP            : (DICT) THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
# 		  -UMI               : (STR)  UMI SEQUENCE OF THE READ
# 		  -CHROM             : (STR)  THE CHROMOSOME MAPPED TO THE READ
# 		  -START             : (INT)  THE START POSITION OF THE READ
# 		  -SEQ               : (STR)  THE SEQUENCE OF THE READ
# 		  -STRAND            : (INT)  THE STRAND OF THE READ (0 = FORWARD | 1 = REVERSE)
# 		  -CURSOR_POS        : (INT)  THE CURRENT POSITION IN THE READ
# 		  -CURSOR_SEQ        : (INT)  THE CURRENT POSITION IN THE SEQUENCE OF THE READ
# 		  -MAXX              : (INT)  THE LENGTH OF THE CIGAR 'M' ELEMENT
# 		  -QUAL              : (STR)  THE QUALITY STRING OF THE READ
# 		  -QUALS             : (DICT) A DICTIONNARY FOR THE CONVERSION OF THE QUALITIES FROM ASCII TO INT
# 		  -MIN_BASE_QUALITY  : (INT)  MINIMUM QUALITY SCORE OF THE BASE FOR IT TO ADDED TO THE PILEUP DICTIONNARY
#		  -ALL_UMIS          : (DICT) A DICTIONNARY CONTAINING UMIS INDEXES
#
# VALUE : -POSITION          : (INT)  THE FINAL POSITION THAT'S BEEN PARSED IN THE READ
#		  -CURSOR_POS        : (INT)  THE CURRENT POSITION IN THE READ
# 		  -CURSOR_SEQ        : (INT)  THE CURRENT POSITION IN THE SEQUENCE OF THE READ


def AddMatches(pileup, umi, chrom, start, seq, strand, cursor_pos, cursor_seq, maxx, qual, quals, MIN_BASE_QUALITY, ALL_UMIS):

	# for each position between limits
	for position in range(start+cursor_pos, start+cursor_pos+maxx):

		# try because base position in read could not be in BED file
		# therefore, the position could not be in the PILEUP dictionnary
		# if position not in BED, skip the read
		try:	

			# get base, base_qual and base_qscore
			base = seq[cursor_seq]
			baseQual = qual[cursor_seq]
			base_qscore = quals[baseQual]

				
			# check if the quality of the base >= minimum base quality
			if base_qscore >= MIN_BASE_QUALITY:
				# increment the corresponding counter in the PILEUP dictionnary
				pileup[chrom][position][base][strand] += 1
				# add the umi to the corresponding set specific to the base in the PILEUP dictionnary
				pileup[chrom][position][base][2].append(umi)
				# increment the qScore of this position
				pileup[chrom][position]['base_error_probability'] += base_qscore
			

		# if position not in BED, skip the read
		except:
			pass	

		# increment the cursor_seq to allow advancing in the read sequence
		cursor_seq += 1


	# increment the cursor_pos with the length of the cigar element when done adding it
	cursor_pos += maxx


	# return position, cursor_seq and cursor_pos to continue from where we left off
	return [position, cursor_seq, cursor_pos]



