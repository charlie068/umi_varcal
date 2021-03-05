
# THIS FUNCTION ALLOWS TO PARSE A CIGAR SEGMENT CONTAINING ONLY DELETED BASES AND INCREMENT
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
# 		  -MAXX              : (INT)  THE LENGTH OF THE CIGAR 'D' ELEMENT
#		  -ALL_UMIS          : (DICT) A DICTIONNARY CONTAINING UMIS INDEXES
#
# VALUE : -POSITION          : (INT)  THE FINAL POSITION THAT'S BEEN PARSED IN THE READ
#		  -CURSOR_POS        : (INT)  THE CURRENT POSITION IN THE READ
# 		  -CURSOR_SEQ        : (INT)  THE CURRENT POSITION IN THE SEQUENCE OF THE READ
#

def AddDeletions(pileup, umi, chrom, start, seq, strand, cursor_pos, cursor_seq, maxx, ALL_UMIS):

	# create a cursor to move in the deleted sequence 
	del_cursor = 0

	# for each position between limits
	for position in range(start+cursor_pos, start+cursor_pos+maxx):
		
		# try because base position in read could not be in BED file
		# therefore, the position could not be in the PILEUP dictionnary
		# if position not in BED, skip the read
		try:
			# test if chrom and position are in the PILEUP dictionnary
			test = pileup[chrom][position]
			
			# if this deletion was already seen at this location
			# maxx-del_cursor corresponds to the length of the deleted sequence
			# this number will serve as an entry in the PILEUP dictionnary for
			# deleted sequences
			try:
				# increment the corresponding counter in the PILEUP dictionnary
				pileup[chrom][position]['del'][maxx-del_cursor][strand] += 1
				# add the umi to the corresponding set specific to the deleted sequence
				pileup[chrom][position]['del'][maxx-del_cursor][2].append(umi)

			# if this is the first time we see this deleted sequence at this position, 
			# an entry should be created with its appropriated structure (2 counters + set)  
			except:
				# create the entry in PILEUP dictionnary
				pileup[chrom][position]['del'][maxx-del_cursor] = [0, 0, [umi]]
				# increment the corresponding counter in the PILEUP dictionnary
				pileup[chrom][position]['del'][maxx-del_cursor][strand] += 1


			
			# advance in the deleted sequence if the deletion length > 1
			del_cursor += 1


		# if position not in BED, skip the read
		except:
			pass	


	# increment the cursor_pos with the length of the cigar element when done adding it
	cursor_pos += maxx

	# return position, cursor_seq and cursor_pos to continue from where we left off
	return [position, cursor_seq, cursor_pos]

