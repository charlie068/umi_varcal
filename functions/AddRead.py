
# THIS FUNCTION ALLOWS TO PARSE A CIGAR SEGMENT CONTAINING ONLY DELETED BASES AND INCREMENT
# A SPECIFIC PILEUP DICTIONNARY ACCORDINGLY
#
# INPUT : -PILEUP            : (DICT) THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
# 		  -UMI               : (STR)  UMI SEQUENCE OF THE READ
# 		  -STRAND            : (INT)  THE STRAND OF THE READ (0 = FORWARD | 1 = REVERSE)
# 		  -CHROM             : (STR)  THE CHROMOSOME MAPPED TO THE READ
# 		  -START             : (INT)  THE START POSITION OF THE READ
#		  -CIGAR             : (STR)  THE CIGAR SEQUENCE OF THE READ
# 		  -SEQ               : (STR)  THE SEQUENCE OF THE READ
# 		  -QUAL              : (STR)  THE QUALITY STRING OF THE READ
# 		  -QUALS             : (DICT) A DICTIONNARY FOR THE CONVERSION OF THE QUALITIES FROM ASCII TO INT
# 		  -MIN_BASE_QUALITY  : (INT)  MINIMUM QUALITY SCORE OF THE BASE FOR IT TO ADDED TO THE PILEUP DICTIONNARY
#		  -ALL_UMIS          : (DICT) A DICTIONNARY CONTAINING UMIS INDEXES
#
# VALUE : NONE
#		  

import re
from AddMatches import *
from AddInsertions import *
from AddDeletions import *

def AddRead(pileup, umi, strand, chrom, start, cigar, seq, qual, quals, MIN_BASE_QUALITY, ALL_UMIS):


	# split the cigar sequence using any of the five characters as a delimiter:
	# M (match/mismatch), S (soft clip), I (insertion), 
	# D (deletion) or H (hard clip)
	# this will return a list with the lengths of the consecutive events
	cigar_lengths = re.split('M|S|I|D|H', cigar)
	
	# try removing the empty element at the end of the list created by splitting
	# if last character is empty
	try:
		# remove it
		cigar_lengths.remove('')

	# if last character is not empty
	# a character not in [M,I,S,D,H] is in the cigar element
	# don't know what to do => skip read 
	except:
		print('\n')
		PrintTime("error", "\t\tUnexpected character found in CIGAR ("+cigar+")... Read skipped !\n")
		return False


	# split the cigar sequence by any number
	# this will return a list with the consecutive events
	cigar_chars = re.split('[0-9]+', cigar)
	
	# remove the empty element at the end of the list created by splitting
	cigar_chars.remove('')
	
	# remove N fragment and its length from cigar if 'N' in cigar:
	if "N" in cigar:
		cigar_lengths.remove(cigar_lengths[cigar_chars.index('N')])
		cigar_chars.remove('N')

	# initialize a counter to increment position only without advancing in the sequence
	cursor_pos = 0
	
	# initialize a counter to advance in the sequence without incrementing the position
	cursor_seq = 0	 

	# the need of two seperate is for cases when the event is an insertion as we need to
	# advance in the sequence of the read without incrementing the position or when the 
	# event is a deletion as we need to increment the position without advancing in the
	# sequence of the read



	# for each cigar event in the cigar events list (M|S|I|D|H)
	for i in range(0, len(cigar_chars)):
		
		# if the event is a match/mismatch
		if cigar_chars[i] == "M":
			
			# get the length of the match/mismatch event
			maxx = int(cigar_lengths[i])

			# call the specific function responsible of parsing a match/mismatch segment of the read 
			value = AddMatches(pileup, umi, chrom, start, seq, strand, cursor_pos, cursor_seq, maxx, qual, quals, MIN_BASE_QUALITY, ALL_UMIS)
			
			# get the returned values to continue from the position where the event ends
			position = value[0]
			cursor_seq = value[1]
			cursor_pos = value[2]
		
		# if the event is an insertion
		elif cigar_chars[i] == "I":
			
			# get the length of the insertion event
			maxx = int(cigar_lengths[i])


			# try to get position
			# normally, an read does not start with an insertion so the position of the insertion
			# is the last position of the previous cigar element 
			try:
				test = position

			# in some cases, the read starts with an insertion, therefore that position is not defined
			# in this case, an exception will be throw in order to attribute the start value to position
			except:
				position = start

			# call the specific function responsible of parsing an inserted segment of the read 
			value = AddInsertions(pileup, umi, chrom, position, seq, strand, cursor_pos, cursor_seq, maxx, qual, quals, MIN_BASE_QUALITY, ALL_UMIS)
			
			# get the returned values to continue from the position where the event ends
			position = value[0]
			cursor_seq = value[1]
			cursor_pos = value[2]

		# if the event is a deletion
		elif cigar_chars[i] == "D":

			# get the length of the deletion event
			maxx = int(cigar_lengths[i])

			# call the specific function responsible of parsing a deleted segment of the read 
			value = AddDeletions(pileup, umi, chrom, start, seq, strand, cursor_pos, cursor_seq, maxx, ALL_UMIS)

			# get the returned values to continue from the position where the event ends
			position = value[0]
			cursor_seq = value[1]
			cursor_pos = value[2]
		
		# if the event is a soft clip
		# soft clipping removes the clipped sequence from the read without correcting read start position
		elif cigar_chars[i] == "S":
			
			# test if the soft clip occurs at the beginning of the read
			try:
				# if not beginning of the read => cursor_pos > 0
				test = 20 / cursor_pos
				# in this case, do nothing since the soft clipping occured at the end of the read
				# since a clipping event can only occurs at the start or at the end of the read 
				continue

			# if beginning of the read => cursor_pos = 0 => exception 
			except:
				# get length of the clipped sequence
				clipped = int(cigar_lengths[i])
				# correct the read start position
				start -= clipped

				# increment the 2 cursors and continue
				cursor_pos += clipped
				cursor_seq += clipped

		# if the event is a hard clip
		# just continue (hard clipping corrects read start position and removes the clipped sequence from the read)
		elif cigar_chars[i] == "H":
			continue

