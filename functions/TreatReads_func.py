
from functions import *

def TreatReads(SAM, BED, MIN_BASE_QUALITY, MIN_READ_QUALITY, MIN_MAPPING_QUALITY, ALL_UMIS):

	quals_str = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"
	quals = {}


	i=0
	for q in quals_str:
		quals[quals_str[i]] = i
		i += 1


	pileup = ParseBED(BED)


	currentLine = 1.0
	lastProgress = 0.0
	totalLines = GetTotalLines(SAM)


	
	for line in open(SAM):

		lastProgress = PrintProgress(currentLine, totalLines, lastProgress)

		currentLine += 1


		line = line.split('\t')
		
		try:
			umi = line[0].split('_')[-1]
		except:
			continue


		flag = int(line[1])
		strand = int("{0:b}".format(int(flag))[-5])
		firstInPair = bool(int("{0:b}".format(int(flag))[-7]))
		chrom = line[2]
		pos = line[3]
		mapq = int(line[4])
		cigar = line[5]
		seq = line[9]
		qual = line[10]

		if ReadIsValid(flag, chrom, pos, umi, cigar, mapq, MIN_MAPPING_QUALITY, qual, quals, MIN_READ_QUALITY):
			
			AddRead(pileup, umi, strand, chrom, int(pos), cigar, seq, qual, quals, MIN_BASE_QUALITY, ALL_UMIS)




	return pileup
