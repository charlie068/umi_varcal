
import os
import sys


from functions import *



def Extract(config):

	# retrieve the configuration parameters
	INPUT = config['input']
	FASTA = config['fasta']
	umi_length = int(config['umi_length'])
	bwa_threads = int(config['bwa_threads'])


	# print the configuration inthe console
	PrintTime("green", "\t\tINPUT file   : "+INPUT)
	PrintTime("green", "\t\tFASTA file   : "+FASTA)
	PrintTime("green", "\t\tUMI length   : "+str(umi_length))
	PrintTime("green", "\t\tBWA Threads  : "+str(bwa_threads))

	PrintTime("console", "\tDone\n")


	# if input is bam => launch samtools view command
	# to convert it to sam
	if ".bam" in INPUT and ".sam" not in INPUT:

		print("\n")
		PrintTime('console', "\tConverting BAM to SAM...")

		SAM = BAMtoSAM(INPUT)

		PrintTime('console', "\tDone")

	else:
		# else => sam = input
		SAM = INPUT




	# get total number of reads
	totalLines = GetTotalLines(SAM)




	# build the output directories and files
	samName = SAM.split('/')[-1].replace('.sam', '')
	raw_sam = 'EXTRACTED/'+samName+"/"+samName+'_extracted.sam'
	raw_bam = 'EXTRACTED/'+samName+"/"+samName+'_extracted.bam'
	sorted_raw_bam = 'EXTRACTED/'+samName+"/"+samName+'_extracted_sorted.bam'
	try:
		os.mkdir('EXTRACTED')
	except:
		pass

	try:
		os.mkdir('EXTRACTED/'+samName)
	except:
		pass
	FASTQ_1 = 'EXTRACTED/'+samName+"/"+samName+'_R1.fastq'
	FASTQ_2 = 'EXTRACTED/'+samName+"/"+samName+'_R2.fastq'

	output1 = open(FASTQ_1, 'w')
	output2 = open(FASTQ_2, 'w')



	# create an empty dictionnary to contain read ids that were already seen
	seen = {}


	nLine = 1.0
	lastProgress = 0.0

	PrintTime('console', "\tReverting to FASTQ...")

	# loop through the sam file
	for line in open(SAM):

		# print progress
		lastProgress = PrintProgress(nLine, totalLines, lastProgress)

		# if line is not in header
		if "@SQ" not in line and "SN:chr" not in line and "LN:" not in line:

			line = line.split('\t')

			# if line corresponds to a read
			if len(line) > 10:

				# get read if and add @ to create fastq read id
				tmp_readID = "@"+line[0]
				
				# extract umi from sequence
				tmp_seq = line[9]
				umi = tmp_seq[:umi_length]
				seq = tmp_seq[umi_length:]
				
				# extract umi quality from quality string
				tmp_qual = line[10]
				qual = tmp_qual[umi_length:]

				# add umi to the end of the read id
				readID = tmp_readID+"_"+str(umi)

				# build the fastq read
				fq_line = readID+"\n"+seq+"\n+\n"+qual+"\n"
				
				# if the read id is already in seen dict
				try:
					# retrieve the value in seen as Read 1
					fq_line1 = seen[readID]
					# write in R1 file the value in seen
					output1.write(fq_line1)
					# write in R2 file the value in the current line
					output2.write(fq_line)
					# delete the entry from seen
					del seen[readID] 

				except:
					# if the read id is not in seen dict
					# add it to the dict
					seen[readID] = fq_line

		# increment the ine counter to track progress
		nLine += 1

	# close the fastq r1 and fastq r2 files
	output1.close()
	output2.close()


	print('\r\n')
	PrintTime('console', "\tDone")
	print('\n\r')


	PrintTime('console', "\tAligning FASTQ to reference genome...")
	# launch bwa to do the alignment : FASTQ_R1 + FASTQ_R2 + FASTA ==> SAM
	os.system("bwa mem -M -t "+str(bwa_threads)+" "+FASTA+" "+FASTQ_1+" "+FASTQ_2+" > "+raw_sam+" 2> /dev/null")
	PrintTime('console', "\tDone")
	print('\n\r')

	PrintTime('console', "\tConverting to BAM...")
	# convert to BAM with samtools
	os.system("samtools view -S -b "+raw_sam+" > "+raw_bam+" 2> /dev/null")
	PrintTime('console', "\tDone")
	print('\n\r')

	PrintTime('console', "\tSorting BAM...")
	# sort bam with samtools 
	os.system("samtools sort "+raw_bam+" > "+sorted_raw_bam+" 2> /dev/null")
	PrintTime('console', "\tDone")
	# remove unsorted bam after sorting
	os.remove(raw_bam)
	os.rename(sorted_raw_bam, raw_bam)
	print('\n\r')

	PrintTime('console', "\tIndexing BAM...")
	# indexing bam with samtools 
	os.system("samtools index "+raw_bam+" 2> /dev/null")
	PrintTime('console', "\tDone")