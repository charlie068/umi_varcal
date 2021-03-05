import os
import sys
import time
import math
import msgpack
from pyfasta import Fasta
from collections import OrderedDict


# import local modules
from functions import *
from PreprocessReads import *


def Call(config, FUNC_PATH):


	# load and define variables from the config
	INPUT                = config['input']
	BED                  = config['bed']
	FASTA                = config['fasta']
	MIN_BASE_QUALITY     = int(config['min_base_quality'])
	MIN_MAPPING_QUALITY  = int(config['min_mapping_quality'])
	MIN_READ_QUALITY     = int(config['min_read_quality'])
	MIN_VARIANT_UMI      = int(config['min_variant_umi'])
	STRAND_BIAS_METHOD   = str(config['strand_bias_method'])
	MAX_STRAND_BIAS      = float(config['max_strand_bias'])
	PILEUP               = config['pileup']
	REBUILD              = False if os.path.isfile(PILEUP) else True
	OUTPUT               = config['output']
	CORES                = config['cores']
	DEFAULT_CORES        = config['default_cores']
	ALPHA                = float(config['alpha'])
	MAX_HP_LENGTH        = int(config['max_hp_length'])
	gVCF                 = config['gvcf']
	KEEP_PILEUP			 = config['keep_pileup']


	# print parameters in the console
	PrintTime("green", "\t\tINPUT file   : "+INPUT)
	PrintTime("green", "\t\tBED file     : "+BED)
	PrintTime("green", "\t\tFASTA file   : "+FASTA)

	if PILEUP != "None":	
		PrintTime("green", "\t\tPILEUP file  : "+PILEUP)
	PrintTime("green", "\t\tOutput       : "+OUTPUT)

	PrintTime("green", "\t\tmin_base_quality      : "+str(MIN_BASE_QUALITY))
	PrintTime("green", "\t\tmin_read_quality      : "+str(MIN_READ_QUALITY))
	PrintTime("green", "\t\tmin_mapping_quality   : "+str(MIN_MAPPING_QUALITY))
	PrintTime("green", "\t\tmin_variant_umi       : "+str(MIN_VARIANT_UMI))
	PrintTime("green", "\t\tstrand_bias_method    : "+str(STRAND_BIAS_METHOD))
	PrintTime("green", "\t\tmax_strand_bias       : "+str(MAX_STRAND_BIAS))
	PrintTime("green", "\t\tmax_hp_length         : "+str(MAX_HP_LENGTH))
	PrintTime("green", "\t\talpha                 : "+str(ALPHA))
	if gVCF:
		PrintTime("green", "\t\tgVCF                  : "+str(gVCF)+" (Experimental)")
	else:
		PrintTime("green", "\t\tgVCF                  : "+str(gVCF))

	if DEFAULT_CORES:
		PrintTime("green", "\t\tcores                 : "+str(CORES)+" (default)")
	else:
		PrintTime("green", "\t\tcores                 : "+str(CORES))

	PrintTime("green", "\t\tkeep_pileup           : "+str(KEEP_PILEUP))

	PrintTime("console", "\tDone\n")





	# load the reference genome file
	f = Fasta(FASTA)




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


	# if a pileup is not given, the pileup has to be build
	if REBUILD:

		print("\n")
		PrintTime('console', "\tBuilding Pileup...")
	 	

		# if multiple cores are used
		# wait until all processes are finished == until all files are generated
		if CORES > 1:

			###############################################################################
			###########################                         ###########################
			########################### PARALLELIZED CODE START ###########################
			###########################                         ###########################
			###############################################################################

			# preprocess reads
			# if more then one core is to be used, separate the input into subfiles
			subFiles = PreprocessReads(SAM, totalLines, CORES)

			# build the empty pileup
			pileup = ParseBED(BED)

			# build the command for the instances to be launched simultanously 
			command = ""
			for subFile in subFiles:
				base = "python3 "+FUNC_PATH+"/TreatReads.py "+subFile+" "+BED+" "+str(MIN_BASE_QUALITY)+" "+str(MIN_READ_QUALITY)+" "+str(MIN_MAPPING_QUALITY)+" "+OUTPUT+" & "
				command += base
				pileupFile = subFile.replace('.sam', '.pileup')


			command = command[:-2]

			# execute the command
			os.system(command)



			# make sure to wait for all cores to finish building their pileups
			# before merging them
			finished = False
			while not finished:

				finished = True
				pileups = []

				for subFile in subFiles:
					samName = subFile.split("/")[-1]
					pileupFile = OUTPUT+"/"+samName.replace('.sam', '.pileup')
					try:
						p = msgpack.unpack(open(pileupFile, 'rb'), encoding="utf-8")
						pileups.append(p)
					except:
						finished = False
						time.sleep(1)



			# remove intermediate sam file
			for subFile in subFiles:
				os.remove(subFile)


			# merge sub pileups to obtain whole pileup
			pileup = MergeSubPileups(pileup, pileups, subFiles, OUTPUT)


			###############################################################################
			############################                       ############################
			############################ PARALLELIZED CODE END ############################
			############################                       ############################
			###############################################################################

		else:

			# get all UMI list
			ALL_UMIS = GetUMIS(SAM)

			# if only one core is to used, launch the function from here since no need to merge
			pileup = TreatReads(SAM, BED, MIN_BASE_QUALITY, MIN_READ_QUALITY, MIN_MAPPING_QUALITY, ALL_UMIS)




		# add depth to pileup
		pileup = AddDepths(pileup)

		# add variant error noise ate each position
		pileup = EstimateNoise(pileup)

		# add refernce bases in the dictionnary
		pileup = AddRefBases(pileup, f)

		# add homopolymers infos
		pileup = AddHomoPolymers(pileup, f)


		# rebuild to SAM original name
		SAM = SAM.replace('_reordered.sam', ".sam")

		# dump pileup in msgpack object
		if KEEP_PILEUP:
			with open(OUTPUT+"/"+SAM.replace(".sam", ".pileup").split("/")[-1], 'wb') as handle:
				msgpack.pack(pileup, handle, encoding="utf-8")



		print("\n")
		PrintTime('console', "\tDone")



	else:

		print("\n")
		PrintTime('console', "\tLoading Pileup...")
		
		# load pileup from msgpack object
		with open(PILEUP, 'rb') as handle:
			pileup = msgpack.unpack(handle, encoding="utf-8")
			pileup = SortPileup(pileup)
			
		PrintTime('console', "\tDone")




	# print(pileup)
	full_pileup = CopyPileup(pileup)


	### Poisson modeling to filter positions
	result = FilterPositions(pileup, ALPHA)
	pileup = result[0]
	potential = result[1]
	foundCandidates = result[2]


	if foundCandidates:
		### call final variants
		finalVariants = CallVariants(pileup, f, STRAND_BIAS_METHOD, MAX_STRAND_BIAS, MIN_VARIANT_UMI, MAX_HP_LENGTH)

		### Writing results to VCF
		final = Output(full_pileup, pileup, finalVariants, INPUT, SAM, BED, FASTA, OUTPUT, MIN_BASE_QUALITY, MIN_READ_QUALITY, MIN_MAPPING_QUALITY, MIN_VARIANT_UMI, STRAND_BIAS_METHOD, MAX_STRAND_BIAS, CORES, ALPHA, MAX_HP_LENGTH, gVCF)

		# calculate and display stats
		CalculateStats(pileup, potential, final)
	else:
		print("\n")
		message = "No candidate positions were found !\n"
		PrintTime('error', "\t"+message)

		PrintTime('console', "\tCalculating statistics...")
		
		# print out stats to console
		message = "Candidate Positions: 0"
		PrintTime('green', "\t\t"+message)
		message = "Final Variants: 0"
		PrintTime('green', "\t\t"+message)

