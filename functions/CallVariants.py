#
# THIS FUNCTION WILL APPLY MULTIPLE FILTERS TO THE VARINATS THAT PASS THE POISSON TEST IN ORDER TO REDUCE FALSE POSITIVES 
# 
# INPUT : 
#         -PILEUP              : (DICT)   THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
#		  -F                   : (DICT)   A DICTIONNARY CONTAINING THE REFERENCE BASES AT ALL POSITIONS OF THE GENOME 
# 		  -SB_METHOD           : (STR)    DEFAULT METHOD for SB CALCULATION OR TORRENT SUITE METHOD
# 		  -MAX_STRAND_BIAS     : (FLOAT)  THRESHOLD FOR A VARIANT TO BE CONSIDERED AS STRAND BIASED
# 		  -MIN_VARIANT_UMI     : (INT)    THRESHOLD for A VARIANT WITH A CERTAIN UMI COUNT TO BE CALLED
# 		  -MAX_HP_LENGTH       : (INT)    HOMOPOLYMER REGION LENGTH THRESHOLD FOR A VARIANT IN IT TO BE CALLED 
#
# VALUE : -FINALVARIANTS       : (DICT)   A DICTIONNARY CONTAINING ONLY THE FINAL VARIANTS THAT SUCCESSFULLY PASSED ALL FILTERS
#

from func import *
from collections import OrderedDict


def CallVariants(pileup, f, SB_METHOD, MAX_STRAND_BIAS, MIN_VARIANT_UMI, MAX_HP_LENGTH):

	# create a dictionnary to contain only final variants
	finalVariants = OrderedDict()

	print("\n")
	PrintTime('console', "\tCalling Variants...")

	# define counters to calculate progress
	currentPos = 1.0
	lastProgress = 0.0
	totalPos = float(GetPileupLength(pileup))

	# loop through the PILEUP dictionnary
	# for each chromosome
	for chrom , infos in pileup.items():
		# try-except block for first variants
		# if chromosome already in finalVariants - first variant - dict
		try:
			# do nothing
			test = finalVariants[chrom+"|v"]
		
		# if chromosome not in finalVariants - first variant - dict
		except:
			# add the chromosome to the finalVariants - first variant - dict
			finalVariants[chrom+"|v"] = OrderedDict()


		# try-except block for second variants
		# if chromosome already in finalVariants - second variant - dict
		try:
			# do nothing
			test = finalVariants[chrom+"|v2"]
		
		# if chromosome not in finalVariants - second variant - dict
		except:
			# add the chromosome to the finalVariants - second variant - dict
			finalVariants[chrom+"|v2"] = OrderedDict()


		# try-except block for third variants
		# if chromosome already in finalVariants - third variant - dict
		try:
			# do nothing
			test = finalVariants[chrom+"|v3"]
		
		# if chromosome not in finalVariants - third variant - dict
		except:
			# add the chromosome to the finalVariants - third variant - dict
			finalVariants[chrom+"|v3"] = OrderedDict()

		# for each position | composition = counts + infos (ref_allele, alt_allele, q-value, ...)
		for position, composition in infos.items():

			# function that calculates and displays progress
			lastProgress = PrintProgress(currentPos, totalPos, lastProgress)

			# retrieve list of reference UMIs from PILEUP dictionnary
			total_ref_umi = composition[composition['ref']][2]
			ref_umi_len = len(total_ref_umi)
			ref_umi = list(set(composition[composition['ref']][2]))


			# make list to contain all the possible alternative alleles
			alt_list = [composition['alt']]
			# another list to contain the indexes ["", "2", "3"]
			alt_indexes = {composition['alt']: ""}

			if composition['alt2'] != None:
				alt_list.append(composition['alt2'])
				alt_indexes[composition['alt2']] = "2"
			if composition['alt3'] != None:
				alt_list.append(composition['alt3'])
				alt_indexes[composition['alt3']] = "3"



			# for each possible alternative allele 
			for alt_element in alt_list:

				# get corresponding index
				alt_index = alt_indexes[alt_element]


				# retrieve list of variant UMIs from PILEUP dictionnary
				# if variant is insertion or deletion
				if alt_element == 'in' or alt_element == 'del':

					
					umi_list = []

					# for each insertion or deletion
					for indel, infos in composition[alt_element].items():
						# append the UMI to the umi list
						umi_list += infos[2]

				# else if variant is substitution
				else:
					# retrieve list of variant UMIs from PILEUP dictionnary
					umi_list = composition[alt_element][2]

				

				# remove variant UMIs that occurs only once and keep unique UMIs only
				alt_umi = list(set(RemoveSingletons(umi_list)))

				

				# create the list of unique noise UMIs
				# noise UMIs are the UMIs found on reads that had neither the reference base nor the variant
				alphabet = ["A", "C", "G", "T", "in", "del"]
				alphabet.remove(composition['ref'])
				alphabet.remove(alt_element)
				noise_umi = set()
				for c in alphabet:
					if c == 'in' or c == 'del':
						for indel, infos in composition[c].items():
							noise_umi.update(infos[2])
					else:
						noise_umi.update(composition[c][2])



				#calculate total + coverage
				pos_covTotal = composition['A'][0]+composition['C'][0]+composition['G'][0]+composition['T'][0]
				
				#calculate total - coverage
				neg_covTotal = composition['A'][1]+composition['C'][1]+composition['G'][1]+composition['T'][1]

				# add insertion + and - coverage to total + and - coverage
				for indel, infos in composition['in'].items():
					pos_covTotal += infos[0]
					neg_covTotal += infos[1]

				# add deletion + and - coverage to total + and - coverage
				for indel, infos in composition['del'].items():
					pos_covTotal += infos[0]
					neg_covTotal += infos[1]

				# retrieve variant + and - coverage
				# if variant is insertion or deletion => go through all insertions / deletions 
				if alt_element == 'in' or alt_element == 'del':
					pos_covAllele = 0
					neg_covAllele = 0
					for indel, infos in composition[alt_element].items():
						pos_covAllele += infos[0]
						neg_covAllele += infos[1]

				# else if variant is substitution, retrieve from PILEUP dictionnary
				else:
					pos_covAllele = composition[alt_element][0]
					neg_covAllele = composition[alt_element][1]

				

				#Strand bias computation ( we add 1 to avoid 0 values)
				Vp=float(pos_covAllele+1)
				Vm=float(neg_covAllele+1)
				
				Cp=float(pos_covTotal+1)
				Cm=float(neg_covTotal+1)
				
				SB = CalculateStrandBias(Cp, Cm, Vp, Vm, SB_METHOD)



				# ref discordant UMIs are UMIs that are found in the reference UMIs list
				# AND in the (variant|noise) UMIs list
				# ref concordant UMIs are UMIs that are found in the reference UMIs list ONLY
				ref_discordant = 0
				for umi in ref_umi:
					if umi in alt_umi or umi in noise_umi:
						ref_discordant += 1

				ref_concordant = len(ref_umi) - ref_discordant

				# alt discordant UMIs are UMIs that are found in the variant UMIs list
				# AND in the (reference|noise) UMIs list
				# alt concordant UMIs are UMIs that are found in the variant UMIs list ONLY
				alt_discordant = 0
				for umi in alt_umi:
					if umi in ref_umi or umi in noise_umi:
						alt_discordant += 1

				alt_concordant = len(alt_umi) - alt_discordant



				###########################################################################
				##################   FOR TESTING PURPOSES - START   #######################
				###########################################################################

				# print ref_umi
				# print alt_umi
				# print noise_umi
				# print "posCovAllele : "+str(pos_covAllele)
				# print "negCovAllele : "+str(neg_covAllele)
				# print "posCovTotal : "+str(pos_covTotal)
				# print "negCovTotal : "+str(neg_covTotal)
				# print "Strand bias : "+str(SB)
				# print "ref_concordant : "+str(ref_concordant)
				# print "ref_discordant : "+str(ref_discordant)
				# print "alt_concordant : "+str(alt_concordant)
				# print "alt_discordant : "+str(alt_discordant)


				# if chrom == "chr1" and position == 27106356:
				# 	print "\n"
				# 	print alt_umi
				# 	print len(alt_umi)
				# 	print alt_discordant
				# 	print alt_concordant
				# 	print SB
				# 	exit()

				###########################################################################
				###################   FOR TESTING PURPOSES - END   ########################
				###########################################################################

				# add alt_concordant and SB informations to pileup
				pileup[chrom][position]['alt_concordant'+alt_index] = alt_concordant
				pileup[chrom][position]['alt_discordant'+alt_index] = alt_discordant
				pileup[chrom][position]['SB'+alt_index] = SB

				# Allelic frequency = AF 
				AF = composition['VAF'+alt_index]
				# DePth = DP
				DP = composition['depth']
				# Allele Observations = AO
				AO = int(round(float(AF)*float(DP), 0))
				# HomoPolymer length = HP
				HP = composition['HP']
				# get error_base_probability
				base_error_probability = composition['base_error_probability']
				# compute confidence level of the variant
				conf = ComputeConfidence(composition['q-value'+alt_index], alt_discordant, alt_concordant, HP, Cp, Cm, Vp, Vm)
				CONF = conf[1]
				CONF_EXTRA = conf[1]+" ("+str(conf[0])+"/5)"




				# if the potential variant at this position passes 4 final filters => true variant
				# filter 1 : number of alt_concordant UMIs >= MIN_VARIANT_UMI value 
				# filter 2 : comptuted strand bias must be <= MAX_STRAND_BIAS value
				# filter 3 : homopolymer length must be <= MAX_HP_LENGTH value
				# filter 4 : total number of umis must be > to concordant UMIS / VAF


				if SB <= MAX_STRAND_BIAS and HP <= MAX_HP_LENGTH and alt_concordant >= MIN_VARIANT_UMI and len(list(list(alt_umi)+list(noise_umi)))+ref_umi_len >= (float(composition['alt_concordant'+alt_index])/composition['VAF'+alt_index]):

					# if the variant is not an indel
					if alt_element != "in" and alt_element != "del":

						# build index
						index = chrom+":"+str(position)
						# get its position
						positionInVCF = str(position)

						# build the variant name
						variantName = index+composition['ref']+">"+alt_element
						# give it the type SNV
						TYPE = 'SNV'

					# else if the variant is an insertion
					elif alt_element == "in":

						# build index (position in index is the position -1)
						index = chrom+":"+str(position-1)
						# gets its position (position in VCF is the position -1)
						positionInVCF = str(position-1)

						# choosing the most frequent insertion
						chosen_ins = ['', 0, 0]
						for ins, infos in pileup[chrom][position]['in'].items():
							if (infos[0]+infos[1]) > (chosen_ins[1]+chosen_ins[2]):
								chosen_ins = [ins, infos[0], infos[1]]

						# build the variant name
						variantName = index+f[chrom][position-1-1]+">"+f[chrom][position-1-1]+chosen_ins[0]
						# give it the type INS
						TYPE = 'INS'

					# else if the variant is a deletion
					else:
						# build index (position in index is the position -1)
						index = chrom+":"+str(position-1)
						# gets its position (position in VCF is the position -1)
						positionInVCF = str(position-1)

						# choosing the most frequent deletion
						chosen_del = ['', 0, 0]
						for dele, infos in pileup[chrom][position]['del'].items():
							if (infos[0]+infos[1]) > (chosen_del[1]+chosen_del[2]):
								chosen_del = [dele, infos[0], infos[1]]


						# get the reference base from the reference dictionnary at the position-1
						lastBasePos = f[chrom][position-1-1].upper()
						# get the deleted seq by going from position -> position+length 
						# of deleted sequence in the reference dictionnary
						dele = int(chosen_del[0])
						del_seq = ""
						while dele > 0:
							del_seq += f[chrom][position+len(del_seq)-1].upper()
							dele -= 1

						# build the variant name
						variantName = index+lastBasePos+del_seq+">"+lastBasePos
						# giv it the type DEL
						TYPE = 'DEL'

					# build the REF and ALT columns of the VCF file for each 
					# variant type (SNV | INS | DEL)
					if TYPE == "SNV":
						REF = composition['ref']
						ALT = alt_element
					elif TYPE == "DEL":
						REF = lastBasePos+del_seq
						ALT = lastBasePos
					else:
						REF = f[chrom][position-1-1]
						ALT = f[chrom][position-1-1]+chosen_ins[0]

				
					# create the INFO line | delimiter is ';'
					INFO = "AF="+str(AF)+";AO="+str(AO)+";DP="+str(DP)+";HP="+str(HP)+";TYPE="+str(TYPE)+";CONF="+CONF.upper()
					


					# calculate total insertion observations on both strands
					n_ins_pos = 0
					n_ins_neg = 0
					for value in composition['in'].values():
						n_ins_pos += value[0]
						n_ins_neg += value[1]
					n_ins = n_ins_pos+n_ins_neg
					
					# calculate total deletion observations on both strands
					n_del_pos = 0
					n_del_neg = 0
					for value in composition['del'].values():
						n_del_pos += value[0]
						n_del_neg += value[1]
					n_del = n_del_pos+n_del_neg

					# build the VCF line for each variant | delimiter = "\t"
					lineVCF = "\t".join([chrom, positionInVCF, variantName, REF, ALT, str(composition['q-value'+alt_index]), "PASS", INFO])
					
					# build the VARIANTS file line for each variant | delimiter = "\t"
					lineExtra = "\t".join([
						chrom, 
						positionInVCF, 
						"-", 
						REF, 
						index.replace(':', '-'), 
						str(composition['A'][0]+composition['A'][1]), 
						str(composition['C'][0]+composition['C'][1]), 
						str(composition['G'][0]+composition['G'][1]), 
						str(composition['T'][0]+composition['T'][1]), 
						"0", 
						"0", 
						str(n_ins), 
						str(n_del), 
						str(DP), 
						str((composition['A'][0]+composition['A'][1])/float(DP)),
						str((composition['C'][0]+composition['C'][1])/float(DP)),
						str((composition['G'][0]+composition['G'][1])/float(DP)),
						str((composition['T'][0]+composition['T'][1])/float(DP)),
						"0",
						str(n_ins/float(DP)),
						str(n_del/float(DP)),
						str((composition[composition['ref']][0]+composition[composition['ref']][1])/float(DP)),
						str(AO/float(DP)),
						"-", "-", "-", "-", "-", "-",
						"-", "-", "-", "-", "-", "-",
						"TRUE",
						alt_element,
						variantName.replace(':', 'g.'),
						"FALSE",
						str(ref_umi_len),
						str(len(umi_list)),
						str(len(ref_umi)), 
						str(len(set(umi_list))), 
						str(len(GetSingletons(total_ref_umi))),
						str(len(GetSingletons(umi_list))),
						str(alt_discordant),
						str(alt_concordant),
						str(base_error_probability),
						str(composition['qScore']),
						str(composition['p-value'+alt_index]),
						str(composition['q-value'+alt_index]),
						TYPE,
						str(SB),
						str(HP),
						CONF_EXTRA,
						str(pos_covAllele),
						str(neg_covAllele),
						str(pos_covTotal),
						str(neg_covTotal),
						])
					
					# add both lines to the finalvariants dictionnary
					finalVariants[chrom+"|v"+alt_index][position] = [lineVCF, lineExtra] 

			# increment to track progress
			currentPos += 1



	print("\n")
	PrintTime('console', "\tDone")

	# return finalVariants dictionnary
	return finalVariants