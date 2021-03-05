
# THIS FUNCTION WILL APPLY THE POISSON TEST ON EACH POSITION TO TRY AND FIND A SIGNIFCANT VARIANTS. THEN, P-VALUES ARE CORRECTED
# WITH THE BENJAMINI-HOCHBERG METHOD TO OBTAIN Q-VALUES. ONLY POSITIONS WITH Q-VALUE > ALPHA ARE KEPT
#
# INPUT : 
#         -PILEUP              : (DICT)   THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
# 		  -ALPHA               : (FLOAT)  TYPE 1 ERROR PROBABLITY OR ALPHA LEVEL
#
# VALUE : -PILEUP              : (DICT)   THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC

import math
import operator
from scipy.stats import poisson
import statsmodels.stats.multitest as smm
from func import *


def FilterPositions(pileup, ALPHA):
	print("\n")
	PrintTime('console', "\tSearching for candidate positions...")

	# define counters to calculate progress
	currentPos = 1.0
	lastProgress = 0.0
	totalPos = GetPileupLength(pileup)


	# create two dicts to contain p-values and q-values
	pValues = {}
	qValues = {}

	# create an empty array to stock positions that are not variants
	toRemove = []

	# loop through pileup
	# for each chromosome
	for chrom, infos in pileup.items():
		# for each position | composition = counts
		for position, composition in infos.items():

			# get estimated error probability for this position
			base_error_probability = composition['base_error_probability']

			# function that calculates and displays progress
			lastProgress = PrintProgress(currentPos, totalPos, lastProgress)

			# calculate the estimaded lambda value for the poisson model
			estimated_lambda=int(float(base_error_probability)*composition["depth"])

			# get insertions + and - counts
			n_ins_pos = 0
			n_ins_neg = 0
			for value in composition['in'].values():
				n_ins_pos += value[0]
				n_ins_neg += value[1]
			
			# get deletions + and - counts
			n_del_pos = 0
			n_del_neg = 0
			for value in composition['del'].values():
				n_del_pos += value[0]
				n_del_neg += value[1]
			
			# create a dictionnary with total counts for only base / ins / del counts for this position 
			calls = { 
				"A": composition["A"][0]+composition["A"][1],
				"C": composition["C"][0]+composition["C"][1],
				"G": composition["G"][0]+composition["G"][1],
				"T": composition["T"][0]+composition["T"][1],
				"in": n_ins_pos+n_ins_neg,
				"del": n_del_pos+n_del_neg
			}

			# remove the reference base counts in the dictionnary
			calls.pop(pileup[chrom][position]["ref"], None)	

			# order the remaining alleles by value
			calls = OrderedDict(sorted(calls.items(), key=operator.itemgetter(1), reverse=True))
			
			# check if there is a possible second alternative allele
			pileup[chrom][position]["alt"] =list(calls.keys())[0] if list(calls.values())[0] > 0 else None
			# check if there is a possible second alternative allele
			pileup[chrom][position]["alt2"]=list(calls.keys())[1] if list(calls.values())[1] > 0 else None
			# check if there is a possible third alternative allele
			pileup[chrom][position]["alt3"]=list(calls.keys())[2] if list(calls.values())[2] > 0 else None


			if composition["depth"] <= 0 or pileup[chrom][position]["alt"] == None:
				# increment counter for progress
				currentPos += 1.0
				# if depth == 0 => remove position from pileup (means position is not covered by BAM/SAM file)
				# if the alternative allele count = 0 => no variants at this position
				toRemove.append(chrom+":"+str(position))
				continue

			
			# get the alternative allele count
			n_alt = calls[pileup[chrom][position]["alt"]]

			# calculate allelic frequency of the variant and add it to the PILEUP dictionnary	
			pileup[chrom][position]["VAF"]=float(n_alt)/composition["depth"]
			
			# calculate the pvalue for the variant to be background noise and add it to the PILEUP dictionnary
			pileup[chrom][position]["p-value"] = float(1-poisson.cdf(n_alt, estimated_lambda))
			
			# add the index and the p value to the pvalues dict 
			pValues[chrom+":"+str(position)+"-"] = pileup[chrom][position]["p-value"]

			if pileup[chrom][position]["alt2"] != None:
				# get the second alternative allele count
				n_alt = calls[pileup[chrom][position]["alt2"]]

				# calculate allelic frequency of the variant and add it to the PILEUP dictionnary	
				pileup[chrom][position]["VAF2"]=float(n_alt)/composition["depth"]
				
				# calculate the pvalue for the variant to be background noise and add it to the PILEUP dictionnary
				pileup[chrom][position]["p-value2"] = float(1-poisson.cdf(n_alt, estimated_lambda))
				
				# append the p value to the pvalues list 
				# pValues.append(pileup[chrom][position]["p-value2"])
				pValues[chrom+":"+str(position)+"-2"] = pileup[chrom][position]["p-value2"]



			if pileup[chrom][position]["alt3"] != None:
				# get the third alternative allele count
				n_alt = calls[pileup[chrom][position]["alt3"]]

				# calculate allelic frequency of the variant and add it to the PILEUP dictionnary	
				pileup[chrom][position]["VAF3"]=float(n_alt)/composition["depth"]
				
				# calculate the pvalue for the variant to be background noise and add it to the PILEUP dictionnary
				pileup[chrom][position]["p-value3"] = float(1-poisson.cdf(n_alt, estimated_lambda))
				
				# append the p value to the pvalues list 
				# pValues.append(pileup[chrom][position]["p-value3"])
				pValues[chrom+":"+str(position)+"-3"] = pileup[chrom][position]["p-value3"]



			currentPos += 1.0






	# remove unwanted positions
	for index in toRemove:
		chrom = index.split(":")[0]
		position = int(index.split(":")[1])

		pileup[chrom].pop(position, None)




	# reset toRemove array
	toRemove = []

	# recalculate pileup length
	totalPos = GetPileupLength(pileup)

	# apply FDR Correction (Benjamini/Hochberg Correction) to correct pvalues 
	# based on multiple test theory
	try:
		corrected = smm.multipletests(list(pValues.values()), alpha=0.05, method='fdr_bh')
		
		# retrieve results of the correction
		# retieve results just in case we would need it
		results = corrected[0]
		qValues_list = corrected[1]
		foundCandidates = True
	except:
		results = []
		qValues_list = []
		foundCandidates = False




	# retrieve qValues
	counter = 0
	for index in pValues.keys():
		qValues[index] = qValues_list[counter]
		counter += 1


	# loop the the qValues dict and the each qValues
	# to its corresponding position in the pileup
	for index, qvalue in qValues.items():
		index = index.split(":")
		chrom = index[0]
		position = int(index[1].split("-")[0])
		alt_index = index[1].split("-")[1]

		pileup[chrom][position]['q-value'+alt_index] = qvalue





	currentPos = 1.0

	# loop through pileup
	# for each chromosome
	for chrom, infos in pileup.items():
		# for each position | composition = counts		
		for position, composition in infos.items():
			# function to calculate and display progress
			lastProgress = PrintProgress(currentPos, totalPos, lastProgress)

			# if qvalue >= alpha(default = 0.05) => consider as artifact
			if pileup[chrom][position]["q-value"] >= ALPHA:
				# remove position from pileup dictionnary
				toRemove.append(chrom+":"+str(position))
			else:
				# if there is a possible second allele variant
				if composition["alt2"] != None:
					# if qvalue2 >= alpha(default = 0.05) => consider as artifact
					if pileup[chrom][position]["q-value2"] >= ALPHA:
						pileup[chrom][position]["alt2"] = None


				# if there is a possible third allele variant
				if composition["alt3"] != None:
					# if qvalue2 >= alpha(default = 0.05) => consider as artifact
					if pileup[chrom][position]["q-value3"] >= ALPHA:
						pileup[chrom][position]["alt3"] = None

			# increment counter for progress
			currentPos += 1




	# remove unwanted positions
	for index in toRemove:
		chrom = index.split(":")[0]
		position = int(index.split(":")[1])

		pileup[chrom].pop(position, None)



	# get total kept positions
	potential = GetTotalVariants(pileup)


	print("\n")
	PrintTime('console', "\tDone")


	# return PILEUP dictionnary
	return [pileup, potential, foundCandidates]