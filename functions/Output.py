
import os
import sys

from functions import * 
from GenerateVCFHeader import *
from GenerateGVCFHeader import *


def Output(full_pileup, pileup, finalVariants, INPUT, SAM, BED, FASTA, OUTPUT, MIN_BASE_QUALITY, MIN_READ_QUALITY, MIN_MAPPING_QUALITY, MIN_VARIANT_UMI, SB_METHOD, MAX_STRAND_BIAS, CORES, ALPHA, MAX_HP_LENGTH, gVCF):

	print("\n")
	PrintTime('console', "\tGenarating VCF & gVCF...") if gVCF else PrintTime('console', "\tGenarating VCF...")

	vcf = open(OUTPUT+"/"+SAM.replace(".sam", ".vcf").split("/")[-1], "w")
	vcf.write(GenerateVCFHeader(INPUT, BED, FASTA, OUTPUT, MIN_BASE_QUALITY, MIN_READ_QUALITY, MIN_MAPPING_QUALITY, MIN_VARIANT_UMI, SB_METHOD, MAX_STRAND_BIAS, CORES, ALPHA, MAX_HP_LENGTH))

	if gVCF:
		gvcf = open(OUTPUT+"/"+SAM.replace(".sam", ".gvcf").split("/")[-1], "w")
		gvcf.write(GenerateGVCFHeader(INPUT, BED, FASTA, OUTPUT, MIN_BASE_QUALITY, MIN_READ_QUALITY, MIN_MAPPING_QUALITY, MIN_VARIANT_UMI, SB_METHOD, MAX_STRAND_BIAS, CORES, ALPHA, MAX_HP_LENGTH))

	variantsFile = open(OUTPUT+"/"+SAM.replace(".sam", ".variants").split("/")[-1], "w")
	header='\t'.join(map(str,["chr","pos","amplicon","reference","idx","A","C","G","T","N","=","+","-","depth","fA","fC","fG","fT","f=","f+","f-","fref","falt","pFisherA","pFisherT","pFisherG","pFisherC","pFisherIns","pFisherDel","pFisherA_adjust","pFisherT_adjust","pFisherG_adjust","pFisherC_adjust","pFisherIns_adjust","pFisherDel_adjust","variantCall","alternative","GR","parsed","n_UMI_wt","n_UMI_mt","n_unique_UMI_wt","n_unique_UMI_mt","n_singleton_UMI_wt","n_singleton_UMI_mt","n_discordant_UMI_mut","n_concordant_UMI_mut","base_error_probability","qScore","p-value","q-value","type","SB","HP","confidence","allele_cov_plus", "allele_cov_minus", "coverage_plus", "coverage_minus"]))
	variantsFile.write(header+"\n")

	currentPos = 1.0
	lastProgress = 0.0
	totalPos = GetPileupLength(finalVariants)

	pos_to_skip = []

	# counter of n variants in VCF
	kept = 0


	variants = {}
	
	for chrom in finalVariants.keys():
		pos_to_skip = []
		for position, info in finalVariants[chrom].items():

			lastProgress = PrintProgress(currentPos, totalPos, lastProgress)

			if "DEL" in info[0]:
				if position not in pos_to_skip:
					vcf.write(info[0]+"\n")
					variantsFile.write(info[1]+"\n")
					kept += 1

					del_seq = info[0].split('\t')[3]
					del_length = len(del_seq)
					pos_to_skip.extend(range(position, position+del_length+1))
					variants[chrom.split("|")[0]+"|"+str(position)] = chrom.split("|")[1]
			else:
				vcf.write(info[0]+"\n")
				variantsFile.write(info[1]+"\n")
				kept += 1
				variants[chrom.split("|")[0]+"|"+str(position)] = chrom.split("|")[1]


			currentPos += 1


	vcf.close()
	variantsFile.close()


	if gVCF:
		# creating gVCF
		gVCF = OrderedDict()
		for chrom in full_pileup.keys():
			gVCF[chrom] = OrderedDict()
			block_counter = 1
			var_counter = 1
			del_positions = []


			notVariants = []
			var = []

			for position in full_pileup[chrom].keys():

				if chrom+"|"+str(position) not in variants.keys():
					if position not in var:
						notVariants.append(position)
				else:
					var_type = finalVariants[chrom+"|"+variants[chrom+"|"+str(position)]][position][0].split("\t")[7].split(";")[4].split("=")[1]
					if var_type == "SNV":
						var.append(position)
					elif var_type == "DEL":
						if position not in var:
							del_len = len(finalVariants[chrom+"|"+variants[chrom+"|"+str(position)]][position][0].split("\t")[3])
							for i in range(position-1, position+del_len-1):
								var.append(i)

							try:
								notVariants.remove(position-1)
							except:
								pass

					elif var_type == "INS":
						var.append(position-1)
						var.append(position)

						try:
							notVariants.remove(position-1)
						except:
							pass




			for x in notVariants:
				if 'block-'+str(block_counter) not in gVCF[chrom].keys():
					gVCF[chrom]['block-'+str(block_counter)] = [x]
				else:
					if len(gVCF[chrom]['block-'+str(block_counter)]) == 1:
						gVCF[chrom]['block-'+str(block_counter)].append(x)
					else:
						if x == gVCF[chrom]['block-'+str(block_counter)][-1]+1:
							gVCF[chrom]['block-'+str(block_counter)].append(x)
						else:
							last_x = gVCF[chrom]['block-'+str(block_counter)][-1] 
							
							if last_x+1 in var:

								gVCF[chrom]['variant-'+str(var_counter)] = [last_x+1]
								var_counter += 1


							block_counter += 1
							gVCF[chrom]['block-'+str(block_counter)] = [x]



				

		for chrom in gVCF.keys():
			for block, positions in gVCF[chrom].items():
				if 'block' in block:
					start = positions[0]
					end = positions[-1]
					qScores = []
					depths = []

					for pos in range(start, end+1):
						qScores.append(full_pileup[chrom][pos]['qScore'])
						depths.append(full_pileup[chrom][pos]['depth'])

					qScore = int(round(sum(qScores)/float(len(qScores)), 1))
					depth  = round(sum(depths) /float(len(depths)), 0)

					line = chrom+"\t"+str(start)+"\t"+"-"+"\t"+full_pileup[chrom][start]['ref']+"\t"+"<NON_REF>"+"\t"+"-"+"\t"+"-"+"\t"+"END="+str(end)+";"+"MEAN_DP="+str(depth)+";"+"MEAN_QSCORE="+str(qScore)+"\n"
					
				
				else:
					position = positions[0]
					try:
						line = finalVariants[chrom+"|"+variants[chrom+"|"+str(position)]][position][0]
					except:
						line = finalVariants[chrom+"|"+variants[chrom+"|"+str(position+1)]][position+1][0]
					
					line = line.split("\t")
					line[4] += ",<NON_REF>"
					line = "\t".join(line)+"\n"

				gvcf.write(line)


		gvcf.close()

	print("\n")
	PrintTime('console', "\tDone")


	return kept