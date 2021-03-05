
# THIS FILE CONTAINS GENERAL FUNCTIONS THAT ARE NEEDED TO 
# SUCCESSFULLY CALL THE VARIANTS IN THE SAM FILE 
#

import os
import sys
import time
import math
import pysam
import datetime
from collections import OrderedDict
from collections import defaultdict


# a class to define colors for each scenario
class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'

# a function to print status messages + the time 
# + a specific color for each status
def PrintTime(nature, message, *args):
		
		if nature == "console":
			c_start = bcolors.BOLD+bcolors.OKBLUE
			c_end = bcolors.ENDC
		elif nature == "warning":
			c_start = bcolors.BOLD+bcolors.WARNING
			c_end = bcolors.ENDC
		elif nature == "error":
			c_start = bcolors.BOLD+bcolors.FAIL
			c_end = bcolors.ENDC
		elif nature == "red":
			c_start = bcolors.BOLD+bcolors.FAIL
			c_end = bcolors.ENDC
		elif nature == "green":
			c_start = bcolors.BOLD+bcolors.OKGREEN
			c_end = bcolors.ENDC
		
		if args:
			message = message % args
			
		print("[ " + time.strftime('%X') + " ]\t"+c_start+ message + c_end)
		sys.stdout.flush()
		sys.stderr.flush()


# this function will return a list with elements that 
# only occured once
def GetSingletons(l):
	d = defaultdict(int)
	for x in l: d[x] += 1
	k = [x for x in l if d[x] == 1]
	return k

# this function will return a list with elements that 
# occured more than once
def RemoveSingletons(l):
	d = defaultdict(int)
	for x in l: d[x] += 1
	k = [x for x in l if d[x] > 1]
	return k

	
# function that calculates and displays progress efficiently
def PrintProgress(currentRead, totalReads, lastProgress,):
	progress = int((currentRead/totalReads*100))
	
	if int(progress) > int(lastProgress):
		### write progress in status file
		sys.stdout.write("\r[ " + time.strftime('%X') + " ]\t\t\tWorking "+str(int(progress))+" %")
		sys.stdout.flush()
		# time.sleep(0.001)

	return progress


# function to retieve intial RAM Usage and the exact launching time of the script 
def Start_dev(path):
	startTime = datetime.datetime.now()
	tmp_name = '.'+str(int(time.time()))
	os.mkdir(tmp_name)
	os.system('python '+path+'/RAM.py '+tmp_name+' &')
	return [tmp_name, startTime]




# function to retieve final RAM Usage and the exact completion time of the script 
# it will then calculate RAM Usage of the program and its execution time
# finally, it will displays the infos in the ocnsole and exit
def Exit_dev(tmp_name, startTime):
	
	# gen total time
	endTime = datetime.datetime.now()
	totalTime = endTime - startTime

	# get ram usage
	usedRAM = 0.0
	for x in os.listdir(tmp_name):
		if ".ram" in x:
			usedRAM = float(x[1:-4])
			break

	stop = open(str(tmp_name)+'/.done', 'w')
	stop.close()

	message = "\tElapsed time :  "+str(totalTime)
	PrintTime("console", message)

	try:
		import psutil

		message = "\tRAM USAGE    :  "+str(round(usedRAM, 2))+" GB"+"\n\n"
		PrintTime("console", message)

		time.sleep(1)
		for x in os.listdir(tmp_name):
			os.remove(tmp_name+"/"+x)
		
		os.rmdir(tmp_name)	
	
	except:
		print("\n")

	exit()


def Start(path):

	startTime = datetime.datetime.now()
	tmp_name = '.'+str(int(time.time()))
	return [tmp_name, startTime]



def Exit(tmp_name, startTime):
	
	endTime = datetime.datetime.now()
	totalTime = endTime - startTime

	message = "\t\tElapsed time :  "+str(totalTime)
	PrintTime("green", message)
	PrintTime("console", "\tEND")
	print("\n")


	exit()



# this function will go through the pileup dictionnary and replace sets with lists
# this is because to successfully dump the pileup as a msppack object, sets are not
# supported and should be replaced with basic lists
def RemoveSets(d):
	for chrom , infos in d.items():
		for position, composition in infos.items():
			d[chrom][position]['A'][2] = list(composition['A'][2])
			d[chrom][position]['C'][2] = list(composition['C'][2])
			d[chrom][position]['G'][2] = list(composition['G'][2])
			d[chrom][position]['T'][2] = list(composition['T'][2])
			
			for key, value in composition['in'].items():
				d[chrom][position]['in'][key][2] = list(value[2])
			for key, value in composition['del'].items():
				d[chrom][position]['del'][key][2] = list(value[2])

	return d



# this function will print the counts in the pileup at a specific genomic location
# it will prints the genomic location, the A,C,G,T,ins,del forward and reverse counts
# as well as the reference base at this location and the corresponding depth
def PrintPileup(pileup, index):
	index = index.split(":")
	try:
		composition = pileup[index[0]][int(index[1])]
	except:
		print("\n")
		PrintTime('error', "\tPosition "+":".join(index)+" not in the pileup")
		return None

	PrintTime('console', "\tPosition "+":".join(index)+" Breakdown :")
	print("\n")

	print(" ref : "+composition['ref']+"\n")

	print("  A  : "+str(composition['A'][0]+composition['A'][1]))
	print("  C  : "+str(composition['C'][0]+composition['C'][1]))
	print("  G  : "+str(composition['G'][0]+composition['G'][1]))
	print("  T  : "+str(composition['T'][0]+composition['T'][1]))
		
	n_ins_pos = 0
	n_ins_neg = 0
	for value in composition['in'].values():
		n_ins_pos += value[0]
		n_ins_neg += value[1]
	
	n_del_pos = 0
	n_del_neg = 0
	for value in composition['del'].values():
		n_del_pos += value[0]
		n_del_neg += value[1]

	print("  +  : "+str(n_ins_pos+n_ins_neg))
	print("  -  : "+str(n_del_pos+n_del_neg))
	print("  #  : "+str(composition['depth'])+"\n")

	print(' base_error_probability : '+str(composition['base_error_probability']))
	print("   HP   : "+str(composition['HP'])+"\n")

	if composition['alt'] == "in":
		AO = n_ins_pos+n_ins_neg
	elif composition['alt'] == "del":
		AO = n_del_pos+n_del_neg
	else:
		AO = composition[composition["alt"]][0]+composition[composition["alt"]][1]
	
	print(' q-value : '+str(round(composition['q-value'], 4)))
	print(" ALT1    : "+composition['alt'])
	print(' AF1     : '+str(round(float(AO)/composition["depth"], 4)))
	print(' n_umi1  : '+str(composition['alt_concordant']))
	print(' SB1     : '+str(round(composition['SB'], 3))+"\n")


	if composition['alt2'] != None:
		if composition['alt2'] == "in":
			AO = n_ins_pos+n_ins_neg
		elif composition['alt2'] == "del":
			AO = n_del_pos+n_del_neg
		else:
			AO = composition[composition["alt2"]][0]+composition[composition["alt2"]][1]
		
		print(' q-value2: '+str(round(composition['q-value2'], 4)))
		print(" ALT2    : "+composition['alt2'])
		print(' AF2     : '+str(round(float(AO)/composition["depth"], 4)))
		print(' n_umi2  : '+str(composition['alt_concordant2']))
		print(' SB2     : '+str(round(composition['SB2'], 3))+"\n")



	if composition['alt3'] != None:
		if composition['alt3'] == "in":
			AO = n_ins_pos+n_ins_neg
		elif composition['alt3'] == "del":
			AO = n_del_pos+n_del_neg
		else:
			AO = composition[composition["alt3"]][0]+composition[composition["alt3"]][1]
		
		print(' q-value3: '+str(round(composition['q-value3'], 4)))
		print(" ALT3    : "+composition['alt3'])
		print(' AF3     : '+str(round(float(AO)/composition["depth"], 4)))
		print(' n_umi3  : '+str(composition['alt_concordant3']))
		print(' SB3     : '+str(round(composition['SB3'], 3))+"\n")


	print("\n")



# this function will take a position as an integer and split it into
# three digits separated by a space. For example : 1245 => 1 245 and
# 125487 => 125 487. This function is used in the PrintVariant function 
# to facilitate the variant position annotation 
def FormatPosition(position):
	n = len(str(position))

	formatted = "" 
	l = 1

	while n > 0:
		a = str(position)[n-1]
		formatted = a+formatted
		if l == 3:
			l = 0
			if n > 1:
				formatted = " "+formatted

		l +=1

		n -= 1

	return formatted



# this function will print the details about a certain variant given under one of these
# three formats : chrX:1215458A>T, chrX:1215458>del or chrX:1215458A>ins.
# it will print the genomic location, the A,C,G,T,ins,del, depth, total counts, the 
# reference and alt base, AF, HP length, qNoise, SB and the number of alt concordant UMI.
def PrintVariant(pileup, variant):
	variant = variant.split(":")
	chrom = variant[0]
	ref = variant[1][variant[1].index('>')-1] if "del" not in variant[1] and "ins" not in variant[1] else ""
	alt = variant[1][variant[1].index('>')+1:]
	position = int(variant[1].replace(ref+'>'+alt, ''))

	try:
		composition = pileup[chrom][position]
	except:
		try:
			composition = pileup[chrom][position+1]
		except:
			print("\n")
			PrintTime('error', "\tVariant "+":".join(variant)+" not in the pileup")
			return None

	print("\n")
	# PrintTime('console', "\tVariant "+chrom+" : "+FormatPosition(position)+" "+ref+" > "+alt+" Breakdown :")
	PrintTime('console', "\tVariant "+":".join(variant)+" Breakdown :")
	print("\n")

	
	#############################################################
	###############   if variant is substitution ################
	#############################################################


	if len(alt) == 1:
		print(" REF : "+composition['ref'])
		print(" ALT : "+alt+"\n")

		print("  A  : "+str(composition['A'][0]+composition['A'][1]))
		print("  C  : "+str(composition['C'][0]+composition['C'][1]))
		print("  G  : "+str(composition['G'][0]+composition['G'][1]))
		print("  T  : "+str(composition['T'][0]+composition['T'][1]))
		
		n_ins_pos = 0
		n_ins_neg = 0
		for value in composition['in'].values():
			n_ins_pos += value[0]
			n_ins_neg += value[1]
		
		n_del_pos = 0
		n_del_neg = 0
		for value in composition['del'].values():
			n_del_pos += value[0]
			n_del_neg += value[1]

		print("  +  : "+str(n_ins_pos+n_ins_neg))
		print("  -  : "+str(n_del_pos+n_del_neg))
		print("  #  : "+str(composition['depth'])+"\n")

		
		if composition['alt'] == alt:
			AO = composition[composition["alt"]][0]+composition[composition["alt"]][1]
	
			print('   AF   : '+str(round(float(AO)/composition["depth"], 4)))
			print(' q-value: '+str(round(composition['q-value'], 4)))
			print('  n_umi : '+str(composition['alt_concordant']))
			print('   SB   : '+str(round(composition['SB'], 3)))

		elif composition['alt2'] == alt:
			AO = composition[composition["alt2"]][0]+composition[composition["alt2"]][1]
		
			print('   AF   : '+str(round(float(AO)/composition["depth"], 4)))
			print(' q-value: '+str(round(composition['q-value2'], 4)))
			print('  n_umi : '+str(composition['alt_concordant2']))
			print('   SB   : '+str(round(composition['SB2'], 3)))

		else:
			AO = composition[composition["alt3"]][0]+composition[composition["alt3"]][1]
		
			print('   AF   : '+str(round(float(AO)/composition["depth"], 4)))
			print(' q-value: '+str(round(composition['q-value3'], 4)))
			print('  n_umi : '+str(composition['alt_concordant3']))
			print('   SB   : '+str(round(composition['SB3'], 3)))
		
		print("   HP   : "+str(composition['HP'])+"\n")


	else:

		#############################################################
		#################   if variant is deletion ##################
		#############################################################


		if alt == "del":

			try:
				composition = pileup[chrom][position+1]
			except:
				print("\n")
				PrintTime('error', "\tVariant "+":".join(variant)+" not in the pileup")
				return None

			print(" REF : "+composition['ref'])
			print(" ALT : - \n")

			print("  A  : "+str(composition['A'][0]+composition['A'][1]))
			print("  C  : "+str(composition['C'][0]+composition['C'][1]))
			print("  G  : "+str(composition['G'][0]+composition['G'][1]))
			print("  T  : "+str(composition['T'][0]+composition['T'][1]))
			
			n_ins_pos = 0
			n_ins_neg = 0
			for value in composition['in'].values():
				n_ins_pos += value[0]
				n_ins_neg += value[1]
			
			n_del_pos = 0
			n_del_neg = 0
			for value in composition['del'].values():
				n_del_pos += value[0]
				n_del_neg += value[1]

			print("  +  : "+str(n_ins_pos+n_ins_neg))
			print("  -  : "+str(n_del_pos+n_del_neg))
			print("  #  : "+str(composition['depth'])+"\n")

			AO = n_del_pos+n_del_neg
			print('   AF   : '+str(round(float(AO)/composition["depth"], 4)))

			if composition['alt'] == alt:
				print(' q-value: '+str(round(composition['q-value'], 4)))
				print('  n_umi : '+str(composition['alt_concordant']))
				print('   SB   : '+str(round(composition['SB'], 3)))

			elif composition['alt2'] == alt:
				print(' q-value: '+str(round(composition['q-value2'], 4)))
				print('  n_umi : '+str(composition['alt_concordant2']))
				print('   SB   : '+str(round(composition['SB2'], 3)))

			else:
				print(' q-value: '+str(round(composition['q-value3'], 4)))
				print('  n_umi : '+str(composition['alt_concordant3']))
				print('   SB   : '+str(round(composition['SB3'], 3)))
			
			print("   HP   : "+str(composition['HP'])+"\n")



		else:

			#############################################################
			#################  if variant is insertion ##################
			#############################################################
			
			alt = "in"

			try:
				composition1 = pileup[chrom][position+1]
			except:
				print("\n")
				PrintTime('error', "\tVariant "+":".join(variant)+" not in the pileup")
				return None

			print(" REF : "+composition['ref'])
			print(" ALT : + \n")

			print("  A  : "+str(composition['A'][0]+composition['A'][1]))
			print("  C  : "+str(composition['C'][0]+composition['C'][1]))
			print("  G  : "+str(composition['G'][0]+composition['G'][1]))
			print("  T  : "+str(composition['T'][0]+composition['T'][1]))
			
			n_ins_pos = 0
			n_ins_neg = 0
			for value in composition1['in'].values():
				n_ins_pos += value[0]
				n_ins_neg += value[1]
			
			n_del_pos = 0
			n_del_neg = 0
			for value in composition['del'].values():
				n_del_pos += value[0]
				n_del_neg += value[1]

			print("  +  : "+str(n_ins_pos+n_ins_neg))
			print("  -  : "+str(n_del_pos+n_del_neg))
			print("  #  : "+str(composition['depth'])+"\n")

			AO = n_ins_pos+n_ins_neg
			print('   AF   : '+str(round(float(AO)/composition["depth"], 4)))

			if composition1['alt'] == alt:
				print(' q-value: '+str(round(composition1['q-value'], 4)))
				print('  n_umi : '+str(composition1['alt_concordant']))
				print('   SB   : '+str(round(composition1['SB'], 3)))

			elif composition1['alt2'] == alt:
				print(' q-value: '+str(round(composition1['q-value'], 4)))
				print('  n_umi : '+str(composition1['alt_concordant2']))
				print('   SB   : '+str(round(composition1['SB2'], 3)))

			else:
				print(' q-value: '+str(round(composition1['q-value'], 4)))
				print('  n_umi : '+str(composition1['alt_concordant3']))
				print('   SB   : '+str(round(composition1['SB3'], 3)))
			
			print("   HP   : "+str(composition1['HP'])+"\n")


	print("\n")



# this function takes the pileup dictionnary as an argument
# for each chrom and at each position, it will calculates the 
# depth. the depth = nA+ + nA- + nC+ + nC- + nG+ + nG- + nT+ + nT- + nDel+ + nDel- 
# insertion counts are not taken into account in depth calculation
def AddDepths(d):
	for chrom , infos in d.items():
		for position, composition in infos.items():
			depth =  composition['A'][0]+composition['A'][1]
			depth += composition['C'][0]+composition['C'][1]
			depth += composition['G'][0]+composition['G'][1]
			depth += composition['T'][0]+composition['T'][1]

			n_ins = 0
			for value in composition['in'].values():
				n_ins += value[0]+value[1]
			n_del = 0
			for value in composition['del'].values():
				n_del += value[0]+value[1]

			# insertions shouldn't be taken into account when calculating depth ???
			depth += n_ins
			depth += n_del

			d[chrom][position]['depth'] = depth

	return d


# this function takes the pileup dictionnary as an argument
# for each chrom and at each position, it will calculate use the
# position depth and the total base quality scores to estimate 
# the noise level at the position
def EstimateNoise(d):
	for chrom, infos in d.items():
		for position, composition in infos.items():
			# get total base quality scores at this position
			qScore = d[chrom][position]['base_error_probability']
			
			# get depth
			depth = d[chrom][position]['depth']
			
			n_del = 0
			for value in composition['del'].values():
				n_del += value[0]+value[1]

			# substract deletion from depth because deletions are not taken into account
			# in the total qscore since deletions don't have qscore
			depth -= n_del
			# divide by the depth to get mean qscore of the position
			mean_qscore = math.floor(float(qScore)/depth) if depth > 0 else 0
			d[chrom][position]['qScore'] = mean_qscore

			# qScore = -10.log(base_error_probability)
			# use this formula to calculate base_error probability
			base_error_probability = 10**(float(mean_qscore)/-10)
			# add the estimated base_error_probability to the position in the pileup
			d[chrom][position]['base_error_probability'] = round(base_error_probability, 6)

	return d



# this function takes the pileup dictionnary as an argument and the refernce genome
# for each chrom and at each position, it will retrieve the reference base in the given
# genome and add it to the pileup dictionnary
def AddRefBases(d, f):
	for chrom , infos in d.items():
		for position, composition in infos.items():
			ref = f[chrom][position-1]
			d[chrom][position]['ref'] = ref.upper()

	return d


# this function takes 3 arguments : the quality string of the read, 
# the quals chracters string and the minimum read quality
# if the mean quality of the read is < minimum read quality, the
# function will return False (meaning the read is invalid) 
def QualIsValid(qual, quals, min_read_qual):
	total = 0
	l_read = len(qual)
	threshold = min_read_qual*l_read

	step = 0
	for q in qual:
		total += quals[q]
		step += 1
		if step == 10:
			if total >= threshold:
				return True
			step=0


	return False



# this function will check if the read is valid
# a valid orphan must not be an orphan, must be aligned to chromosome
# at a certain position, must have an UMI of length 12 (customizable),
# must have a read quality > min_read_quality and a mapping quality 
# different than 255. if all these conditions are met, the read is valid
# and the function returns True
def ReadIsValid(flag, chrom, pos, umi, cigar, mapq, min_mapping_qual, qual, quals, min_read_qual):

	orphan = not bool(int("{0:b}".format(int(flag))[-2]))

	# test if orphan
	if orphan:
		return False


	try:
		# test chromosome format
		test = chrom.index('chr')

		# test pos format
		test = 25/int(pos)

		# test umi length
		# lUMI = len(umi)
		# test = 25/(lUMI-1)
		# test = 25/(lUMI-2)
		# test = "xx12xx".index(str(len(umi)))
		
		# test cigar format
		test = 25/(1-len(cigar))

		# test read quality
		if not QualIsValid(qual, quals, min_read_qual):
			return False

		# test mapping quality
		if mapq == 255 or mapq < min_mapping_qual:
			return False



		return True

	except:

		return False


# this function will calculate the total number of keys in the given dictionnary
# and return it as a value
def GetPileupLength(pileup):
	total = 0
	
	for chrom , infos in pileup.items():
		for position, composition in infos.items():
			total += 1

	return total


# this function will calculate the total number of variants in the pileup dictionnary
# and return it as a value
def GetTotalVariants(pileup):
	total = 0
	
	for chrom , infos in pileup.items():
		for position, composition in infos.items():
			total += 1
			if composition["alt2"] != None:
				total += 1
			if composition["alt3"] != None:
				total += 1

	return total


# this function will return the number of lines in a file
def GetTotalLines(filename):
	with open(filename) as f:
		for i, l in enumerate(f):
			pass
	return i + 1


# this function takes the pileup dictionnary as an argument and the refernce genome
# for each chrom and at each position, it will calculate the length of the homopolymer 
# and add it to the pileup dictionnary
def AddHomoPolymers(d, f):
	for chrom , infos in d.items():
		for position, composition in infos.items():
			ref = f[chrom][position-1]

			hp = 1
			# check hp left - max - 20
			for i in range(2, 22):
				x = f[chrom][position-i]
				if x == ref:
					hp += 1
				else:
					break

			# check hp right - max - 20
			for i in range(0, 20):
				x = f[chrom][position+i]
				if x == ref:
					hp += 1
				else:
					break


			d[chrom][position]['HP'] = hp

	return d



# this function takes the pileup dictionnary as an argument
# it will sort the dictionnary by chromosome first, and by position
# second and return the sorted dictionnary.
def SortPileup(pileup):

	num_chroms = []
	alpha_chroms = []
	new_pileup = OrderedDict()
	final_pileup = OrderedDict()

	for chrom in pileup.keys():
		try:
			num_chroms.append(int(chrom.replace('chr', '')))
		except:
			alpha_chroms.append(chrom.replace('chr', ''))

	num_chroms = sorted(num_chroms)
	alpha_chroms = sorted(alpha_chroms)
	sorted_chroms = num_chroms + alpha_chroms
	
	for chrom in sorted_chroms:
		chrom = 'chr'+str(chrom)
		poss = pileup[chrom].keys()
		sorted_poss = sorted(poss)
		new_pileup[chrom] = OrderedDict()
		for pos in sorted_poss:
			new_pileup[chrom][pos] = pileup[chrom][pos]

	return new_pileup



# this function aims to create a distinct copy of the given dict 
def CopyPileup(pileup):
	new_pileup = OrderedDict()
	for chrom in pileup.keys():
		new_pileup[chrom] = OrderedDict()
		for position in pileup[chrom].keys():
			new_pileup[chrom][position] = pileup[chrom][position]

	return new_pileup



# this function will call saltools to convert the BAM file into a SAM file
def BAMtoSAM(bam):
	sam = bam.replace('.bam', '.sam')
	samFile = open(sam, "w")
	samFile.close()
	pysam.view('-o', sam, bam, save_stdout=sam)
	return sam




# this function will give a variant a confidence level based
# on its q-value, alt_concordant_umi, strand bias and hp length.
def ComputeConfidence(qval, alt_discordant, alt_concordant, hp, Cp, Cm, Vp, Vm):

	sb = CalculateStrandBias(Cp, Cm, Vp, Vm, "default")
	confs = ['low', 'average', 'high', 'strong', 'certain']
	levels = []

	# check q-value
	if qval == 0:
		levels.append(5)
	elif qval < 0.00005:
		levels.append(4)
	elif qval < 0.0005:
		levels.append(3)
	elif qval < 0.005:
		levels.append(2)
	elif qval < 0.01:
		levels.append(1)
	else:
		levels.append(0)



	# check disordant/concordant ratio
	if alt_discordant == 0:
		levels.append(5)
	else:
		if alt_concordant > 30:
			levels.append(5)
		elif alt_concordant > 20:
			levels.append(4)
		elif alt_concordant > 15:
			levels.append(3)
		elif alt_concordant > 10:
			levels.append(2)
		else:
			levels.append(1)


	# check SB
	if sb < 0.5:
		levels.append(5)
	elif sb < 0.60:
		levels.append(4)
	elif sb < 0.70:
		levels.append(3)
	elif sb < 0.80:
		levels.append(2)
	elif sb < 0.90:
		levels.append(1)
	else:
		levels.append(0)

	
	# check HP
	if hp < 3:
		levels.append(5)
	elif hp == 3:
		levels.append(4)
	elif hp == 4:
		levels.append(3)
	elif hp == 5:
		levels.append(2)
	elif hp == 6:
		levels.append(1)
	else:
		levels.append(0)

	conf = int(math.floor(float(sum(levels))/len(levels)))
	return [conf, confs[conf-1]]



# this function will calculate the SB of the variant depending on the method the user selected
def CalculateStrandBias(Cp, Cm, Vp, Vm, method):
	if method == "default":
		return abs(((float(Vp)/(Cp+Vp))-(float(Vm)/(Cm+Vm)))/((float(Vp+Vm))/(Cp+Vp+Cm+Vm)))
	else:
		return float(max([Vp*Cm,Vm*Cp]) / (Vp*Cm + Vm*Cp))



def PrintProgramName():
	print(
"""
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
||||||||| ||||||| |||  ||||||  |||       ||||||||||||  |||||||  |||       ||||       |||||    ||||       ||| ||||||||||||||||||
||||||||| ||||||| ||| | |||| | |||||| |||||||||||||||  |||||||  ||| ||||| |||| ||||| |||| |||| ||| ||||| ||| ||||||||||||||||||
||||||||| ||||||| ||| || || || |||||| |||||||||||||||  |||||||  ||| ||||| |||| ||||| ||| ||||||||| ||||| ||| ||||||||||||||||||
||||||||| ||||||| ||| |||  ||| |||||| ||||||      ||||  |||||  ||||       ||||       ||| |||||||||       ||| ||||||||||||||||||
||||||||| ||||||| ||| |||||||| |||||| |||||||||||||||||  |||  ||||| ||||| |||| ||  ||||| ||||||||| ||||| ||| ||||||||||||||||||
|||||||||  |||||  ||| |||||||| |||||| ||||||||||||||||||  |  |||||| ||||| |||| |||  ||||| |||| ||| ||||| ||| ||||||||||||||||||
||||||||||       |||| |||||||| |||       ||||||||||||||||   ||||||| ||||| |||| ||||  |||||    |||| ||||| |||       ||||||||||||
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||V.S||||
|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
""")






# this function will return a dictionnary with all UMIS and their indexes
def GetUMIS(SAM):

	# value = set()

	# for line in open(SAM):

	# 	line = line.split('\t')
		
	# 	try:
	# 		value.add(line[0].split('_')[-1])

	# 	except:
	# 		continue

	
	# c = 0
	# final = {}

	# for el in value:
	# 	final[el] = c
	# 	c+=1

	# return final
	return 1