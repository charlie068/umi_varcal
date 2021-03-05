
# THIS FUNCTION ALLOWS TO CONFIGURE THE PARAMETERS FROM THE COMMAND LINE
# AND PROPERLY PASS THEM TO THE CORRESPONDING FUNCTION CALLS
#
# INPUT : 
# 		  NONE
#
# VALUE : -CONFIG              : (DICT)  A DICTIONNARY THAT CONTAINS ALL NECESSARY THE PARAMETERS AND THEIR VALUES
#		

import os
import sys
import time

### import functions from the script in the same dir
from func import *
from PrintHelp import *


def Config():
	
	config={}


	if len(sys.argv) == 1 or "--help" in sys.argv or "-h" in sys.argv:

		if "extract" in sys.argv:
			PrintHelp("extract")
		elif "call" in sys.argv:
			PrintHelp("call")
		else:
			PrintHelp("general")



	# configure the extraction tool
	elif sys.argv[1] == "extract":
		if len(sys.argv) == 2:
			PrintHelp("extract")


		### just in case script path contains '-' or '--'
		sys.argv[0] = "None"
		sys.argv.remove('extract')


		# default parameters
		config["bwa_threads"] = 1


		### required parameters
		# input       : path to the bam/sam file
		# fasta       : path to the fasta file
		# umi length  : length of the umi

		required = ['--input | -i', '--fasta | -f', '--umi_length | -l']
		for param in required:
			missing = False
			for arg in param.split(" | "):
				try:
					test = sys.argv.index(arg)
					missing = False
					break
				except:
					missing = True



			if missing:
				PrintTime('error', "\tThe parameter "+param+" is required but missing !\n\t\t\tExiting...")
				exit()
		


		print("\n")
		PrintTime("console", "\tConfiguring Parameters...")

		# differentiate the '-' in fasta of the '-' for arguments
		try:
			fastaIdx = sys.argv.index("--fasta")

			fastaIdx+=1
			sys.argv[fastaIdx] = sys.argv[fastaIdx].replace("-", ":")
		except:
			try:
				fastaIdx = sys.argv.index("-f")

				fastaIdx+=1
				sys.argv[fastaIdx] = sys.argv[fastaIdx].replace("-", ":")
			except:
				pass



		# differentiate the '-' in input of the '-' for arguments
		try:
			inputIdx = sys.argv.index("--input")

			inputIdx+=1
			sys.argv[inputIdx] = sys.argv[inputIdx].replace("-", ":")
		except:
			try:
				inputIdx = sys.argv.index("-i")

				inputIdx+=1
				sys.argv[inputIdx] = sys.argv[inputIdx].replace("-", ":")
			except:
				pass










		### verify that no arguments are empty
		params_in_config = ['input', 'fasta', 'umi_length']
		params_mini = ['i', 'f', 'l']
		params = ['--input | -i', '--fasta | -f', '--umi_length | -l']

		pointer = 1
		while pointer < len(sys.argv):
			param = sys.argv[pointer]
			try:
				value = sys.argv[pointer+1]
			except:
				if "--" in param:
					if param in ['--input', '--fasta', '--umi_length']:
						PrintTime('error', "\tThe parameter "+params[params_in_config.index(param.replace("--", ""))]+" is required and cannot be empty !\n\t\t\tExiting...")
						exit()
					else:
						if param.replace("--", "") in params_in_config:
							PrintTime('error', "\tThe parameter's "+params[params_in_config.index(param.replace("--", ""))]+" value cannot be empty !\n\t\t\tExiting...")
							exit()
						else:
							PrintTime('error', "\tThe parameter "+param+" is unknown !\n\t\t\tExiting...")
							exit()
				elif "-" in param and "--" not in param:
					if param in ['-i', '-f', '-l']:
						PrintTime('error', "\tThe parameter "+params[params_mini.index(param.replace("-", ""))]+" is required and cannot be empty !\n\t\t\tExiting...")
						exit()
					else:
						if param.replace("-", "") in params_mini:
							PrintTime('error', "\tThe parameter's "+params[params_mini.index(param.replace("-", ""))]+" value cannot be empty !\n\t\t\tExiting...")
							exit()
						else:
							PrintTime('error', "\tThe parameter "+param+" is unknown !\n\t\t\tExiting...")
							exit()
				else:
					PrintTime('error', "\tThe parameter "+param+" is unknown !\n\t\t\tExiting...")
					exit()


			if "-" in value and len(value) == 2:
				PrintTime('error', "\tThe parameter "+param+" is required and cannot be empty !\n\t\t\tExiting...")
				exit()

			if "--" in value and len(value) > 3:
				PrintTime('error', "\tThe parameter "+param+" is required and cannot be empty !\n\t\t\tExiting...")
				exit()

			pointer += 2





		args=""
		for arg in sys.argv:
			args+=" "+arg


		args = args.replace("--", "|")
		args = args.replace("-", "|")

		args = args.split("|")
		del args[0]


		for arg in args:
			param = arg.split(" ")[0]
			value = arg.split(" ")[1]

			if param == "input" or param == "i":
				param = "input"
				if value[0] != "/":
					value = os.getcwd()+"/"+value.replace(":", "-")
				config[param]=value.replace(":", "-")

			if param == "fasta" or param == "f":
				param = "fasta"
				if value[0] != "/":
					value = os.getcwd()+"/"+value.replace(":", "-")
				config[param]=value.replace(":", "-")


			if param == "umi_length" or param == "l":
				param = "umi_length"
				try:
					config[param]=int(value)
				except:
					PrintTime('error', "\tThe parameter's --"+param+" value should be an integer !\n\t\t\tExiting...")
					exit()


			if param == "bwa_threads" or param == "t":
				param = "bwa_threads"
				try:
					config[param]=int(value)
				except:
					PrintTime('error', "\tThe parameter's --"+param+" value should be an integer !\n\t\t\tExiting...")
					exit()
			


		return config





	# configure the variant calling tool
	elif sys.argv[1] == "call":

		if len(sys.argv) < 3:

			PrintHelp("call")

		
		### just in case script path contains '-' or '--'
		sys.argv[0] = "None"
		sys.argv.remove('call')


		### default values for parameters
		config["min_base_quality"]    = 10
		config["min_read_quality"]    = 20
		config["min_mapping_quality"] = 20
		config["min_variant_umi"]     = 5
		config["strand_bias_method"]  = "default"
		config["output"]              = os.getcwd()
		config["pileup"]              = "None"
		config["cores"]               = 1
		config["default_cores"]       = True
		config["max_hp_length"]       = 7
		config["alpha"]               = 0.05
		config["gvcf"]                = False
		config["keep_pileup"]         = True



		### required parameters
		# input : path to the bam/sam file
		# fasta  : path to the fasta file
		# bed  : path to the bed file

		required = ['--input | -i', '--fasta | -f', '--bed | -b']
		for param in required:
			missing = False
			for arg in param.split(" | "):
				try:
					test = sys.argv.index(arg)
					missing = False
					break
				except:
					missing = True



			if missing:
				PrintTime('error', "\tThe parameter "+param+" is required but missing !\n\t\t\tExiting...")
				exit()
		

		# print program name
		PrintProgramName()
		PrintTime("console", "\tConfiguring Parameters...")

		# differentiate the '-' in fasta of the '-' for arguments
		try:
			fastaIdx = sys.argv.index("--fasta")

			fastaIdx+=1
			sys.argv[fastaIdx] = sys.argv[fastaIdx].replace("-", ":")
		except:
			try:
				fastaIdx = sys.argv.index("-f")

				fastaIdx+=1
				sys.argv[fastaIdx] = sys.argv[fastaIdx].replace("-", ":")
			except:
				pass




		# differentiate the '-' in bed of the '-' for arguments
		try:
			bedIdx = sys.argv.index("--bed")

			bedIdx+=1
			sys.argv[bedIdx] = sys.argv[bedIdx].replace("-", ":")
		except:
			try:
				bedIdx = sys.argv.index("-b")

				bedIdx+=1
				sys.argv[bedIdx] = sys.argv[bedIdx].replace("-", ":")
			except:
				pass
		



		# differentiate the '-' in input of the '-' for arguments
		try:
			inputIdx = sys.argv.index("--input")

			inputIdx+=1
			sys.argv[inputIdx] = sys.argv[inputIdx].replace("-", ":")
		except:
			try:
				inputIdx = sys.argv.index("-i")

				inputIdx+=1
				sys.argv[inputIdx] = sys.argv[inputIdx].replace("-", ":")
			except:
				pass





		# differentiate the '-' in output of the '-' for arguments
		try:
			outputIdx = sys.argv.index("--output")
		
			outputIdx+=1
			sys.argv[outputIdx] = sys.argv[outputIdx].replace("-", ":")
		except:
			try:
				outputIdx = sys.argv.index("-o")

				outputIdx+=1
				sys.argv[outputIdx] = sys.argv[outputIdx].replace("-", ":")
			except:
				pass



		# differentiate the '-' in pileup of the '-' for arguments
		try:
			pileupIdx = sys.argv.index("--pileup")

			pileupIdx+=1
			sys.argv[pileupIdx] = sys.argv[pileupIdx].replace("-", ":")
		except:
			try:
				pileupIdx = sys.argv.index("-p")

				pileupIdx+=1
				sys.argv[pileupIdx] = sys.argv[pileupIdx].replace("-", ":")
			except:
				pass







		### verify that no arguments are empty		
		params_in_config = ['input', 'bed', 'fasta', 'min_base_quality', 'min_read_quality', 'min_mapping_quality', 'min_variant_umi', 'strand_bias_method', 'max_strand_bias', 'output', 'pileup', 'cores', 'alpha', 'max_hp_length', 'gvcf', 'keep_pileup']
		params_mini = ['i', 'b', 'f', 'x', 'x', 'x', 'x', 'x', 'x', 'o', 'p', 'c', 'x', 'x', 'x', 'x']
		params = ['--input | -i', '--bed | -b', '--fasta | -f', '--min_base_quality', '--min_read_quality', '--min_mapping_quality', '--min_variant_umi', '--strand_bias_method', '--max_strand_bias', '--output | -o', '--pileup | -p', '--cores | -c', '--alpha' , '--max_hp_length', '--gvcf', '--keep_pileup']

		pointer = 1
		while pointer < len(sys.argv):
			param = sys.argv[pointer]
			try:
				value = sys.argv[pointer+1]
			except:
				if "--" in param:
					if param in ['--input', '--bed', '--fasta']:
						PrintTime('error', "\tThe parameter "+params[params_in_config.index(param.replace("--", ""))]+" is required and cannot be empty !\n\t\t\tExiting...")
						exit()
					else:
						if param.replace("--", "") in params_in_config:
							PrintTime('error', "\tThe parameter's "+params[params_in_config.index(param.replace("--", ""))]+" value cannot be empty !\n\t\t\tExiting...")
							exit()
						else:
							PrintTime('error', "\tThe parameter "+param+" is unknown !\n\t\t\tExiting...")
							exit()
				elif "-" in param and "--" not in param:
					if param in ['-i', '-b', '-f']:
						PrintTime('error', "\tThe parameter "+params[params_mini.index(param.replace("-", ""))]+" is required and cannot be empty !\n\t\t\tExiting...")
						exit()
					else:
						if param.replace("-", "") in params_mini:
							PrintTime('error', "\tThe parameter's "+params[params_mini.index(param.replace("-", ""))]+" value cannot be empty !\n\t\t\tExiting...")
							exit()
						else:
							PrintTime('error', "\tThe parameter "+param+" is unknown !\n\t\t\tExiting...")
							exit()
				else:
					PrintTime('error', "\tThe parameter "+param+" is unknown !\n\t\t\tExiting...")
					exit()


			if "-" in value and len(value) == 2:
				PrintTime('error', "\tThe parameter "+param+" is required and cannot be empty !\n\t\t\tExiting...")
				exit()

			if "--" in value and len(value) > 3:
				PrintTime('error', "\tThe parameter "+param+" is required and cannot be empty !\n\t\t\tExiting...")
				exit()

			pointer += 2





		args=""
		for arg in sys.argv:
			args+=" "+arg


		args = args.replace("--", "|")
		args = args.replace("-", "|")

		args = args.split("|")
		del args[0]

		sb_set = False

		for arg in args:
			param = arg.split(" ")[0]
			value = arg.split(" ")[1]

			if param == "input" or param == "i":
				param = "input"
				if value[0] != "/":
					value = os.getcwd()+"/"+value.replace(":", "-")
				config[param]=value.replace(":", "-")

			if param == "fasta" or param == "f":
				param = "fasta"
				if value[0] != "/":
					value = os.getcwd()+"/"+value.replace(":", "-")
				config[param]=value.replace(":", "-")

			if param == "bed" or param == "b":
				param = "bed"
				if value[0] != "/":
					value = os.getcwd()+"/"+value.replace(":", "-")
				config[param]=value.replace(":", "-")


			if param == "pileup" or param == "p":
				param = "pileup"
				if value[0] != "/":
					value = os.getcwd()+"/"+value.replace(":", "-")
				config[param]=value.replace(":", "-")


			if param == "output" or param == "o":
				param = "output"
				if value[0] != "/":
					value = os.getcwd()+"/"+value.replace(":", "-")
				config[param]=value.replace(":", "-")


			if param == "cores" or param == "c":
				param = "cores"
				try:
					config[param]=int(value)
					config["default_cores"] = False
				except:
					PrintTime('error', "\tThe parameter's --"+param+" value should be an integer !\n\t\t\tExiting...")
					exit()


			if param == "min_base_quality":
				try:
					config[param]=float(value)
				except:
					PrintTime('error', "\tThe parameter's --"+param+" value should be a float or an integer !\n\t\t\tExiting...")
					exit()

			if param == "min_read_quality":
				try:
					config[param]=float(value)
				except:
					PrintTime('error', "\tThe parameter's --"+param+" value should be a float or an integer !\n\t\t\tExiting...")
					exit()

			if param == "min_mapping_quality":
				try:
					config[param]=float(value)
				except:
					PrintTime('error', "\tThe parameter's --"+param+" value should be a float or an integer !\n\t\t\tExiting...")
					exit()

			if param == "min_variant_umi":
				try:
					config[param]=float(value)
				except:
					PrintTime('error', "\tThe parameter's --"+param+" value should be a float or an integer !\n\t\t\tExiting...")
					exit()

			if param == "strand_bias_method":
				value = str(value).replace("'", '').replace('"', "")
				try:
					test = value.index('default')
					config[param]=value
				except:
					try:
						test = value.index('torrent_suite')
						config[param]=value
					except:
						PrintTime('error', "\tThe parameter --"+param+" can only be set to \"default\" or \"torrent_suite\" !\n\t\t\tExiting...")
						exit()

			if param == "max_strand_bias":
				try:
					config[param]=float(value)
					sb_set = True
				except:
					PrintTime('error', "\tThe parameter's --"+param+" value should be a float or an integer !\n\t\t\tExiting...")
					exit()	

			if param == "alpha":
				try:
					config[param]=float(value)
				except:
					PrintTime('error', "\tThe parameter's --"+param+" value should be a float or an integer !\n\t\t\tExiting...")
					exit()

			if param == "max_hp_length":
				try:
					config[param]=float(value)
				except:
					PrintTime('error', "\tThe parameter's --"+param+" value should be a float or an integer !\n\t\t\tExiting...")
					exit()

			if param == "gvcf":

				value = value.lower()
				if value in ['true', 'false']:
					value = True if value == 'true' else False
					config[param]=value
				else:
					PrintTime('error', "\tThe parameter's --"+param+" value should be a boolean (True/False) !\n\t\t\tExiting...")
					exit()

			if param == "keep_pileup":

				value = value.lower()
				if value in ['true', 'false']:
					value = True if value == 'true' else False
					config[param]=value
				else:
					PrintTime('error', "\tThe parameter's --"+param+" value should be a boolean (True/False) !\n\t\t\tExiting...")
					exit()



		if not sb_set:
			if config["strand_bias_method"] == "default":
				config["max_strand_bias"] = 1.0
			else:
				config["max_strand_bias"] = 0.743

		return config



	else:
		PrintTime('error', "\tPlease precise the tool you want to use first!\n\t\t\tExiting...")
		exit()
