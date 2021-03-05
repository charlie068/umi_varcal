
#!/usr/bin/python
import os
import sys

# add functions folders the the system path
# sys.path.append(os.getcwd()+'/functions')

# get true path
FUNC_PATH = sys.argv[0].replace('umi-varcal.py', 'functions')
sys.path.append(FUNC_PATH)

# import local modules
# from ParseBED import *
# from PreprocessReads import *
# from Config import *
# from MergeSubPileups import *
# from func import *
# from PileupToCSV import *
# from FilterVariants import *
# from CallFinalVariants import *
# from Output import *


from functions import *


config = Config()


value = Start(FUNC_PATH)
tmp_name = value[0]  
startTime = value[1]  



if 'umi_length' in config:
	from Extract import *
	Extract(config)
else:
	from Call import *
	Call(config, FUNC_PATH) 




Exit(tmp_name, startTime)

