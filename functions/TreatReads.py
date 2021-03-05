import os
import sys
import time
import math
import msgpack
import multiprocessing
from pyfasta import Fasta


# import local modules
from functions import *

# get all UMI list
ALL_UMIS = GetUMIS(sys.argv[1])

# launch the function that analyzes a set of reads
# this function will create a pileup from these reads
pileup = TreatReads(sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), ALL_UMIS)

# dump pileup in msgpack object
with open(sys.argv[6]+"/"+sys.argv[1].replace(".sam", ".pileup").split("/")[-1], 'wb') as handle:
	msgpack.pack(pileup, handle)