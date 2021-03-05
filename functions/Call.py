import os
import sys
import time
import math
import msgpack
from pyfasta import Fasta
from collections import OrderedDict
import glob
from subprocess import Popen, call
from timeit import default_timer as timer

# import local modules
from functions import *
from PreprocessReads import *


def Call(config, FUNC_PATH):

    start = timer()
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
    KEEP_PILEUP             = config['keep_pileup']


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
    #print(type(f))
    #print(len(f))







    # get total number of reads

    bedpath=os.path.dirname(BED)
    #print('bedpath: {}'.format(bedpath))
    bedfiles=glob.glob(bedpath+'/*.bed')
    bedfiles.sort()
    #print(bedfiles)
    PrintTime('console', "PreprocessReads")
    
    bednb=0
    for bedfile in bedfiles:
        bednb+=1
    # if a pileup is not given, the pileup has to be build
        if REBUILD:
            
                #split bed in chromosome
            #list of bed file per chromosome
            li_chr=[]
            with open(bedfile,'r+') as fbed:
                for line in fbed:
                    li_chr.append(line.split('\t')[0])
            li_chr=list(set(li_chr))
            li_chr=' '.join(li_chr)
            #split bam:
            #print(INPUT)
            
            if ".bam" in INPUT and ".sam" not in INPUT:
                print("\n")
                PrintTime('console', "\tConverting BAM to SAM...")
                SAM=INPUT.replace('.bam',"_"+str(bednb)+".sam")
                command="samtools view -@"+str(CORES)+" "+INPUT+" "+li_chr+" -o "+SAM
                #print(command)
                #sys.exit() 
                call(command, shell=True)
                
                # if input is bam => launch samtools view command
                # to convert it to sam
            
        
                print("\n")
                PrintTime('console', "\tConverting BAM to SAM...")
                #command=
                #SAM = BAMtoSAM(bamoutname)
                #SAM= bamoutname
        
                PrintTime('console', "\tDone")
            else:
                PrintTime('console', "No indexed BAM provided, Please provide an indexed bam file")
                sys.exit()

            totalLines = GetTotalLines(SAM)
            subFiles = PreprocessReads(SAM, totalLines, CORES)
            
            print("\n")
            PrintTime('console', "\tBuilding Pileup...")
            PrintTime('console', "File {} out of {}".format(bednb,len(bedfiles)))

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
                

                # build the empty pileup
                pileup = ParseBED(bedfile)

                # build the command for the instances to be launched simultanously 
                command = []
                for subFile in subFiles:
                    command.append("python3 "+FUNC_PATH+"/TreatReads.py "+subFile+" "+bedfile+" "+str(MIN_BASE_QUALITY)+" "+str(MIN_READ_QUALITY)+" "+str(MIN_MAPPING_QUALITY)+" "+OUTPUT)
                    #command += base
                #pileupFile = subFile.replace('.sam', '.pileup')
                
                #base_list=base.split('&')
                #pipe=[]
                #for bas in base_list:
                #    pipe.append(subprocess.call(bas, shell=True))
                #print (command)
                print("\n")
                PrintTime('console', "TreatReads")
                procs = [ Popen(i, shell=True) for i in command ]
                for p in procs:
                   p.wait()
                procs=None



                # remove intermediate sam file
                #for subFile in subFiles:
            #        os.remove(subFile)

                PrintTime('console', "MergeSubPileups")
                # merge sub pileups to obtain whole pileup
                pileup = MergeSubPileups(pileup,subFiles, OUTPUT)


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


            #try:
            PrintTime('console', "AddDepths")
            # add depth to pileup
            pileup = AddDepths(pileup)
            PrintTime('console', "Estimate Noise")
            # add variant error noise ate each position
            pileup = EstimateNoise(pileup)
            PrintTime('console', "AddRefBases")
            # add refernce bases in the dictionnary
            pileup = AddRefBases(pileup, f)
            PrintTime('console', "homopolymers")
            # add homopolymers infos
            pileup = AddHomoPolymers(pileup, f)

            PrintTime('console', "SAMreplace")
            # rebuild to SAM original name
            SAM = SAM.replace('_reordered.sam', ".sam")
            #except:
        #        print('problem with file {}'.format(bednb))
            
            # dump pileup in msgpack object
            if KEEP_PILEUP:
                PrintTime('console', "dump")
                with open(OUTPUT+"/"+SAM.replace(".sam", "_"+str(bednb)+".pileup").split("/")[-1], 'wb') as handle:
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



        PrintTime('console', "CopyPileup")
        # print(pileup)
        full_pileup = CopyPileup(pileup)

        PrintTime('console', "Poisson")
        ### Poisson modeling to filter positions
        #result = FilterPositions(pileup, ALPHA)
        #pileup = result[0]
        #potential = result[1]
        #foundCandidates = result[2]
        pileup, potential, foundCandidates = FilterPositions(pileup, ALPHA)
        #result=None

        filestats=OUTPUT+"/"+os.path.basename(INPUT).replace(".bam", ".stats")
        if foundCandidates:
            PrintTime('console', "finalVariants")
            ### call final variants

            finalVariants = CallVariants(pileup, f, STRAND_BIAS_METHOD, MAX_STRAND_BIAS, MIN_VARIANT_UMI, MAX_HP_LENGTH)
            
            PrintTime('console', "Writing results")
            ### Writing results to VCF
            #namefilesout=SAM.replace(".sam","")
            #namefilesout=namefilesout+'_'+str(bednb)+'.sam'
            final = Output(full_pileup, pileup, finalVariants, INPUT, SAM, BED, FASTA, OUTPUT, MIN_BASE_QUALITY, MIN_READ_QUALITY, MIN_MAPPING_QUALITY, MIN_VARIANT_UMI, STRAND_BIAS_METHOD, MAX_STRAND_BIAS, CORES, ALPHA, MAX_HP_LENGTH, gVCF)
            PrintTime('console', "calculate stats")
            # calculate and display stats
            coverage, avg_depth, uniformity, potential, final=CalculateStats(pileup, potential, final)
            PrintTime('console', "Writing stats file")
            with open(filestats,'a+') as fileout:
                fileout.write('file {}, BAM/BED coverage: {}, Average design depth: {}, Uniformity: {}, Candidate Positions: {}, Final Variants: {} \n'.format(bednb, coverage, avg_depth, uniformity, potential, final))
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
            
            PrintTime('console', "Writing stats file")
            with open(filestats,'a+') as fileout:
                fileout.write('file {}, No candidate positions were found !\n'.format(bednb))
                
        pileup=None
        full_pileup=None
        potential=None
        foundCandidates=None
        finalVariants=None
        # remove intermediate sam file
        for subFile in subFiles:
            os.remove(subFile)
        os.remove(SAM)
    end = timer()
    etime=time.strftime("%H:%M:%S", time.gmtime(end - start))
    with open(filestats,'a+') as fileout:
                fileout.write('Elapse time:   {}\n'.format(etime))
    PrintTime('console', "Elapse time:   {}".format(etime))

