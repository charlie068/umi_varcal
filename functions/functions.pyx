import os
import sys
import re
import time
import math
import pysam
import msgpack
from pyfasta import Fasta
from collections import OrderedDict
from collections import defaultdict
import datetime
from func import PrintTime, PrintProgress
import operator
from scipy.stats import poisson
import statsmodels.stats.multitest as smm
import psutil





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

#
# THIS FUNCTION WILL APPLY MULTIPLE FILTERS TO THE VARINATS THAT PASS THE POISSON TEST IN ORDER TO REDUCE FALSE POSITIVES 
# 
# INPUT : 
#         -PILEUP              : (DICT)   THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
#          -F                   : (DICT)   A DICTIONNARY CONTAINING THE REFERENCE BASES AT ALL POSITIONS OF THE GENOME 
#           -SB_METHOD           : (STR)    DEFAULT METHOD for SB CALCULATION OR TORRENT SUITE METHOD
#           -MAX_STRAND_BIAS     : (FLOAT)  THRESHOLD FOR A VARIANT TO BE CONSIDERED AS STRAND BIASED
#           -MIN_VARIANT_UMI     : (INT)    THRESHOLD for A VARIANT WITH A CERTAIN UMI COUNT TO BE CALLED
#           -MAX_HP_LENGTH       : (INT)    HOMOPOLYMER REGION LENGTH THRESHOLD FOR A VARIANT IN IT TO BE CALLED 
#
# VALUE : -FINALVARIANTS       : (DICT)   A DICTIONNARY CONTAINING ONLY THE FINAL VARIANTS THAT SUCCESSFULLY PASSED ALL FILTERS
#


def CallVariants(pileup, f, SB_METHOD, MAX_STRAND_BIAS, MIN_VARIANT_UMI, MAX_HP_LENGTH):

    # create a dictionnary to contain only final variants
    finalVariants = {}

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
            finalVariants[chrom+"|v"] = {}


        # try-except block for second variants
        # if chromosome already in finalVariants - second variant - dict
        try:
            # do nothing
            test = finalVariants[chrom+"|v2"]
        
        # if chromosome not in finalVariants - second variant - dict
        except:
            # add the chromosome to the finalVariants - second variant - dict
            finalVariants[chrom+"|v2"] = {}


        # try-except block for third variants
        # if chromosome already in finalVariants - third variant - dict
        try:
            # do nothing
            test = finalVariants[chrom+"|v3"]
        
        # if chromosome not in finalVariants - third variant - dict
        except:
            # add the chromosome to the finalVariants - third variant - dict
            finalVariants[chrom+"|v3"] = {}

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
                #     print "\n"
                #     print alt_umi
                #     print len(alt_umi)
                #     print alt_discordant
                #     print alt_concordant
                #     print SB
                #     exit()

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
    pileup=None
    f=None
    # return finalVariants dictionnary
    return finalVariants


# THIS FUNCTION ALLOWS TO PARSE A CIGAR SEGMENT CONTAINING ONLY DELETED BASES AND INCREMENT
# A SPECIFIC PILEUP DICTIONNARY ACCORDINGLY
#
# INPUT : -PILEUP            : (DICT) THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
#           -UMI               : (STR)  UMI SEQUENCE OF THE READ
#           -CHROM             : (STR)  THE CHROMOSOME MAPPED TO THE READ
#           -START             : (INT)  THE START POSITION OF THE READ
#           -SEQ               : (STR)  THE SEQUENCE OF THE READ
#           -STRAND            : (INT)  THE STRAND OF THE READ (0 = FORWARD | 1 = REVERSE)
#           -CURSOR_POS        : (INT)  THE CURRENT POSITION IN THE READ
#           -CURSOR_SEQ        : (INT)  THE CURRENT POSITION IN THE SEQUENCE OF THE READ
#           -MAXX              : (INT)  THE LENGTH OF THE CIGAR 'D' ELEMENT
#          -ALL_UMIS          : (DICT) A DICTIONNARY CONTAINING UMIS INDEXES
#
# VALUE : -POSITION          : (INT)  THE FINAL POSITION THAT'S BEEN PARSED IN THE READ
#          -CURSOR_POS        : (INT)  THE CURRENT POSITION IN THE READ
#           -CURSOR_SEQ        : (INT)  THE CURRENT POSITION IN THE SEQUENCE OF THE READ
#

def AddDeletions(pileup, umi, chrom, start, seq, strand, cursor_pos, cursor_seq, maxx, ALL_UMIS):

    # create a cursor to move in the deleted sequence 
    del_cursor = 0

    # for each position between limits
    for position in range(start+cursor_pos, start+cursor_pos+maxx):
        
        # try because base position in read could not be in BED file
        # therefore, the position could not be in the PILEUP dictionnary
        # if position not in BED, skip the read
        try:
            # test if chrom and position are in the PILEUP dictionnary
            test = pileup[chrom][position]
            
            # if this deletion was already seen at this location
            # maxx-del_cursor corresponds to the length of the deleted sequence
            # this number will serve as an entry in the PILEUP dictionnary for
            # deleted sequences
            try:
                # increment the corresponding counter in the PILEUP dictionnary
                pileup[chrom][position]['del'][maxx-del_cursor][strand] += 1
                # add the umi to the corresponding set specific to the deleted sequence
                pileup[chrom][position]['del'][maxx-del_cursor][2].append(umi)

            # if this is the first time we see this deleted sequence at this position, 
            # an entry should be created with its appropriated structure (2 counters + set)  
            except:
                # create the entry in PILEUP dictionnary
                pileup[chrom][position]['del'][maxx-del_cursor] = [0, 0, [umi]]
                # increment the corresponding counter in the PILEUP dictionnary
                pileup[chrom][position]['del'][maxx-del_cursor][strand] += 1


            
            # advance in the deleted sequence if the deletion length > 1
            del_cursor += 1


        # if position not in BED, skip the read
        except:
            pass    


    # increment the cursor_pos with the length of the cigar element when done adding it
    cursor_pos += maxx

    # return position, cursor_seq and cursor_pos to continue from where we left off
    return [position, cursor_seq, cursor_pos]




# THIS FUNCTION ALLOWS TO PARSE A CIGAR SEGMENT CONTAINING ONLY DELETED BASES AND INCREMENT
# A SPECIFIC PILEUP DICTIONNARY ACCORDINGLY
#
# INPUT : -PILEUP            : (DICT) THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
#           -UMI               : (STR)  UMI SEQUENCE OF THE READ
#           -STRAND            : (INT)  THE STRAND OF THE READ (0 = FORWARD | 1 = REVERSE)
#           -CHROM             : (STR)  THE CHROMOSOME MAPPED TO THE READ
#           -START             : (INT)  THE START POSITION OF THE READ
#          -CIGAR             : (STR)  THE CIGAR SEQUENCE OF THE READ
#           -SEQ               : (STR)  THE SEQUENCE OF THE READ
#           -QUAL              : (STR)  THE QUALITY STRING OF THE READ
#           -QUALS             : (DICT) A DICTIONNARY FOR THE CONVERSION OF THE QUALITIES FROM ASCII TO INT
#           -MIN_BASE_QUALITY  : (INT)  MINIMUM QUALITY SCORE OF THE BASE FOR IT TO ADDED TO THE PILEUP DICTIONNARY
#          -ALL_UMIS          : (DICT) A DICTIONNARY CONTAINING UMIS INDEXES
#
# VALUE : NONE
#          


def AddRead(pileup, umi, strand, chrom, start, cigar, seq, qual, quals, MIN_BASE_QUALITY, ALL_UMIS):

    # split the cigar sequence using any of the five characters as a delimiter:
    # M (match/mismatch), S (soft clip), I (insertion), 
    # D (deletion) or H (hard clip)
    # this will return a list with the lengths of the consecutive events
    cigar_lengths = re.split('M|S|I|D|H', cigar)
    
    # try removing the empty element at the end of the list created by splitting
    # if last character is empty
    try:
        # remove it
        cigar_lengths.remove('')

    # if last character is not empty
    # a character not in [M,I,S,D,H] is in the cigar element
    # don't know what to do => skip read 
    except:
        print('\n')
        PrintTime("error", "\t\tUnexpected character found in CIGAR ("+cigar+")... Read skipped !\n")
        return False


    # split the cigar sequence by any number
    # this will return a list with the consecutive events
    cigar_chars = re.split('[0-9]+', cigar)
    
    # remove the empty element at the end of the list created by splitting
    cigar_chars.remove('')
    
    # remove N fragment and its length from cigar if found     if 'N' in cigar:
    if "N" in cigar:
        cigar_lengths.remove(cigar_lengths[cigar_chars.index('N')])
        cigar_chars.remove('N')
    
    # initialize a counter to increment position only without advancing in the sequence
    cursor_pos = 0
    
    # initialize a counter to advance in the sequence without incrementing the position
    cursor_seq = 0     

    # the need of two seperate is for cases when the event is an insertion as we need to
    # advance in the sequence of the read without incrementing the position or when the 
    # event is a deletion as we need to increment the position without advancing in the
    # sequence of the read



    # for each cigar event in the cigar events list (M|S|I|D|H)
    for i in range(0, len(cigar_chars)):
        
        # if the event is a match/mismatch
        if cigar_chars[i] == "M":
            
            # get the length of the match/mismatch event
            maxx = int(cigar_lengths[i])

            # call the specific function responsible of parsing a match/mismatch segment of the read 
            value = AddMatches(pileup, umi, chrom, start, seq, strand, cursor_pos, cursor_seq, maxx, qual, quals, MIN_BASE_QUALITY, ALL_UMIS)
            
            # get the returned values to continue from the position where the event ends
            position = value[0]
            cursor_seq = value[1]
            cursor_pos = value[2]
        
        # if the event is an insertion
        elif cigar_chars[i] == "I":
            
            # get the length of the insertion event
            maxx = int(cigar_lengths[i])

            # try to get position
            # normally, an read does not start with an insertion so the position of the insertion
            # is the last position of the previous cigar element 
            try:
                test = position

            # in some cases, the read starts with an insertion, therefore that position is not defined
            # in this case, an exception will be throw in order to attribute the start value to position
            except:
                position = start

            # call the specific function responsible of parsing an inserted segment of the read 
            value = AddInsertions(pileup, umi, chrom, position, seq, strand, cursor_pos, cursor_seq, maxx, qual, quals, MIN_BASE_QUALITY, ALL_UMIS)
            
            # get the returned values to continue from the position where the event ends
            position = value[0]
            cursor_seq = value[1]
            cursor_pos = value[2]

        # if the event is a deletion
        elif cigar_chars[i] == "D":

            # get the length of the deletion event
            maxx = int(cigar_lengths[i])

            # call the specific function responsible of parsing a deleted segment of the read 
            value = AddDeletions(pileup, umi, chrom, start, seq, strand, cursor_pos, cursor_seq, maxx, ALL_UMIS)

            # get the returned values to continue from the position where the event ends
            position = value[0]
            cursor_seq = value[1]
            cursor_pos = value[2]
        
        # if the event is a soft clip
        # soft clipping removes the clipped sequence from the read without correcting read start position
        elif cigar_chars[i] == "S":
            
            # test if the soft clip occurs at the beginning of the read
            try:
                # if not beginning of the read => cursor_pos > 0
                test = 20 / cursor_pos
                # in this case, do nothing since the soft clipping occured at the end of the read
                # since a clipping event can only occurs at the start or at the end of the read 
                continue

            # if beginning of the read => cursor_pos = 0 => exception 
            except:
                # get length of the clipped sequence
                clipped = int(cigar_lengths[i])
                # correct the read start position
                start -= clipped

                # increment the 2 cursors and continue
                cursor_pos += clipped
                cursor_seq += clipped

        # if the event is a hard clip
        # just continue (hard clipping corrects read start position and removes the clipped sequence from the read)
        elif cigar_chars[i] == "H":
            continue




# THIS FUNCTION ALLOWS TO MERGE MULTIPLE SUBPILEUPS TO CREATE ONE WHOLE PILEUP.
#
# INPUT : -PILEUP            : (DICT) THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
#       : -PILEUPS           : (LIST) A LIST CONTAINING ALL THE SUBPILEUPS TO BE MERGED
#       : -SUBFILES          : (LIST) THE LIST OF SUBFILES NAMES
#
# VALUE : -POSITION          : (INT)  THE FINAL POSITION THAT'S BEEN PARSED IN THE READ
#          -CURSOR_POS        : (INT)  THE CURRENT POSITION IN THE READ
#           -CURSOR_SEQ        : (INT)  THE CURRENT POSITION IN THE SEQUENCE OF THE READ
#    



def MergeSubPileups(pileup, subFiles, OUTPUT):

    # check that the number of subfiles == the number of given cub pileups
    # the 2 number must be equal. If not, an error occured in the creation
    # of some pileups ==> fatal error ==> script exists with error
    #if len(pileups) != len(subFiles):
    #    print('\n')
    #    PrintTime("error", "\tError while attempting to merge pileups : lengths differ!\n\t\t\tExiting...")

    # if the pileup lists == 1 ==> no merging has to be done ==> return the pileup directly
    #if len(pileups) == 1:
    #    return pileups[0]
    

    ### remove subpileups after loading them
    for subFile in subFiles:
        samName = subFile.split("/")[-1]
        pileupFile = OUTPUT+"/"+samName.replace('.sam', '.pileup')
        #os.remove(pileupFile)
    


    # loop through the sub pileups list
    
    #for p in pileups:
    for subFile in subFiles:
        samName = subFile.split("/")[-1]
        pileupFile = OUTPUT+"/"+samName.replace('.sam', '.pileup')
        #print(pileupFile)

        with open(pileupFile, 'rb') as fp:
            p = msgpack.unpack(fp, encoding="utf-8")
        os.remove(pileupFile)
        #perc_memory=list(psutil.virtual_memory())[2]
        # for each chromosome
        for chrom in p.keys():

            # for each position and its counters
            for position, composition in p[chrom].items():

            # for each position and its counters

                
                # increment the A forward counters
                pileup[chrom][position]['A'][0] += p[chrom][position]['A'][0]
                # increment the A reverse counters
                pileup[chrom][position]['A'][1] += p[chrom][position]['A'][1]
                # add the umis to the A unique umis set 
                pileup[chrom][position]['A'][2] += p[chrom][position]['A'][2]


                # increment the A forward counters
                pileup[chrom][position]['C'][0] += p[chrom][position]['C'][0]
                # increment the C reverse counters
                pileup[chrom][position]['C'][1] += p[chrom][position]['C'][1]
                # add the umis to the C unique umis set    
                pileup[chrom][position]['C'][2] += p[chrom][position]['C'][2]


                # increment the G forward counters
                pileup[chrom][position]['G'][0] += p[chrom][position]['G'][0]
                # increment the G reverse counters
                pileup[chrom][position]['G'][1] += p[chrom][position]['G'][1]
                # add the umis to the G unique umis set    
                pileup[chrom][position]['G'][2] += p[chrom][position]['G'][2]


                # increment the T forward counters
                pileup[chrom][position]['T'][0] += p[chrom][position]['T'][0]
                # increment the T reverse counters
                pileup[chrom][position]['T'][1] += p[chrom][position]['T'][1]
                # add the umis to the T unique umis set    
                pileup[chrom][position]['T'][2] += p[chrom][position]['T'][2]


                # add the insertions on the sub pileup to the big pileup
                # for each insertion in the insertion dictionnary
                for ins in p[chrom][position]['in'].keys():
                    # try-except block
                    try:
                        # if the insertion is already present ==> increment this insertion counters
                        pileup[chrom][position]['in'][ins][0] += p[chrom][position]['in'][ins][0]
                        pileup[chrom][position]['in'][ins][1] += p[chrom][position]['in'][ins][1]
                        pileup[chrom][position]['in'][ins][2] += p[chrom][position]['in'][ins][2]
                    except:
                        # if insertion is not present in the insertions dict ==> it has to be inserted
                        # with its own counters and unique umis set
                        pileup[chrom][position]['in'][ins] = [p[chrom][position]['in'][ins][0], p[chrom][position]['in'][ins][1], p[chrom][position]['in'][ins][2]]


                # add the deletions on the sub pileup to the big pileup
                # for each deletion in the deletion dictionnary
                for dele in p[chrom][position]['del'].keys():
                    # try-except block
                    try:
                        # if the deletion is already present ==> increment this deletion counters
                        pileup[chrom][position]['del'][dele][0] += p[chrom][position]['del'][dele][0]
                        pileup[chrom][position]['del'][dele][1] += p[chrom][position]['del'][dele][1]
                        pileup[chrom][position]['del'][dele][2] += p[chrom][position]['del'][dele][2]
                    except:
                        # if deletion is not present in the deletions dict ==> it has to be inserted
                        # with its own counters and unique umis set
                        pileup[chrom][position]['del'][dele] = [p[chrom][position]['del'][dele][0], p[chrom][position]['del'][dele][1], p[chrom][position]['del'][dele][2]]


                # increment the total base quality scores 
                pileup[chrom][position]['base_error_probability'] += p[chrom][position]['base_error_probability']
    
    # return the final whole pileup
    p=None
    return pileup




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
        gVCF = {}
        for chrom in full_pileup.keys():
            gVCF[chrom] = {}
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

    full_pileup=None
    pileup=None
    FASTA=None
    
    return kept




# THIS FUNCTION WILL WRITE THE PILEUP DICTIONNARY IN A CSV FILE IN THE OUTPUT DIRECTORY 
#
# INPUT : 
#         -PILEUP            : (DICT) THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
#         -SAM               : (STR)  THE PATH OF THE SAM FILE
#         -OUTPUT            : (OUTPUT) THE DIRECTORY IN WHICH THE CSV FILE IS TO BE CREATED
#
# VALUE : NONE
#    



def PileupToCSV(pileup, SAM, OUTPUT):
    
    print("\n")
    PrintTime('console', "\tConverting PILEUP to CSV...")

    # create a csv file in the output dir with the same name as the SAM file
    csv = open(OUTPUT+"/"+SAM.replace(".sam", ".pileup.csv").split("/")[-1], "w")
    # create the header of the file
    header = "index,chr,position,A,C,G,T,in,del,total_reads,qScore,base_error_probability,nb_unique_umi,unique_umi_list"
    # write the header into the file
    csv.write(header+"\n")

    # counters to track progress
    currentPos = 1.0
    lastProgress = 0.0
    totalPos = GetPileupLength(pileup)


    # loop through the PILEUP dictionnary
    # for each chromosome
    for chrom , infos in pileup.items():
        # for each position | composition = counts
        for position, composition in infos.items():

            # function that displays progress efficiently
            lastProgress = PrintProgress(currentPos, totalPos, lastProgress)

            # create a set that will contain unique umis
            unique_umis = set()
            
            # get n insertions at this location for the + and the - strands
            # update the unique_umis set with the insertions umis
            n_ins_pos = 0
            n_ins_neg = 0
            for value in composition['in'].values():
                n_ins_pos += value[0]
                n_ins_neg += value[1]
                unique_umis.update(value[2])
            

            # get n deletions at this location for the + and the - strands
            # update the unique_umis set with the deletions umis
            n_del_pos = 0
            n_del_neg = 0
            for value in composition['del'].values():
                n_del_pos += value[0]
                n_del_neg += value[1]
                unique_umis.update(value[2])

            # get total bases count on the + strand
            total_reads_pos = composition['A'][0]+composition['C'][0]+composition['G'][0]+composition['T'][0]
            
            # get total bases count on the - strand
            total_reads_neg = composition['A'][1]+composition['C'][1]+composition['G'][1]+composition['T'][1]
            
            # get total bases count
            total_reads = total_reads_pos+total_reads_neg
            
            # update the unique_umis set with the 4 bases umis
            for base in ['A', 'C', 'G', 'T']:
                unique_umis.update(composition[base][2])

            
            # create the line to be written in the CSV file for each position  and write it to the file
            line = [chrom+"|"+str(position), chrom, str(position), str(composition['A'][0])+"/"+str(composition['A'][1]), str(composition['C'][0])+"/"+str(composition['C'][1]), str(composition['G'][0])+"/"+str(composition['G'][1]), str(composition['T'][0])+"/"+str(composition['T'][1]),str(n_ins_pos)+"/"+str(n_ins_neg), str(n_del_pos)+"/"+str(n_del_neg), str(total_reads), str(composition['qScore']), str(composition['base_error_probability']), str(len(unique_umis)), ";".join(unique_umis)]
            line = ",".join(line)
            csv.write(line+"\n")

            # increment counter to track progress            
            currentPos += 1.0


    # close csv file
    csv.close()

    # print done
    print("\n")
    PrintTime('console', "\tDone")







def Call(config):


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

        # preprocess reads
        # if more then one core is to be used, separate the input into subfiles
        subFiles = PreprocessReads(SAM, totalLines, CORES)


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
                base = "python3 functions/TreatReads.py "+subFile+" "+BED+" "+str(MIN_BASE_QUALITY)+" "+str(MIN_READ_QUALITY)+" "+str(MIN_MAPPING_QUALITY)+" & "
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
                    pileupFile = "pileups/"+samName.replace('.sam', '.pileup')
                    try:
                        p = msgpack.unpack(open(pileupFile, 'rb'), encoding="utf-8")
                        pileups.append(p)
                    except:
                        finished = False
                        time.sleep(1)


            # merge sub pileups to obtain whole pileup
            pileup = MergeSubPileups(pileup, pileups, subFiles)


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
            #pileup = SortPileup(pileup)

        PrintTime('console', "\tDone")






    full_pileup = CopyPileup(pileup)



    ### Poisson modeling to filter positions
    result = FilterPositions(pileup, ALPHA)
    pileup = result[0]
    potential = result[1]



    ### call final variants
    finalVariants = CallVariants(pileup, f, STRAND_BIAS_METHOD, MAX_STRAND_BIAS, MIN_VARIANT_UMI, MAX_HP_LENGTH)



    ### Writing results to VCF
    final = Output(full_pileup, pileup, finalVariants, INPUT, SAM, BED, FASTA, OUTPUT, MIN_BASE_QUALITY, MIN_READ_QUALITY, MIN_MAPPING_QUALITY, MIN_VARIANT_UMI, STRAND_BIAS_METHOD, MAX_STRAND_BIAS, CORES, ALPHA, MAX_HP_LENGTH, gVCF)



    # calculate and display stats
    CalculateStats(pileup, potential, final)































































# THIS FUNCTION ALLOWS TO PARSE A CIGAR SEGMENT CONTAINING ONLY DELETED BASES AND INCREMENT
# A SPECIFIC PILEUP DICTIONNARY ACCORDINGLY
#
# INPUT : 
#           -BED               : (STR)  THE PATH OF THE BED FILE
#
# VALUE : -PILEUP            : (DICT) A DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
#        



def ParseBED(bed):    

    # create an empty dictionnary for pileup counters
    pileup = {}

    # read the bed file line by line
    with open(bed,'r+') as fbed:
        for line in fbed:
            if not (line[0] in ['#','@']): 

                    # split the line using tab as delimiter
                line = line.split('\t')

                if len(line) >= 3:
                    # first element is the chromosome
                    chrom = line[0]

                    # second element is the start of the region
                    # third element is the end of the region
                    # create a list containing start and end 
                    # sort to ensure that start < end
                    limits = [int(line[1]), int(line[2])]
                    limits.sort()

                    # for each position in the [start;end] interval
                    for pos in range(limits[0], limits[1]+1):
                        
                        # try to add the position to the pileup for the specified chromosome
                        # if the chromosome was already aded to the pileup
                        try:
                            # add the position to the chrom with the default counters
                            pileup[chrom][pos] = { 'A': [0, 0, []], 'C': [0, 0, []], 'G': [0, 0, []], 'T': [0, 0, []], 'in': {}, 'del': {}, 'base_error_probability': 0 } 
                        except:
                            # add the chrom to the pileup with the first position and the default counters
                            pileup[chrom] = { pos: { 'A': [0, 0, []], 'C': [0, 0, []], 'G': [0, 0, []], 'T': [0, 0, []], 'in': {}, 'del': {}, 'base_error_probability': 0 } }



    # the pileup dictionnary contains :
    # for each chromosome:
    #         for each position:
    #        {
    #            'A'   : 0 | 0 | list()
    #            'C'   : 0 | 0 | list()
    #            'T'   : 0 | 0 | list()
    #            'G'   : 0 | 0 | list()
    #            'in'  :     {}
    #            'del' :     {}
    #        }
    # for the four normal bases, this structure 0 | 0 | set() is a list containing 
    # a counter for the forward strand, a counter for the reverse strand, and a set
    # that will contain the UMIs of the reads. 
    # a set was used instead of a list to contain unique UMIs only  
    #
    # the insertion dictionnary will be updated each time an insertion is encountered
    # and will end up looking like this:
    # {
    #     'inserted_seq_1' : N1_+ | N1_- | list()
    #     'inserted_seq_2' : N2_+ | N2_- | list()
    #     'inserted_seq_n' : Nn_+ | Nn_- | list()
    # }
    #
    # the deletion dictionnary will be updated each time a deletion is encountered
    # and will end up looking like this:
    # {
    #     len of the deleted segment_1 : N1_+ | N1_- | list()
    #     len of the deleted segment_2 : N2_+ | N2_- | list()
    #     len of the deleted segment_n : Nn_+ | Nn_- | list()
    # }    
    # 
    # the base_error_probability by bases qscores at each position
    # at the end, this value will be divided by depth to get the mean qscore a the location
    # the the mean qscore will be used to estimate the noise at the location

    # sort pileup
    #pileup = SortPileup(pileup)
    
    # return the pileup
    return pileup





# THIS FUNCTION ALLOWS TO PARSE A CIGAR SEGMENT CONTAINING ONLY INSERETD BASES AND INCREMENT
# A SPECIFIC PILEUP DICTIONNARY ACCORDINGLY
#
# INPUT : -PILEUP            : (DICT) THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
#           -UMI               : (STR)  UMI SEQUENCE OF THE READ
#           -CHROM             : (STR)  THE CHROMOSOME MAPPED TO THE READ
#           -POSITION          : (INT)  THE POSITION OF THE INSERTION SITE IN THE READ
#           -SEQ               : (STR)  THE SEQUENCE OF THE READ
#           -STRAND            : (INT)  THE STRAND OF THE READ (0 = FORWARD | 1 = REVERSE)
#           -CURSOR_POS        : (INT)  THE CURRENT POSITION IN THE READ
#           -CURSOR_SEQ        : (INT)  THE CURRENT POSITION IN THE SEQUENCE OF THE READ
#           -MAXX              : (INT)  THE LENGTH OF THE CIGAR 'I' ELEMENT
#           -QUAL              : (STR)  THE QUALITY STRING OF THE READ
#           -QUALS             : (DICT) A DICTIONNARY FOR THE CONVERSION OF THE QUALITIES FROM ASCII TO INT
#           -MIN_BASE_QUALITY  : (INT)  MINIMUM QUALITY SCORE OF THE BASE FOR IT TO ADDED TO THE PILEUP DICTIONNARY
#          -ALL_UMIS          : (DICT) A DICTIONNARY CONTAINING UMIS INDEXES
#
# VALUE : -POSITION          : (INT)  THE FINAL POSITION THAT'S BEEN PARSED IN THE READ
#          -CURSOR_POS        : (INT)  THE CURRENT POSITION IN THE READ
#           -CURSOR_SEQ        : (INT)  THE CURRENT POSITION IN THE SEQUENCE OF THE READ

def AddInsertions(pileup, umi, chrom, position, seq, strand, cursor_pos, cursor_seq, maxx, qual, quals, MIN_BASE_QUALITY, ALL_UMIS):
    
    # create an empty string to contain the inserted sequence
    # create a mean qscore for the inserted sequence 
    inserted_seq = ""
    inserted_qscore = 0

    # try because position in read could not be in BED file
    # therefore, the position could not be in the PILEUP dictionnary
    # if position not in BED, skip the read
    try:
        # testing if position + 1 if in the BED file
        # position is the position of the last matched/mismatched position
        # returned by the AddMatches function, therefore the insertion 
        # occurs at position+1
        test = pileup[chrom][position+1]

        # looping through the inserted sequence
        # maxx = length of the "I" cigar element => length of the inserted_sequence
        for i in range(0, maxx):
            # get base at i
            base = seq[cursor_seq]
            # get base qual at i
            baseQual = qual[cursor_seq]
            # convert base qual to base qscore
            base_qscore = quals[baseQual]

            # add the base at i in seq to the inserted sequence
            inserted_seq += base
            # increment the inserted_seq qscore with the quality 
            # of the inserted base
            inserted_qscore += base_qscore
            # advance in the read sequence
            cursor_seq += 1


        # calculate qscore for the inserted sequence corresponding
        # to the mean of the qscores of the bases in the inserted seq
        inserted_qscore = inserted_qscore/len(inserted_seq)


            
        # check if the quality of the inserted sequence >= minimum base quality
        if inserted_qscore >= MIN_BASE_QUALITY:

            # if this insertion was already seen at this location
            try:
                # try to increment the corresponding counter in the PILEUP dictionnary
                pileup[chrom][position+1]['in'][inserted_seq][strand] += 1
                # add the umi to the corresponding set specific to the inserted sequence
                pileup[chrom][position+1]['in'][inserted_seq][2].append(umi)
                # increment the qScore of this position
                pileup[chrom][position+1]['base_error_probability'] += inserted_qscore

            # if this is the first time we see this inserted sequence at this position, 
            # an entry should be created with its appropriated structure (2 counters + set)  
            except:
                # create the entry in PILEUP dictionnary
                pileup[chrom][position+1]['in'][inserted_seq] = [0, 0, [umi]]
                # increment the corresponding counter in the PILEUP dictionnary
                pileup[chrom][position+1]['in'][inserted_seq][strand] += 1
                # increment the qScore of this position
                pileup[chrom][position+1]['base_error_probability'] += inserted_qscore


            cursor_pos += 1
            cursor_seq += 1

    
    # if position not in BED, skip the read
    except:
        cursor_seq += maxx

    # return position, cursor_seq and cursor_pos to continue from where we left off
    return [position, cursor_seq, cursor_pos]





def TreatReads(SAM, BED, MIN_BASE_QUALITY, MIN_READ_QUALITY, MIN_MAPPING_QUALITY, ALL_UMIS):

    quals_str = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"
    quals = {}


    i=0
    for q in quals_str:
        quals[quals_str[i]] = i
        i += 1


    pileup = ParseBED(BED)


    currentLine = 1.0
    lastProgress = 0.0
    totalLines = GetTotalLines(SAM)

    lastLine = {}
    lastUMI = ""


    for line in open(SAM):

        lastProgress = PrintProgress(currentLine, totalLines, lastProgress)

        currentLine += 1


        line = line.split('\t')

        flag = int(line[1])
        strand = int("{0:b}".format(int(flag))[-5])
        firstInPair = bool(int("{0:b}".format(int(flag))[-7]))
        chrom = line[2]
        pos = line[3]
        mapq = int(line[4])
        cigar = line[5]
        seq = line[9]
        qual = line[10]

        #added UMI based on position as well. Better for hybrid capture with low number of UMIs
        #UMI will be filtered based on UMI sequence AND Position of the read
        #Original umi = line[0].split('_')[-1]+str(pos)[-3:]
        try:
            umi = line[0].split('_')[-1]+str(pos)[-3:]
        except:
            continue


        if ReadIsValid(flag, chrom, pos, umi, cigar, mapq, MIN_MAPPING_QUALITY, qual, quals, MIN_READ_QUALITY):

            AddRead(pileup, umi, strand, chrom, int(pos), cigar, seq, qual, quals, MIN_BASE_QUALITY, ALL_UMIS)




    return pileup



# THIS FUNCTION ALLOWS TO CONFIGURE THE PARAMETERS FROM THE COMMAND LINE
# AND PROPERLY PASS THEM TO THE CORRESPONDING FUNCTION CALLS
#
# INPUT : 
#           NONE
#
# VALUE : -CONFIG              : (DICT)  A DICTIONNARY THAT CONTAINS ALL NECESSARY THE PARAMETERS AND THEIR VALUES
#        





# THIS FUNCTION ALLOWS TO CONFIGURE THE PARAMETERS FROM THE COMMAND LINE
# AND PROPERLY PASS THEM TO THE CORRESPONDING FUNCTION CALLS
#
# INPUT : 
#           NONE
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





def GenerateVCFHeader(INPUT, BED, FASTA, OUTPUT, MIN_BASE_QUALITY, MIN_READ_QUALITY, MIN_MAPPING_QUALITY, MIN_VARIANT_UMI, SB_METHOD, MAX_STRAND_BIAS, CORES, ALPHA, MAX_HP_LENGTH):
    
    date = datetime.datetime.now().strftime("%Y-%m-%d")

    l1 = "##fileformat=VCFv4.2\n"
    l2 = "##fileDate="+date+"\n"
    l3 = '##source="UMI-VarCal"\n'
    l4 = '##command="input=\''+INPUT+'\', bed=\''+BED+'\', out=\''+OUTPUT+'\', min_base_quality='+str(MIN_BASE_QUALITY)+', min_read_quality='+str(MIN_READ_QUALITY)+', min_mapping_quality='+str(MIN_MAPPING_QUALITY)+', min_variant_umi='+str(MIN_VARIANT_UMI)+', strand_bias_method=\''+str(SB_METHOD)+'\', max_strand_bias='+str(MAX_STRAND_BIAS)+', alpha='+str(ALPHA)+', max_hp_length='+str(MAX_HP_LENGTH)+', cores='+str(CORES)+'\"\n'
    l5 = '##reference="'+FASTA+'"\n'
    rest = """##FILTER=<ID=ID,Description="Generic filter">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency based on allele count/depth">
##INFO=<ID=AO,Number=2,Type=Integer,Description="Alternate allele observation count">
##INFO=<ID=DP,Number=3,Type=Integer,Description="Total read depth at the locus">
##INFO=<ID=HP,Number=4,Type=Integer,Description="Homopolymer length at the locus">
##INFO=<ID=TYPE,Number=5,Type=String,Description="The variation type: either SNV, INS or DEL">
##INFO=<ID=CONF,Number=6,Type=String,Description="Confidence of the variant. It has 5 levels: low(1/5), average(2/5), high(3/5), strong(4/5) or certain(5)">
#CHROM    POS    ID    REF    ALT    Q-VALUE    FILTER    INFO
"""

    return l1+l2+l3+l4+l5+rest





def GenerateGVCFHeader(INPUT, BED, FASTA, OUTPUT, MIN_BASE_QUALITY, MIN_READ_QUALITY, MIN_MAPPING_QUALITY, MIN_VARIANT_UMI, SB_METHOD, MAX_STRAND_BIAS, CORES, ALPHA, MAX_HP_LENGTH):
    
    date = datetime.datetime.now().strftime("%Y-%m-%d")

    l1 = "##fileformat=VCFv4.2\n"
    l2 = "##fileDate="+date+"\n"
    l3 = '##source="UMI-VarCal"\n'
    l4 = '##command="input=\''+INPUT+'\', bed=\''+BED+'\', out=\''+OUTPUT+'\', min_base_quality='+str(MIN_BASE_QUALITY)+', min_read_quality='+str(MIN_READ_QUALITY)+', min_mapping_quality='+str(MIN_MAPPING_QUALITY)+', min_variant_umi='+str(MIN_VARIANT_UMI)+', strand_bias_method=\''+str(SB_METHOD)+'\', max_strand_bias='+str(MAX_STRAND_BIAS)+', alpha='+str(ALPHA)+', max_hp_length='+str(MAX_HP_LENGTH)+', cores='+str(CORES)+'\"\n'
    l5 = '##reference="'+FASTA+'"\n'
    rest = """##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
##FILTER=<ID=ID,Description="Generic filter">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency based on allele count/depth">
##INFO=<ID=AO,Number=2,Type=Integer,Description="Alternate allele observation count">
##INFO=<ID=DP,Number=3,Type=Integer,Description="Total read depth at the locus">
##INFO=<ID=HP,Number=4,Type=Integer,Description="Homopolymer length at the locus">
##INFO=<ID=TYPE,Number=5,Type=String,Description="The variation type: either SNV, INS or DEL">
##INFO=<ID=CONF,Number=6,Type=String,Description="Confidence of the variant. It has 5 levels: low(1/5), average(2/5), high(3/5), strong(4/5) or certain(5)">
##INFO=<ID=END,Number=7,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=MEAN_DP,Number=8,Type=Integer,Description="Average depth on the interval">
##INFO=<ID=MEAN_QSCORE,Number=9,Type=Float,,Description="Average base quality score on the interval">
#CHROM    POS    ID    REF    ALT    Q-VALUE    FILTER    INFO
"""

    return l1+l2+l3+l4+l5+rest





# THIS FUNCTION ALLOWS TO PARSE A CIGAR SEGMENT CONTAINING ONLY MATCHED/MISMATCHED BASES AND INCREMENT
# A SPECIFIC PILEUP DICTIONNARY ACCORDINGLY
#
# INPUT : -PILEUP            : (DICT) THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
#           -UMI               : (STR)  UMI SEQUENCE OF THE READ
#           -CHROM             : (STR)  THE CHROMOSOME MAPPED TO THE READ
#           -START             : (INT)  THE START POSITION OF THE READ
#           -SEQ               : (STR)  THE SEQUENCE OF THE READ
#           -STRAND            : (INT)  THE STRAND OF THE READ (0 = FORWARD | 1 = REVERSE)
#           -CURSOR_POS        : (INT)  THE CURRENT POSITION IN THE READ
#           -CURSOR_SEQ        : (INT)  THE CURRENT POSITION IN THE SEQUENCE OF THE READ
#           -MAXX              : (INT)  THE LENGTH OF THE CIGAR 'M' ELEMENT
#           -QUAL              : (STR)  THE QUALITY STRING OF THE READ
#           -QUALS             : (DICT) A DICTIONNARY FOR THE CONVERSION OF THE QUALITIES FROM ASCII TO INT
#           -MIN_BASE_QUALITY  : (INT)  MINIMUM QUALITY SCORE OF THE BASE FOR IT TO ADDED TO THE PILEUP DICTIONNARY
#          -ALL_UMIS          : (DICT) A DICTIONNARY CONTAINING UMIS INDEXES
#
# VALUE : -POSITION          : (INT)  THE FINAL POSITION THAT'S BEEN PARSED IN THE READ
#          -CURSOR_POS        : (INT)  THE CURRENT POSITION IN THE READ
#           -CURSOR_SEQ        : (INT)  THE CURRENT POSITION IN THE SEQUENCE OF THE READ


def AddMatches(pileup, umi, chrom, start, seq, strand, cursor_pos, cursor_seq, maxx, qual, quals, MIN_BASE_QUALITY, ALL_UMIS):

    # for each position between limits
    for position in range(start+cursor_pos, start+cursor_pos+maxx):

        # try because base position in read could not be in BED file
        # therefore, the position could not be in the PILEUP dictionnary
        # if position not in BED, skip the read
        try:    

            # get base, base_qual and base_qscore
            base = seq[cursor_seq]
            baseQual = qual[cursor_seq]
            base_qscore = quals[baseQual]
                
            # check if the quality of the base >= minimum base quality
            if base_qscore >= MIN_BASE_QUALITY:
                # increment the corresponding counter in the PILEUP dictionnary
                pileup[chrom][position][base][strand] += 1
                # add the umi to the corresponding set specific to the base in the PILEUP dictionnary
                pileup[chrom][position][base][2].append(umi)
                # increment the qScore of this position
                pileup[chrom][position]['base_error_probability'] += base_qscore
            

        # if position not in BED, skip the read
        except:
            pass    

        # increment the cursor_seq to allow advancing in the read sequence
        cursor_seq += 1


    # increment the cursor_pos with the length of the cigar element when done adding it
    cursor_pos += maxx


    # return position, cursor_seq and cursor_pos to continue from where we left off
    return [position, cursor_seq, cursor_pos]




# THIS FUNCTION IS CALLED TO SPLIT THE INITIAL SAM FILE INTO N SAM SUBFILES, N CORRESPONDING
# TO THE NUMBER OF CORES USED TO EXECUTE THE PROGRAM. THEREFORE, EACH CORE WILL ANALYZE A 
# SUBFILE, ALLOWING THEREFORE ALL THE SUBFILES TO BE ANALYZED AT THE SAME TIME BY ALL THE 
# PRECISED CORES TO BE USED.
#
# INPUT : 
#           -FILENAME                     : (STR)    THE PATH OF THE INITIAL SAM FILE
#           -TOTALLINES                   : (FLOAT)  TOTAL LINES IN THE FILE
#           -CORES                        : (INT)    NUMBER OF CORES TO BE USED
#
# VALUE : -SUBFILES                     : (LIST)   A LIST CONTAINING THE NAMES OF THE CREATED SUBFILES
#        
#




def PreprocessReads(fileName, totalLines, cores):

    # if only one core is to be used, no need to split the file
    # return the initial file
    if cores == 1:
        return [fileName]




    print("\n")
    PrintTime('console', "\tPreprocessing Reads...")


    # create an empty list to contain the created sub files names
    subFiles = []

    # calculate at the interval at which a new subfile is to be created
    # this allows for the subfiles to be created with approximately the
    # same size
    interval =  int(totalLines/cores)
    interval = interval if interval % 2 == 0 else interval - 1


    # set counters
    i = 0
    nLine = 0.0
    lastProgress = 0.0

    # read the file line by line
    for line in open(fileName, 'r'):

        # display progress
        lastProgress = PrintProgress(nLine+1, totalLines, lastProgress)

        # if interval value is reached
        if nLine % interval == 0:
            # if this subfile is not the last subfile
            if i != cores:
                # define the name of the subfile
                outName = fileName.replace('.sam', '_'+str(i)+'.sam')
                # create the output file
                output = open(outName, 'w')
                # add the subfile name to the subfiles list
                subFiles.append(outName)
                # increment the subfiles counter
                i += 1
        # write the line in the i subfile
        output.write(line)
        # increment the line counter
        nLine += 1


    
    print("\n")
    PrintTime('console', "\tDone")

    # return the subfiles list 
    return subFiles


# THIS FUNCTION WILL APPLY THE POISSON TEST ON EACH POSITION TO TRY AND FIND A SIGNIFCANT VARIANTS. THEN, P-VALUES ARE CORRECTED
# WITH THE BENJAMINI-HOCHBERG METHOD TO OBTAIN Q-VALUES. ONLY POSITIONS WITH Q-VALUE > ALPHA ARE KEPT
#
# INPUT : 
#         -PILEUP              : (DICT)   THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
#           -ALPHA               : (FLOAT)  TYPE 1 ERROR PROBABLITY OR ALPHA LEVEL
#
# VALUE : -PILEUP              : (DICT)   THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC



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

    #pileup=None
    # return PILEUP dictionnary

    return (pileup, potential, foundCandidates)

    

# THIS FUNCTION WILL CALCULATE THREE STATISTICS : THE PERCENTAGE OF THE BED COVERED BY THE BAM/SAM FILE, 
# THE AVERAGE DEPTH OF THE BAM/SAM FILE ACROSS ALL POSITIONS AND THE UNIFORMITY OF THE SEQUENCING (THE 
# PERCENTAGE OF THE LOCATIONS THAT HAVE A DEPTH > 0.2*AVERGAE DEPTH)
#
# INPUT : -PILEUP            : (DICT) THE DICTIONNARY CONTAINING COUNTERS THAT ARE CHROMOSOME-POSITION-BASE-STRAND SPECIFIC
#          -POTENTIAL         : (INT)  NUMBER OF POTENTIAL VARIANTS FOUND
#          -FINAL                : (INT)  NUMBER OF FINAL VARIANTS FOUND
#
# VALUE : NONE
#    


def CalculateStats(pileup, potential, final):
    print("\n")
    PrintTime('console', "\tCalculating statistics...")

    # initiate a counter for not covered positions
    not_covered = 0
    # get total number of positions in pileup
    totalPos = GetPileupLength(pileup)
    # create an empty list to stock depths
    depths = []

    # loop through the pileup items
    # for each chromosome
    for chrom , infos in pileup.items():
        # for each position | composition = counts
        for position, composition in infos.items():
            # try-except block
            # if position is covered => depth must be > 0
            try:        
                # calculate insertion counts
                n_ins_pos = 0
                n_ins_neg = 0
                for value in composition['in'].values():
                    n_ins_pos += value[0]
                    n_ins_neg += value[1]
                
                # calculate deletion counts
                n_del_pos = 0
                n_del_neg = 0
                for value in composition['del'].values():
                    n_del_pos += value[0]
                    n_del_neg += value[1]
                
                # insertions shouldn't be taken into account when calculating depth
                # n_reads = composition['A'][0]+composition['C'][0]+composition['G'][0]+composition['T'][0]+composition['A'][1]+composition['C'][1]+composition['G'][1]+composition['T'][1]+n_ins_pos+n_ins_neg+n_del_pos+n_del_neg
                n_reads = composition['A'][0]+composition['C'][0]+composition['G'][0]+composition['T'][0]+composition['A'][1]+composition['C'][1]+composition['G'][1]+composition['T'][1]+n_del_pos+n_del_neg
                # make a test to check that depth > 0
                test = 25 / n_reads
                # if test succeed => position covered => append depth to list of depths
                depths.append(n_reads)

            # if position not covered => increment not_covered counter
            except:
                not_covered += 1
    pileup=None

    # calculate coverage
    coverage = round(float(totalPos-not_covered)/float(totalPos), 5)*100
    # calculate average depth
    avg_depth = int(round(float(sum(depths))/float(len(depths)), 0))

    # calculate uniformity
    uniform = 0
    for depth in depths:
        if depth >= 0.2*avg_depth:
            uniform += 1
    uniformity = round(float(uniform)/float(len(depths)), 2)*100


    # print out stats to console
    message = "BAM/BED coverage: "+ str(coverage)+" %"
    PrintTime('green', "\t\t"+message)
    message = "Average design depth: "+ str(avg_depth)+"x"
    PrintTime('green', "\t\t"+message)
    message = "Uniformity: "+ str(uniformity)+" %"
    PrintTime('green', "\t\t"+message)
    message = "Candidate Positions: "+ str(potential)
    PrintTime('green', "\t\t"+message)
    message = "Final Variants: "+ str(final)
    PrintTime('green', "\t\t"+message)
    return (str(coverage), str(avg_depth), str(uniformity), str(potential), str(final))


# THIS FILE CONTAINS GENERAL FUNCTIONS THAT ARE NEEDED TO 
# SUCCESSFULLY CALL THE VARIANTS IN THE SAM FILE 
#


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

        message = "\tRAM USAGE    :  "+str(round(usedRAM, 2))+" GB"+"\n\n"
        PrintTime("console", message)

        time.sleep(1)
        for x in os.listdir(tmp_name):
            os.remove(tmp_name+"/"+x)
        
        os.rmdir(tmp_name)    
    
    except:
        print("\n")

    exit()



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
    new_pileup = {}
    final_pileup = {}

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
        new_pileup[chrom] = {}
        for pos in sorted_poss:
            new_pileup[chrom][pos] = pileup[chrom][pos]
    

    return new_pileup



# this function aims to create a distinct copy of the given dict 
def CopyPileup(pileup):
    new_pileup = {}
    for chrom in pileup.keys():
        new_pileup[chrom] = {}
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

#    value = set()

#    for line in open(SAM):

#        line = line.split('\t')
        
#        try:
#            value.add(line[0].split('_')[-1])

#        except:
#            continue

    
#    c = 0
#    final = {}

#    for el in value:
#        final[el] = c
#        c+=1

#    return final
    return 1




def PrintHelp(tool):

    tools = ["general", "extract", "call"]

    messages= []
    message = """
GENERAL DESCRIPTION
===============================================

Software for analyzing BAM/SAM files produced by an Illumina MiSeq sequencer using UMI (Unique Molecular Identifiers).
The software starts by building a pileup based on the regions specified in the BED file. This software integrates an  
extraction tool that can be used to extract UMI tags from the start of the reads. The software has also a calling tool
for calling variants and that outputs the results under VCF format.
Although SAM files generated by an Illumina MiSeq sequencer were used to develop this tool, it should work on BAM/SAM 
files generated by any other sequencer as long as they respect the Sequence Alignment/Map Format Specification.


:Author: Vincent SATER, LITIS EA 4108 / Centre Henri Becquerel, Rouen, France
:Author: Pierre-Julien VIAILLY, Centre Henri Becquerel, Rouen, France
:Version: 2.1.5
:Last updated: 30-01-2021
:Tags: Genomics, UMI, Illumina, Miseq, Variant Calling


There are 2 tools available in this software:

    -extract
    -call


To get help on a specific tool, type:

    python3 umi-varcal.py <tool>
    or 
    python3 umi-varcal.py <tool> --help
    or 
    python3 umi-varcal.py <tool> -h

To use a specific tool, type:

    python3 umi-varcal.py <tool> [tool arguments]\n\n"""

    messages.append(message)





    message = """
TOOL DESCRIPTION : extract
===============================================

The extraction tool allows to extract the UMI sequence from the SAM/BAM sequence 
and add it to the end of the read ID, preceeded by an '_'. 
The UMI is always extracted from the start of the sequence.
Extracting the UMI from within the sequence is not handled by this tool.


:Author: Vincent SATER, LITIS EA 4108 / Centre Henri Becquerel, Rouen, France
:Author: Pierre-Julien VIAILLY, Centre Henri Becquerel, Rouen, France
:Version: 2.1.5
:Last updated: 30-01-2021
:Tags: Genomics, UMI, Illumina, Miseq, Variant Calling, Alignment


USAGE
===============================================

python3 umi-varcal.py extract [required_arguments] [optional_arguments]


DEPENDENCIES
===============================================

-bwa
-Python module: pysam


REQUIRED ARGUMENTS
===============================================

-f, --fasta=FILE (case sensitive)
    The reference genome in FASTA format with its index present
    in the same directory.
    It is recommended that you enter the FASTA file full path
    For example : /home/my_dir/my_subdir/hg19.fasta

-i, --input=FILE (case sensitive)
    The input file SAM/BAM generated by the sequencer.
    It is recommended that you enter the input file full path
    For example : /home/my_dir/my_subdir/ind1.sam

-l, --umi_length=INT 
    umi_length : the length of the UMI sequence to extract
    The value must be a valid integer
    If not precised, the default value (umi_length = 12) will be used


OPTIONAL ARGUMENTS
===============================================

-t, --bwa_threads=INT 
    The number of cores to be used for the bwa alignment. Using more
    cores will significantly reduce the alignment time.
    The value of this parameter defaults to 1


EXAMPLE
===============================================

python3 umi-varcal.py extract -i /home/.../test.sam -f /home/.../hg19.fa -l 12 -t 5


OUTPUT
===============================================

The tool will automatically create a directory in the current directory 
called EXTRACTED. It will contain a folder named after the sample input.
In this directory, you can find the extracted R1 FASTQ, the extracted R2
FASTQ, the extracted BAM (with its index) and the extracted SAM file.\n\n"""

    messages.append(message)






    message = """
TOOL DESCRIPTION : call
===============================================

A tool for calling variants in a SAM/BAM file from which the UMI tags were extracted.
The calling tool starts by building a pileup based on the regions specified in the BED file. Given that the UMI tags
had already been extracted correctly, this tool will call variants and output the results under VCF format.
Although BAM/SAM files generated by an Illumina MiSeq sequencer were used to develop this tool, it should work on 
BAM/SAM files generated by any other sequencer as long as they respect the Sequence Alignment/Map Format Specification.


:Author: Vincent SATER, LITIS EA 4108 / Centre Henri Becquerel, Rouen, France
:Author: Pierre-Julien VIAILLY, Centre Henri Becquerel, Rouen, France
:Version: 2.1.5
:Last updated: 30-01-2021
:Tags: Genomics, UMI, Illumina, Miseq, Variant Calling, Alignment


USAGE
===============================================

python3 umi-varcal.py call [required_arguments] [optional_arguments]


DEPENDENCIES
===============================================

-Python modules:   
                    msgpack
                    pyfasta
                    pysam
                    scipy.stats
                    statsmodels.stats.multitest


REQUIRED ARGUMENTS
===============================================

-f, --fasta=FILE (case sensitive)
    The reference genome in FASTA format with its index present
    in the same directory.
    It is recommended that you enter the FASTA file full path
    For example : /home/my_dir/my_subdir/hg19.fasta

-b, --bed=FILE (case sensitive)
    The panel design file in a BED format.
    It is recommended that you enter the BED file full path
    For example : /home/my_dir/my_subdir/design.bed

-i, --input=FILE (case sensitive)
    The input file SAM/BAM generated by the sequencer and from which
    the UMI tags have already been extracted from the sequences and
    and added at the end of the read ID, preceded by an '_'. This format
    is the only format accepted by the variant caller.
    It is recommended that you enter the input file full path
    For example : /home/my_dir/my_subdir/ind1.sam



OPTIONAL ARGUMENTS
===============================================

-p, --pileup=FILE (case sensitive)
    The PILEUP file for the SAM file previously dumped during
    a older analysis with the same BED and FASTA files. This 
    is not the pileup format generated by SAMtools. If specified,
    the tool will attempt te retrieve the counts of each position
    from the PILEUP file instead of rebuilding the pileup from 
    scratch, a process that is relatively fast and can save time.
    It is recommended that you enter the PILEUP file full path
    For example : /home/my_dir/my_subdir/ind1.pileup

-o, --output=DIR (case sensitive)
    The output directory in which your files will be generated.
    It is recommended that you enter the output dir full path
    For example : /home/my_dir/my_subdir/out_dir

-c, --cores=INT
    The number of cores to be used for the analysis. Using more
    cores will significantly reduce analysis time but RAM usage
    will be higher (up to 3x higher)
    The value of this parameter defaults to 1

--min_base_quality=INT
    The minimum value for which a base with a certain quality is
    considered as "good". Only bases with a quality >= min_base_quality
    will be considered.
    The value of this parameter defaults to 10 

--min_read_quality=INT
    The minimum value for which a read with a certain mean quality is
    considered as "good". Only reads with a quality >= min_read_quality
    will be considered.
    The value of this parameter defaults to 20 

--min_mapping_quality=INT
    The minimum value for which a read with a certain mapping 
    quality is considered as "good". Only bases with a mapping
    quality >= min_mapping_quality will be considered.
    The value of this parameter defaults to 20

--min_variant_depth=INT
    The minimum count for a variant with a certain raw count to be
    called. Only variants with a raw count >= min_variant_depth 
    will be called.
    The value of this parameter defaults to 5

--alpha=INT/FLOAT
    This is the type I error probability known as alpha level.
    Positions with a p-value/q-value >= alpha will be considered as
    noisy and variants in these positions won't be called.
    The value of this parameter defaults to 0.05 (5%)

--min_variant_umi=INT
    The minimum count for a variant with a certain unique UMIs 
    count to be called. Only variants with a unique UMI count >= 
    min_variant_umi will be called.
    The value of this parameter defaults to 5

--strand_bias_method=STR
    This parameter lets you specifiy the method you want to use
    to calculate variant strand bias. The value can only be 
    "default" or "torrent_suite". Choosing any of the methods 
    can have minimal effect on the final variants called and the SB
    value differs: the torrent suite method generates SB values [0.5;1]
    but the default method generates SB values [0;+Inf].
    Obviously, this parameter defaults to the value "default"

--max_strand_bias=INT/FLOAT
    The minimum strand bias value for a variant with a certain  
    negative and positive coverage to be called. Only variants 
    with a strand bias <= max_strand_bias will be called.
    The value of this parameter defaults to 1.0 for the default
    method and to 0.743 for the torrent suite method

--max_hp_length=INT
    The maximum homopolymer length for a variant at a certain  
    position in an homopolymer region to be called. Only variants 
    with a homopolymer length <= max_hp_length will be called.
    The value of this parameter defaults to 10

--gvcf=BOOL
    This is an experimental feature that is turned off by default
    and that, if set to True, will generate a gVCF file along
    the normal VCF file

--keep_pileup=BOOL
    This lets the user keep or delete the generated pileup file.
    By default, this parameter is set to True


EXAMPLE
===============================================

python3 umi-varcal.py call -i /home/.../test.sam -b /home/.../design.bed -f /home/.../hg19.fa -p /home/.../test.pileup -o /home/.../out -c 3 --min_base_quality 15 --min_read_quality 25 --min_mapping_quality 40 --min_variant_depth 8 --min_variant_umi 7 --strand_bias_method default --max_strand_bias 0.8 --alpha 0.01 --max_hp_length 5 --gvcf True --keep_pileup False


OUTPUT
===============================================

1-  A VCF file (<INPUT_NAME>.vcf) containing all the variants that
    were successfully called.

2-  A gVCF file (<INPUT_NAME>.gvcf) containing all the variants that
    were successfully called + informations for the blocks in which
    no variant was detected.
    
3-  A VARIANTS file (<INPUT_NAME>.variants) containing all the variants
    that were successfully called. This file contains the same variants of
    the VCF file, in addition to detailed metrics about the variants.

4-  A PILEUP file (<INPUT_NAME>.pileup) corresponding to the entire pileup dumped.
    This file is produced in the output directory. This file can be used to skip 
    the pileup building and load it instead if the analysis was already done.\n\n"""

    messages.append(message)



    PrintProgramName()
    print(messages[tools.index(tool)])

    exit()

