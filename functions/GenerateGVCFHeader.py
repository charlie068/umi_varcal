
import datetime

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
#CHROM	POS	ID	REF	ALT	Q-VALUE	FILTER	INFO
"""

	return l1+l2+l3+l4+l5+rest
