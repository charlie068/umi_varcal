#Splitting your bed file
	gatk SplitIntervals -R fasta_file \
	--subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW -L targets.bed --scatter-count number_of_files -O ./

#Sometimes, you have 3 files from the sequencer, Freads, Rreads and fastq with UMIs.
#To concatenate UMI sequences on the sequences headers, one way is to :
#fastq1, fastq3 are Freads and  Rreadsm 
#fastq2 is the fastq with UMI
#This gives fastqumi1 and fastqumi2

	awk -v FS="\t" -v OFS="\t" 'NR==FNR {split($1, id, " "); umi[id[1]]=$2; next;} {split($1, id, " "); $1=id[1]"_"umi[id[1]]" "id[2]; print $0}'  <(zcat $fastq2|paste - - - -) <(zcat $fastq1|paste - - - -)|tr "\t" "\n" | bgzip -@4 > $fastqumi1 & \
				 awk -v FS="\t" -v OFS="\t" 'NR==FNR {split($1, id, " "); umi[id[1]]=$2; next;} {split($1, id, " "); $1=id[1]"_"umi[id[1]]" "id[2]; print $0}'  <(zcat $fastq2|paste - - - -) <(zcat $fastq3|paste - - - -)|tr "\t" "\n" | bgzip -@4 > $fastqumi2

#bwamem2 can be used for the alignement:
	bwa-mem2 mem -M -t 24  $reffile $fastqumi1 $fastqumi2 > ${pathsample}/${samplename}_bwa.sam
	
	#Sort based on ccordinates and index bam file
	samtools sort -T ${pathsample} -@24 ${pathsample}/${samplename}_bwa.sam -o ${pathsample}/${samplename}_sorted.bam
	samtools index ${pathsample}/${samplename}_sorted.bam
	rm ${pathsample}/${samplename}_bwa.sam 


#Varcal
python3 /home/jeancharles/bioinf/umi_varcal2/umi-varcal.py call -i ${pathsample}/${samplename}_sorted.bam -b $t0 \
-f $reffile -o ${pathsample}/out -c 24 --keep_pileup False --min_variant_umi 5 --min_mapping_quality 30  \
--min_base_quality 20  --min_read_quality 30 

#You can assemble the files with a little script like
	assemble_vcf_files.py

#and sort with 
	grep '^#' original_vcf_file  > ${pathout}/${samplename}.outsorted.vcf
	cat $f  | sort -u -k1,1V -k2,2n >> ${pathout}/${samplename}.outsorted.vcf
