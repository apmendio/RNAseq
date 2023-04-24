#Trimgalore (does trimming and fastqc)
#Fastqscreen (alignment quant)
#Bowtie2 (alignment)
#Phantomqualpeaks (quantifying peaks)
#MACS2 (peak calling)
##	-S 02-BOWTIE2/${i}/${i}.sam \ used for the sam file output, not required when using a pipe to generate a bam
# | - pipes allow you to remove intermediate files and variables
# know the resources to identify which parts can be piped

module load trim_galore/0.6.6
source activate cutadapt-2.10

trimChIP(){
	i=$1
	mkdir 01-TRIMGALORE
	trim_galore \
	--paired \
	--illumina \
	--basename ${i} \
	--cores 8 \
	--fastqc \
	--output_dir 01-TRIMGALORE/ \
	/share/lasallelab/Aron/wgbs-hmedip/hmedip/hmedip_1.30.22/raw_sequences/${i}_1.fq.gz \
	/share/lasallelab/Aron/wgbs-hmedip/hmedip/hmedip_1.30.22/raw_sequences/${i}_2.fq.gz 
}
export -f trimChIP
cat task_samples.txt | parallel --will-cite --jobs 8 trimChIP

module load bowtie2
module load samtools

alignChIP(){
	i=$1
	mkdir 02-BOWTIE2
	(bowtie2 \
	-x /share/lasallelab/genomes/mm10/Bowtie2/mm10 \
	-1 01-TRIMGALORE/${i}_val_1.fq.gz \
	-2 01-TRIMGALORE/${i}_val_2.fq.gz \
#	--no-unal \
	-p 20) \
	2>02-BOWTIE2/${i}.log | \
	samtools view -bS - > \
	02-BOWTIE2/${i}.bam 
}
export -f alignChIP
cat task_samples.txt | parallel --will-cite --jobs 8 alignChIP

module load picard-tools
module load samtools

dedupChIP(){
	i=$1
	mkdir 03-PICARD
	picard SortSam \
	INPUT=02-BOWTIE2/${i}.bam \
	OUTPUT=03-PICARD/${i}_sorted.bam \
	SORT_ORDER=coordinate
	
	picard MarkDuplicates \
	I=03-PICARD/${i}_sorted.bam \
	O=03-PICARD/${i}_sorted_marked_duplicates.bam \
	M=03-PICARD/${i}_marked_dup_metrics.txt \
	REMOVE_DUPLICATES=true \
	ASSUME_SORT_ORDER=coordinate
	
	picard CollectInsertSizeMetrics \
	INPUT=03-PICARD/${i}_sorted_marked_duplicates.bam \
	OUTPUT=03-PICARD/${i}.insert.txt \
	HISTOGRAM_FILE=03-PICARD/${i}.insert.histogram.pdf \
	ASSUME_SORTED=TRUE
}
export -f dedupChIP
cat task_samples.txt | parallel --will-cite --jobs 8 dedupChIP

module load picard-tools
mkdir 03-PICARD_dedup
for i in `cat task_samples.txt`
do 
	picard SortSam \
	INPUT=02-BOWTIE2/${i}.bam \
	OUTPUT=03-PICARD_dedup/${i}_sorted.bam \
	SORT_ORDER=coordinate
	
	picard MarkDuplicates \
	I=03-PICARD_dedup/${i}_sorted.bam \
	O=03-PICARD_dedup/${i}_sorted_marked_duplicates.bam \
	M=03-PICARD_dedup/${i}_marked_dup_metrics.txt \
	REMOVE_DUPLICATES=true \
	ASSUME_SORT_ORDER=coordinate
	
	picard CollectInsertSizeMetrics \
	INPUT=03-PICARD_dedup/${i}_sorted_marked_duplicates.bam \
	OUTPUT=03-PICARD_dedup/${i}.insert.txt \
	HISTOGRAM_FILE=03-PICARD_dedup/${i}.insert.histogram.pdf \
	ASSUME_SORTED=TRUE
done

indexChIP(){
	i=$1
	samtools index \
	03-PICARD/${i}_sorted.bam
}
export -f indexChIP
cat task_samples.txt | parallel --will-cite --jobs 8 indexChIP

module load macs2

callPeaks(){
	i=$1
	mkdir 04.3-PEAKCALL
	macs2 callpeak \
	-t 03-PICARD/${i}_sorted.bam \
	-f BAMPE \
	-B \
	--nomodel \
#	--extsize  \
#	--bw 100 \
#	--SPMR \
	-n ${i} \
	-g mm \
	--cutoff-analysis \
	--outdir 04.3-PEAKCALL
}
export -f callPeaks
cat task_samples.txt | parallel --will-cite --jobs 8 callPeaks

callPeaks(){
	i=$1
	mkdir 04.1-PEAKCALL
	macs3 callpeak \
	-t 03-PICARD/${i}_sorted.bam \
	-f BAMPE \
	-B \
	--nomodel \
	-n ${i} \
	-g mm \
	--outdir 04.1-PEAKCALL/
}
export -f callPeaks
cat task_samples.txt | parallel --will-cite --jobs 8 callPeaks

mkdir 05-DEEPTOOLS
module load anaconda3
module load deeptools

conBigWig(){
	i=$1
	mkdir 05-DEEPTOOLS
	bamCoverage \
	-b  03-PICARD/${i}_sorted_marked_duplicates.bam \
	-o 05-DEEPTOOLS/${i}_coverage.bw
}
export -f conBigWig
cat task_males.txt | parallel --will-cite --jobs 8 conBigWig

#compare samples to input control
compareBam(){
	i=$1
	bamCompare \
	-
}

rsync -azPnvr apmendio@epigenerate.genomecenter.ucdavis.edu:/share/lasallelab/Aron/wgbs-hmedip/hmedip/03-sortedBAM "/Users/aron/Downloads"

rsync -azPnvr apmendio@epigenerate.genomecenter.ucdavis.edu:/share/lasallelab/Aron/wgbs-hmedip/hmedip/05-PEAKCALL "/Users/aron/Downloads"

rsync -azPnvr apmendio@epigenerate.genomecenter.ucdavis.edu:/share/lasallelab/Aron/wgbs-hmedip/hmedip/05.2-PEAKCALL "/Users/aron/Downloads"

rsync -azP blaufer@epigenerate.genomecenter.ucdavis.edu:/share/lasallelab/Archive "/Volumes/lasalle lab/Sequence Data Archive"

rsync -azP apmendio@epigenerate.genomecenter.ucdavis.edu:/share/lasallelab/Aron/RNA_seq_analysis/Analysis_Aug12022/multiqc_report.html "/Users/aron/Downloads"

