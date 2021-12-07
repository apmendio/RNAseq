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
	--cores 8 \
	--fastqc \
	--output_dir 01-TRIMGALORE/ \
	/share/lasallelab/Aron/wgbs-hmedip/hmedip/Un_DTSA432/Project_JLAM_Hm5DIP002/${i}*_R1_*.fastq.gz \
	/share/lasallelab/Aron/wgbs-hmedip/hmedip/Un_DTSA432/Project_JLAM_Hm5DIP002/${i}*_R2_*.fastq.gz 
}
export -f trimChIP
cat task_males.txt | parallel --will-cite --jobs 8 trimChIP

module load bowtie2
module load samtools

alignChIP(){
	i=$1
	mkdir 02-BOWTIE2
	(bowtie2 \
	-x /share/lasallelab/genomes/mm10/Bowtie2/mm10 \
	-1 01-TRIMGALORE/${i}*_R1_*_1.fq.gz \
	-2 01-TRIMGALORE/${i}*_R2_*_2.fq.gz \
#	--no-unal \
	-p 20) \
	2>02-BOWTIE2/${i}.log | \
	samtools view -bS - > \
	02-BOWTIE2/${i}.bam 
}
export -f alignChIP
cat task_males.txt | parallel --will-cite --jobs 8 alignChIP

#conver SAM to BAM
mkdir 02.5-BAM
for i in `cat task_samples2.txt`
do
	samtools view -S -b 02-ALIGNMENT/${i}_trimmed_val_1.fq.gz_aligned.sam > 02.5-BAM/${i}.bam
done
mkdir 03-sortedBAM

parallel -j 6 "samtools sort -@ 4 {1} -o 03-sortedBAM/{1/.}_sorted.bam" ::: 02.5-BAM/*bam

module load picard-tools
mkdir 04-PICARD

parallel -j 6 "picard MarkDuplicates I={1} O=04-PICARD/{1/.}.markeddupes.bam M={1/.}_markeddupes.txt" ::: 03-sortedBAM/*.bam
parallel -j 6 "picard MarkDuplicates --I={1} O=04-PICARD/{1/.}.markeddupes.bam M={1/.}_markeddupes.txt" ::: 03-sortedBAM/*.bam

module load macs2
mkdir 05-PEAKCALL
callPeaks(){
	i=$1
	macs2 callpeak \
	-t 04-PICARD/${i}_sorted.markeddupes.bam \
	-f BAMPE \
	-B \
	--nomodel \
	--extsize 200 \
	--bw 150 \
	--SPMR \
	-n ${i} \
	-g mm \
	--cutoff-analysis \
	--outdir 05.2-PEAKCALL/
}
export -f callPeaks
cat task_males.txt | parallel --will-cite --jobs 6 callPeaks

mkdir 05-DEEPTOOLS
module load anaconda3
module load deeptools

conBigWig(){
	i=$1
	bamCoverage \
	-b  04-PICARD/${i}_sorted.markeddupes.bam \
	-o 05-DEEPTOOLS/${i}_coverage.bw
}
export -f conBigWig
cat task_males.txt | parallel --will-cite --jobs 4 conBigWig

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

