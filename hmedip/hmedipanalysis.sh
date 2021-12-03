raw_sequences

module load trim_galore/0.6.6
source activate cutadapt-2.10

mkdir 01-TRIMGALORE2
for i in `cat task_samples.txt`
do 
	trim_galore \
	--paired \
	--basename ${i}_trimmed \
	--cores 4 \
	--fastqc \
	--output_dir 01-TRIMGALORE2 \
	raw_sequences/${i}*1*.fq.gz \
	raw_sequences/${i}*2*.fq.gz
done

#files that didnt trim
for i in `cat task_samples2.txt`
do 
	trim_galore \
	--paired \
	--basename ${i}_trimmed \
	--cores 4 \
	--fastqc \
	--output_dir 01.1TRIMGALORE \
	raw_sequences/${i}_1.fq.gz \
	raw_sequences/${i}_2.fq.gz
done

module load bowtie2

mkdir 02-ALIGNMENT
for i in `cat task_samples.txt`
do
	bowtie2 \
	-t \
	-x /share/lasallelab/genomes/mm10/Bowtie2 \
	-1 01-TRIMGALORE2/${i}*R1*trimmed.fq.gz \
	-2 01-TRIMGALORE2/${i}*R2*trimmed.fq.gz \
	-S 02-ALIGNMENT/${i}_aligned.sam
done

module load trim_galore/0.6.6
source activate cutadapt-2.10

mkdir 01-TEST
for i in `cat test_samples.txt`
do 
	parallel --plus --dry-run --jobs 4 --will-cite
	'trim_galore \
	--paired \
	--basename ${i}_trimmed \
	--cores 4 \
	--fastqc \
	--output_dir 01-TRIMGALORE/${i} \
	TEST/${i}*1*.fq.gz \
	TEST/${i}*2*.fq.gz' \
	::: task_samples.txt
done

mkdir 02-TEST
for i in `cat test_samples.txt`
do
	bowtie2 \
	-t \
	-x /share/lasallelab/genomes/mm10/Bowtie2/mm10 \
	-1 test/${i}_trimmed*1*.fq.gz \
	-2 test/${i}_trimmed*2*.fq.gz \
	-S 02-TEST/${i}_aligned.sam
done

mkdir 02-TEST
for i in `cat test_samples.txt`
do
	bowtie2 \
	-t \
	-x /share/lasallelab/genomes/mm10/Bowtie2/mm10 \
	-1 test/${i}_trimmed*1*.fq.gz \
	-2 test/${i}_trimmed*2*.fq.gz \
	-S 02-TEST/${i}_aligned.sam
done
}


parallel -j 2 "bowtie2 -t --threads 4 -x /share/lasallelab/genomes/mm10/Bowtie2/mm10 -1 {1} -2 {2} -S 02-TEST/{1/.}_aligned.sam" ::: test/*1*.fq.gz :::+ test/*2*.fq.gz

#bowtie 2 parallelization
mkdir 02-ALIGNMENT
parallel -j 4 "bowtie2 -t --threads 4 -x /share/lasallelab/genomes/mm10/Bowtie2/mm10 -1 {1} -2 {2} -S 02-ALIGNMENT/{1/}_aligned.sam" ::: 01-TRIMGALORE2/*trimmed*1*.fq.gz :::+ 01-TRIMGALORE2/*trimmed*2*.fq.gz

parallel -j 4 "bowtie2 -t --threads 4 -x /share/lasallelab/genomes/mm10/Bowtie2/mm10 -1 {1} -2 {2} -S 02-ALIGNMENT/{1/}_aligned.sam" ::: 01.1TRIMGALORE/*trimmed*1*.fq.gz :::+ 01.1TRIMGALORE/*trimmed*2*.fq.gz

module load samtools

#conver SAM to BAM
mkdir 02.5-BAM
for i in `cat task_samples.txt`
do
	samtools view -S -b 02-ALIGNMENT/${i}*.sam > 02.5-BAM/${i}.bam
done
mkdir 03-sortedBAM
parallel -j 4 "samtools sort -@ 4 {1} -o 03-sortedBAM/{1/.}_sorted.bam" ::: 01.1TRIMGALORE/*trimmed*1*.fq.gz :::+ 01.1TRIMGALORE/*trimmed*2*.fq.gz