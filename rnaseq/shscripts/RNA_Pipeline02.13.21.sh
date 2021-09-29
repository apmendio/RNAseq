#RNA_seq_Pipeline

#Navigate to main directory
#Make directories for samples
#Transfer or copy files into sample directories 
#A samples are male and B samples are female
#Trim
#Align
#Map
#Quantify

#Makes folders for raw files
cd /share/lasallelab/Aron/RNA_seq_analysis/Analysis1_Feb112021
mkdir 00-RawData

for i in `cat task_samples.txt`
do 
	mkdir 00-RawData/${i}
done

#testing file transfer
#ran trial and it seems like it's working 

for i in `cat task_samples.txt`
do
rsync -azPnv /share/lasallelab/Aron/RNA_seq/Dec222020/52.7.5.48/1735-880-12212020/${i}_*.fastq.gz /share/lasallelab/Aron/RNA_seq_analysis/Analysis1_Feb112021/00-RawData/${i}
done

#File Transfer
	
for i in `cat task_samples.txt`
do
rsync -azP /share/lasallelab/Aron/RNA_seq/Dec222020/52.7.5.48/1735-880-12212020/${i}_*.fastq.gz /share/lasallelab/Aron/RNA_seq_analysis/Analysis1_Feb112021/00-RawData/${i}
done

#Trimming using trim_galore
#will try --hardtrim and --quality to see differences in sequence quality
#cut adapt needs to be loaded this way in davis nodes
mkdir 01-TRIMGALORE
module load trim_galore/0.6.6
source activate cutadapt-2.10

#Trimming
# removing adapters

mkdir 01-TRIMGALORE
for i in `cat trial_sampleset.txt`
do 
	mkdir 01-TRIMGALORE/${i}
	trim_galore \
	--paired \
	--basename ${i}_trimmed \
	--cores 8 \
	--fastqc \
	--output_dir 01-TRIMGALORE/${i} \
	00-RawData/${i}/${i}*R1*.fastq.gz \
	00-RawData/${i}/${i}*R2*.fastq.gz
done

#Trim did ont produce fastqc files although --fastqc was specified
#Need to run fastqc

for i in `cat trial_sampleset.txt`
do
	fastqc \
	01-TRIMGALORE/${i}/${i}*R1*.fq.gz
done
#Building Star index
	STAR \
	--runMode genomeGenerate \
	--runThreadN 25 \
	--genomeDir /share/lasallelab/Aron/mouserefgenome/mm10/star \
	--genomeFastaFiles /share/lasallelab/Aron/mouserefgenome/mm10/genomes/GRCm38.primary_assembly.genome.fa \
	--sjdbGTFfile /share/lasallelab/Aron/mouserefgenome/mm10/genomes/gencode.vM25.primary_assembly.annotation.gtf
#STAR Alignmient
mkdir 02-alignment

for i in `cat trial_sampleset.txt` 
do 
	mkdir 02-alignment/${i}
	STAR \
	--runMode alignReads \
	--runThreadN 25 \
	--outSAMtype BAM \
	Unsorted --readFilesCommand zcat \
	--genomeDir /share/lasallelab/Aron/mouserefgenome/mm10/star \
	--outFileNamePrefix 02-alignment/${i}/${i}  \
	--readFilesIn 01-TRIMGALORE/${i}/${i}*1*.fq.gz 01-TRIMGALORE/${i}/${i}*2*.fq.gz
done

#Salmon quantification
mkdir 03-quantification

for i in `cat trial_sampleset.txt`
do
	salmon \
	quant \
	-i /share/lasallelab/Aron/mouserefgenome/mm10/CDS_GRCm39/salmon_index \
	-p 24 \
	-l A \
	-1 /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/01-TRIMGALORE/${i}/${i}*1*.fq.gz \
	-2 /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/01-TRIMGALORE/${i}/${i}*2*.fq.gz \
	--validateMappings \
	-o /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/03-quantification/${i}
done 
	

#Does not work, may be due to incorrect index
mkdir 03.1-bamquant

for i in `cat trial_sampleset.txt`
do	
	salmon \
	quant \
	-t /share/lasallelab/Aron/mouserefgenome/mm10/CDS_GRCm39/gentrome.fa.gz \
	-l A \
	-a /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/02-alignment/${i}/${i}Aligned.out.bam \
	-o /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/03.1-bamquant/${i}
done

#featureCounts - before using reads must first be sorted using samtools
mkdir 02.5-alignsort
for i in `cat trial_sampleset.txt`
do
mkdir 02.5-alignsort/${i}
	samtools sort \
	-m 2G \
	-@ 30 \
	-O BAM \
	-o 02.5-alignsort/${i}/${i}_sorted.bam \
	02-alignment/${i}/${i}Aligned.out.bam
done


	
mkdir 04-featurecounts
for i in `cat trial_sampleset.txt`
do
mkdir 04-featurecounts/${i}
featureCounts -T 12 -s 1 -p \
  -a /share/lasallelab/Aron/mouserefgenome/mm10/genomes/gencode.vM25.primary_assembly.annotation.gtf \
  -o /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/04-featurecounts/${i}/${i}_Ffeatures \
  02.5-alignsort/${i}/${i}_sorted.bam
done

for i in `cat trial_sampleset.txt`
do
mkdir 04-featurecounts/R${i}
featureCounts -T 12 -s 2 -p \
  -a /share/lasallelab/Aron/mouserefgenome/mm10/genomes/gencode.vM25.primary_assembly.annotation.gtf \
  -o /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/04-featurecounts/R${i}/${i}_Rfeatures \
  02.5-alignsort/${i}/${i}_sorted.bam
done
#use this for featurecounts 1 is forward 2 is reverse

mkdir 04-featurecounts

for i in `cat trial_sampleset.txt`
do
mkdir 04-featurecounts/${i}
featureCounts -T 12 -s 1 -p -P -d 50 -D 200 \
  -a /share/lasallelab/Aron/mouserefgenome/mm10/genomes/gencode.vM25.primary_assembly.annotation.gtf \
  -o /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/04-featurecounts/${i}/${i}counts.txt \
  02.5-alignsort/${i}/${i}_sorted.bam
done

for i in `cat trial_sampleset.txt`
do
mkdir 04-featurecounts/R${i}
featureCounts -T 12 -s 1 -p -P -d 50 -D 200 \
  -a /share/lasallelab/Aron/mouserefgenome/mm10/genomes/gencode.vM25.primary_assembly.annotation.gtf \
  -o /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/04-featurecounts/R${i}/R${i}counts.txt \
  02.5-alignsort/${i}/${i}_sorted.bam
done

featureCounts -T 12 -s 1 -p -P -d 100 -D 200 \
  -a /share/lasallelab/Aron/mouserefgenome/mm10/genomes/gencode.vM25.primary_assembly.annotation.gtf \
  -o /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/04-featurecounts/100counts.txt \
  02.5-alignsort/*/*.bam
  
featureCounts -T 12 -s 2 -p -P -d 100 -D 200 \
  -a /share/lasallelab/Aron/mouserefgenome/mm10/genomes/gencode.vM25.primary_assembly.annotation.gtf \
  -o /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/04-featurecounts/R100counts.txt \
  02.5-alignsort/*/*.bam
  
for i in `cat trial_sampleset.txt`
do 
mv 02.5-alignsort/${i}/*bam 02.5-alignsort/allsortedbam
done
