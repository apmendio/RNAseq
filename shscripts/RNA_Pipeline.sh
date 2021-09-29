#RNA_seq_Pipeline

#Navigate to main directory
#Make directories for samples
#Transfer or copy files into sample directories 
# A samples are male and B samples are female

cd /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis
mkdir 00-RawData

for i in `cat trial_sampleset.txt`
do 
	mkdir 00-Raw/${i}
done

#testing file transfer
#ran trial and it seems like it's working 

for i in `cat trial_sampleset.txt`
do
rsync -azPnv /share/lasallelab/Aron/RNA_seq/Dec12020/52.7.5.48/${i}_*.fastq.gz /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/00-Raw/${i}
done
	
#will proceed with a smaller test set A1 and A7
#test run worked and transferred samples to appropriate folder

for i in `cat trial.txt`
do
rsync -azP /share/lasallelab/Aron/RNA_seq/Dec12020/52.7.5.48/${i}_*.fastq.gz /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/00-Raw/${i}
done

#actual run
for i in `cat trial_sampleset.txt`
do
rsync -azP /share/lasallelab/Aron/RNA_seq/Dec12020/52.7.5.48/${i}_*.fastq.gz /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/00-Raw/${i}
done

#Renaming samples in each folder for easier analysis
for i in `cat Read.txt`, k in `cat sample.txt` 
do
mv ${i} ${k}_R1.fastq.gz
done

#Trimming using trim_galore
#will try --hardtrim and --quality to see differences in sequence quality
mkdir 01-TRIMGALORE
module load trim_galore
module load cutadapt

# removing corrupted 50 bp reads 
# this step is not required for the actual samples only used for Dec122020 sequences

for i in `cat trial_sampleset.txt`
do 
	mkdir 01-TRIMGALORE/${i}
	trim_galore \
	--paired \
	--hardtrim5 100 \
	--cores 8 \
	--fastqc \
	--output_dir 01-TRIMGALORE/${i} \
	00-RawData/${i}/${i}*R1*.fastq.gz \
	00-RawData/${i}/${i}*R2*.fastq.gz 
done

# removing adapters
for i in `cat trial.txt`
do 
	trim_galore \
	--paired \
	--cores 8 \
	--fastqc \
	--output_dir 01-TRIMGALORE/${i} \
	01-TRIMGALORE/${i}/${i}*R1*.fq.gz \
	01-TRIMGALORE/${i}/${i}*R2*.fq.gz 
done

#Trim seemed successful

for i in `cat trial_sampleset.txt`
do 
	mkdir 01-TRIMGALORE/${i}
	trim_galore \
	--paired \
	--quality 100 \
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
	--readFilesIn 01-TRIMGALORE/${i}/${i}*R1*.fq.gz 01-TRIMGALORE/${i}/${i}*R2*.fq.gz
done

#Salmon quantification
for i in `cat trial_sampleset.txt`
do
	salmon \
	quant \
	-i /share/lasallelab/Aron/mouserefgenome/mm10/CDS_GRCm39/salmon_index \
	-p 24 \
	-l A \
	-1 /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/01-TRIMGALORE/${i}/${i}*R1*.fq.gz \
	-2 /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/01-TRIMGALORE/${i}/${i}*R2*.fq.gz \
	--validateMappings \
	-o /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/03-quantification/${i}
done 
	

#Does not work, may be due to incorrect index
	
salmon quant -t /share/lasallelab/Aron/mouserefgenome/mm10/CDS_GRCm39/gentrome.fa.gz -l A -a /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/02-alignment/A1/A1Aligned.out.bam -o /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/03-quantification/A1/salmon_quant 

#featureCounts

mkdir 04-featurecounts
for i in `cat trial_sampleset.txt`
do
mkdir 04-featurecounts/${i}
featureCounts -T 12 -s 2 -p \
  -a /share/lasallelab/Aron/mouserefgenome/mm10/genomes/gencode.vM25.primary_assembly.annotation.gtf \
  -o /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/04-featurecounts/${i}/${i}_features_Rev \
  /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/02-alignment/${i}/${i}Aligned.out.bam
done
