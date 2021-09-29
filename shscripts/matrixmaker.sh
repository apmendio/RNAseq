# Annotation File
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz

# first column is the transcript ID, second column is the gene ID, third column is the gene symbol
zcat gencode.vM25.annotation.gtf.gz | awk -F "\t" 'BEGIN{OFS="\t"}{if($3=="transcript"){split($9, a, "\""); print a[4],a[2],a[8]}}' > tx2gene.gencode.vm25.csv

gsub("\\..*","", tx2gene.edit.gencode.vm25.csv)

wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.basic.annotation.gtf.gz
zcat gencode.vM25.basic.annotation.gtf.gz | awk -F "\t" 'BEGIN{OFS="\t"}{if($3=="transcript"){split($9, a, "\""); print a[4],a[2],a[8]}}' > tx2.1gene.gencode.vm25.csv

wget http://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
zcat Mus_musculus.GRCm38.102.gtf.gz | awk -F "\t" 'BEGIN{OFS="\t"}{if($3=="transcript"){split($9, a, "\""); print a[4],a[2],a[8]}}' > tx2gene.ensembl.102.csv

# UCSC
wget rsync://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz
zcat mm10.refGene.gtf.gz | awk -F "\t" 'BEGIN{OFS="\t"}{if($3=="transcript"){split($9, a, "\""); print a[4],a[2],a[8]}}' > tx2gene.refGeve.csv
### Identify Counts per read per sample 

grep -v "N_" A1_ReadsPerGene.out.tab | awk '{unst+=$2;forw+=$3;rev+=$4}END{print unst,forw,rev}'

### Generates total reads per sample
### Column 1 Total
### Column 2 Forward Strand
### Column 3 Reverse Strand

for i in `cat task_samples.txt`
do grep -v "N_" ${i}_ReadsPerGene.out.tab | awk '{unst+=$2;forw+=$3;rev+=$4}END{print unst,forw,rev}' > ${i}_counts_total.txt
done

### Generates count information with gene_id

for i in `cat task_samples.txt`
do grep -v "N_" ${i}_ReadsPerGene.out.tab > ${i}_Raw_gene_counts.txt
done

# retrieve the 4th column of each "ReadsPerGene.out.tab" file + the first column that contains the gene IDs
paste *_ReadsPerGene.out.tab | grep -v "_" | awk '{printf "%s\t", $1}{for (i=4;i<=NF;i+=4) printf "%s\t", $i; printf "\n" }' > tmp


for i in `cat task_samples.txt`
do paste ${i}_ReadsPerGene.out.tab | grep -v "_" | awk '{printf "%s\t", $1}{for (i=4;i<=NF;i+=4) printf "%s\t", $i; printf "\n" }' > tmp
done

# add header: "gene_name" + the name of each of the counts file
sed -e "1igene_name\t$(ls *_ReadsPerGene.out.tab | tr '\n' '\t' | sed 's/ReadsPerGene.out.tab//g')" tmp | cut -f1-7 > reverse_strand_counts.txt

for i in `cat task_samples.txt`
do sed -e "1iGene_ID\t$(ls ${i}_ReadsPerGene.out.tab | tr '\n' '\t' | sed 's/ReadsPerGene.out.tab//g')" tmp | cut -f1-7 > reverse_strand_counts.txt
done

# extract reverse reads

for i in *ReadsPerGene.out.tab
do echo $i
# retrieve the first (gene name) and fourth column (raw reads)
cut -f1,4 $i | grep -v "_" > `basename $i ReadsPerGene.out.tab`_counts.txt
done

# extract forward reads

for i in *ReadsPerGene.out.tab
do echo $i
# retrieve the first (gene name) and fourth column (raw reads)
cut -f1,3 $i | grep -v "_" > `basename $i ReadsPerGene.out.tab`_Fcounts.txt
done