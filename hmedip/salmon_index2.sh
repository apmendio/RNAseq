#generating genome from GRCm38
grep "^>" <(gunzip -c GRCm38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat gencode.vM25.transcripts.fa.gz GRCm38.primary_assembly.genome.fa.gz > gentrome.fa.gz
salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode 

grep -v "N_" /share/lasallelab/Aron/RNA_seq_analysis/trial_analysis/02-alignment/A1/A1SJ.out.tab | awk '{unst+=$2;forw+=$3;rev+=$4}END{print unst,forw,rev}'

grep "^>" <(gunzip -c GRCm38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat gencode.vM25.transcripts.fa.gz GRCm38.primary_assembly.genome.fa.gz > gentrome.fa.gz
salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode 

mkdir 01-TRIMGALORE
for i in `cat task_samples.txt`
do 
	mkdir 01-TRIMGALORE/${i}
	parallel --dry-run --jobs 4 --will-cite
	'trim_galore \
	--paired \
	--basename ${i}_trimmed \
	--cores 4 \
	--fastqc \
	--output_dir 01-TRIMGALORE/${i} \
	00-Raw/${i}/${i}*1*.fq.gz \
	00-Raw/${i}/${i}*2*.fq.gz' \
	::: task_samples.txt
done