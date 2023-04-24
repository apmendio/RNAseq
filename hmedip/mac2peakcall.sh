module load macs2
mkdir 04-PEAKCALL
callPeaks(){
	i=$1
	mkdir 04-PEAKCALL_updated/${i}
	macs2 callpeak \
	-t 03-PICARD/${i}/${i}_sorted.bam \
	-c 03-PICARD/14_inputcontrol/14_inputcontrol_sorted.bam \
	-f BAMPE \
	-B \
	--nomodel \
	--exitsize 150 \
	--SPMR \
	-n ${i} \
	-g mm \
	--cutoff-analysis \
	--outdir 04-PEAKCALL_updated/${i}
}
export -f callPeaks
cat task_samples.txt | parallel --will-cite --jobs 4 callPeaks

module load sratoolkit
fetch(){
	i=$1
	prefetch -v ${i}
}
export -f fetch
cat task_samples.txt | parallel --will-cite --jobs 9 fetch

dump(){
	i=$1
	fastq-dump -O . --split-files ${i}/${i}.sra
}
export -f dump
cat task_samples.txt | parallel --will-cite --jobs 9 dump