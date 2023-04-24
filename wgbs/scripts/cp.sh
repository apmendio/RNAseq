cd /share/lasallelab/Aron/WGBS_pws/wt_cr

/share/lasallelab/Aron/WGBS_pws/wt_cr/wgbs_lane1/Un_DTSA432/Project_JLAM_Nova393P_Mendiola
/share/lasallelab/Aron/WGBS_pws/wt_cr/wgbs_lane2/Un_DTSA466/Project_JLAM_Nova393P_Mendiola
/share/lasallelab/Aron/WGBS_pws/wt_cr/wgbs_lane3/Un_DTSA487/Project_JLAM_Nova393P_Mendiola
/share/lasallelab/Aron/WGBS_pws/wt_cr/wgbs_lane4/Un_DTSA524/Project_JLAM_Nova393_Mendiola

mkdir merge_lanes 
rsync -azPnv /share/lasallelab/Aron/WGBS_pws/wt_cr/wgbs_lane1/Un_DTSA432/Project_JLAM_Nova393P_Mendiola/*.fastq.gz /share/lasallelab/Aron/WGBS_pws/wt_cr/merge_lanes
rsync -azP /share/lasallelab/Aron/WGBS_pws/wt_cr/wgbs_lane1/Un_DTSA432/Project_JLAM_Nova393P_Mendiola/*.fastq.gz /share/lasallelab/Aron/WGBS_pws/wt_cr/merge_lanes

rsync -azPnv /share/lasallelab/Aron/WGBS_pws/wt_cr/wgbs_lane2/Un_DTSA466/Project_JLAM_Nova393P_Mendiola/*.fastq.gz /share/lasallelab/Aron/WGBS_pws/wt_cr/merge_lanes
rsync -azP /share/lasallelab/Aron/WGBS_pws/wt_cr/wgbs_lane2/Un_DTSA466/Project_JLAM_Nova393P_Mendiola/*.fastq.gz /share/lasallelab/Aron/WGBS_pws/wt_cr/merge_lanes

rsync -azPnv /share/lasallelab/Aron/WGBS_pws/wt_cr/wgbs_lane3/Un_DTSA487/Project_JLAM_Nova393P_Mendiola/*.fastq.gz /share/lasallelab/Aron/WGBS_pws/wt_cr/merge_lanes
rsync -azP /share/lasallelab/Aron/WGBS_pws/wt_cr/wgbs_lane3/Un_DTSA487/Project_JLAM_Nova393P_Mendiola/*.fastq.gz /share/lasallelab/Aron/WGBS_pws/wt_cr/merge_lanes

rsync -azPnv /share/lasallelab/Aron/WGBS_pws/wt_cr/wgbs_lane4/Un_DTSA524/Project_JLAM_Nova393_Mendiola/*.fastq.gz /share/lasallelab/Aron/WGBS_pws/wt_cr/merge_lanes
rsync -azP /share/lasallelab/Aron/WGBS_pws/wt_cr/wgbs_lane4/Un_DTSA524/Project_JLAM_Nova393_Mendiola/*.fastq.gz /share/lasallelab/Aron/WGBS_pws/wt_cr/merge_lanes

#renaming lanes

#trial (echo performs renaming without actually renaming the files, removing echo will perform the taskq)
for i in *L001*.fastq.gz do echo mv ${i/L001/L003} done

#actual renaming
for i in *L001*.fastq.gz do mv ${i/L001/L003} done

if md5sum -c \md5
then
    echo
    echo "All files have the correct md5sum"
else
    echo "ERROR: Some files are corrupt or missing"
    exit 1
fi

rsync -azPnv /share/lasallelab/Rochelle/RC_Prad_Circadian/RC_RNA_seq_circ/raw_sequences/htcounts ""

sbatch --array=1-12 /share/lasallelab/programs/CpG_Me/Paired-end/CpG_Me_PE_controller.sh mm10

sbatch --array=13-91 /share/lasallelab/programs/CpG_Me/Paired-end/CpG_Me_PE_controller.sh mm10
