# Load the sequenza_conda.yaml environment

# Paths
hg38_ref_path = "" # Path to the reference genome FASTA
ref_dir_path  = "" # Path to GC content reference
bam_loc_path  = "" # Location of the bam files
results_path  = "" # Where you'd like the results to go

# calculate GC content (once)
sequenza-utils gc_wiggle -f $hg38_ref_path \
-o ${ref_dir_path}/Homo_sapiens_assembly38_gc_wiggle_sequenza_50bp.txt.gz

# Shell script for generating Sequenza input
for sample in sample_name1 sample_name2 sample_name3
do
	# create seqz file
	sequenza-utils bam2seqz -n ${bam_loc_path}/H7S14.recal.bam \ # this is the parental bam file
	-t ${bam_loc_path}/${sample}.recal.bam \
	-gc ${ref_dir_path}/Homo_sapiens_assembly38_gc_wiggle_sequenza_50bp.txt.gz \
	-F $hg38_ref_path \
	-o ${results_path}/NCNB_${sample}_parental_normal.seqz.gz --het 0.4

	# bin the data
	sequenza-utils seqz_binning -s ${results_path}/NCNB_${sample}_parental_normal.seqz.gz \
	-o ${results_path}/NCNB_${sample}_parental_normal_binned.seqz.gz
done
