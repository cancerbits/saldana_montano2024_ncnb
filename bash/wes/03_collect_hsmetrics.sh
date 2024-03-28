#!/bin/bash

# Usage
# from base dir run: bash bash/wes/03_collect_hsmetrics.sh metadata/Twist_Comprehensive_Exome_Covered_Targets_hg38.bed sarek_run/results/Preprocessing/*/DuplicatesMarked/*.bam


interval_file="$1"
bam_files=("${@:2}")
ref_base=igenomes_ref
ref_path=${ref_base}/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/

# Check if the provided file is in BED format and convert it to interval file if needed
if [[ "${interval_file: -4}" == ".bed" ]]; then
    echo "Converting BED to interval file..."

    interval_list="${interval_file%.bed}.interval_list"
    if [[ -f "$interval_list" ]]; then
        echo "interval_list file already exists. Please run the command again with $interval_list"
        exit 1
    else
        mkdir -p hsmetrics
        picard BedToIntervalList \
        I=${interval_file} \
        O="hsmetrics/$(basename ""${interval_file%.bed}.interval_list")" \
        SD=${ref_path}Homo_sapiens_assembly38.dict
        interval_list="hsmetrics/$(basename ""${interval_file%.bed}.interval_list")"
        echo "Done converting BED to interval_list"
    fi
elif [[ "${interval_file: -14}" == ".interval_list" ]]; then
    interval_list=$interval_file
else
    echo "Unrecognized file format. Please provide a BED or an interval_list file."
    exit 1
fi

# Loop through each BAM file and run Picard tool
mkdir -p hsmetrics
for bam_file in "${bam_files[@]}"; do
    echo "Processing $bam_file..."
    name="$(basename "${bam_file%.bam}")"

    picard CollectHsMetrics \
        -I ${bam_file} \
        -O hsmetrics/${name}.hs_metrics \
        -R ${ref_path}Homo_sapiens_assembly38.fasta \
        -BI ${interval_list} \
        -TI ${interval_list} \
        --PER_TARGET_COVERAGE hsmetrics/${name}.hs_target_cov \
        --THEORETICAL_SENSITIVITY_OUTPUT hsmetrics/${name}.hs_theor_sensitivity 2> hsmetrics/${name}_hsmetrics.log
    echo "CollectHsMetrics for $name.bam..."
done
