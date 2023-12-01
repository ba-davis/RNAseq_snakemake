#!/bin/bash

# PE log files

# Collect Trimmomatic stats from a given directory containing trimmomatic log files
#  these are often the ".err" files after a SLURM submission

# Usage:

dir=$1            # path to Trimmomatic log files
raw_file_path=$2  # path to raw fastq files
suffix=$3         # ending of fastq files to remove, leaving sample name
outfile=$4        # summary txt outfile

for i in $dir/trimmomatic_*.err; do

    filename=$i

    # Get Sample Name
    grep "phred" $filename | cut -d " " -f3 | sed -e "s|$raw_file_path||g" | sed -e 's/\///g' | sed -e "s/$suffix//g" >> tmp_sample_name

    # Get Input Read Pairs
    grep "Input Read Pairs" $filename | cut -d " " -f4 >> tmp_input_read_pairs
    
    # Get Both Surviving
    grep "Input Read Pairs" $filename | cut -d " " -f7 >> tmp_both_surviving

    # Get Both Surviving Percent
    grep "Input Read Pairs" $filename | cut -d " " -f8 | tr -d '()' >> tmp_both_surviving_percent

    # Get Forward Only Surviving
    grep "Input Read Pairs" $filename | cut -d " " -f12 >> tmp_forward_surviving

    # Get Forward Only Surviving Percent
    grep "Input Read Pairs" $filename | cut -d " " -f13 | tr -d '()' >> tmp_forward_surviving_percent

    # Get Reverse Only Surviving
    grep "Input Read Pairs" $filename | cut -d " " -f17 >> tmp_reverse_surviving

    # Get Reverse Only Surviving Percent
    grep "Input Read Pairs" $filename | cut -d " " -f18 | tr -d '()'>> tmp_reverse_surviving_percent

    # Get Dropped
    grep "Input Read Pairs" $filename | cut -d " " -f20 >> tmp_drop

    # Get Dropped Percent
    grep "Input Read Pairs" $filename | cut -d " " -f21 | tr -d '()' >> tmp_drop_percent;
done

paste -d "\t" tmp_sample_name tmp_input_read_pairs tmp_both_surviving tmp_both_surviving_percent tmp_forward_surviving tmp_forward_surviving_percent tmp_reverse_surviving tmp_reverse_surviving_percent tmp_drop tmp_drop_percent > tmp_out.txt

# sort on the first column (sample name)
sort -k1,1 tmp_out.txt > $outfile

sed -i '1 i\SampleName\tInputReadPairs\tBothSurviving\tBothSurvivingPerc\tForwardSurviving\tForwardSurvivingPerc\tReverseSurviving\tReverseSurvivingPerc\tReadPairsDropped\tReadPairsDroppedPerc' $outfile

rm tmp_*
