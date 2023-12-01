#!/bin/bash

# Collect STAR alignment stats from a given directory containing star ".Log.final.out" files

# Usage: ./collect_star_stats.sh path/to/star/log_files

dir=$1
outfile=$2

for i in $dir/*.Log.final.out; do

  filename=$i

  # get the base sample name from the filename
  # grep -l : returns filename
  # use sed to replace the suffix with nothing
  # use sed to replace the inpath directory with nothing
  # use sed to replace the final "/" with nothing
  grep -l "Number of input reads" $filename | sed -e 's/.Log.final.out//g' | sed -e "s|$dir||g" | sed -e 's/\///g' >> tmp_sample_name

  # Get Number of Input Reads
  grep "Number of input reads" $filename | cut -f2 >> tmp_input_reads

  # Get Uniquely Mapped Reads Number
  grep "Uniquely mapped reads number" $filename | cut -f2 >> tmp_unq_number

  # Get Uniquely Mapped Reads Percent
  grep "Uniquely mapped reads %" $filename | cut -f2 >> tmp_unq_percent

  # Get Multimapped Reads Number
  grep "Number of reads mapped to multiple loci" $filename | cut -f2 >> tmp_multi_number

  # Get Multimapped Reads Percent
  grep "% of reads mapped to multiple loci" $filename | cut -f2 >> tmp_multi_percent

  # Get Reads Mapped to Too Many Loci Number
  grep "Number of reads mapped to too many loci" $filename | cut -f2 >> tmp_too_many_loci_number

  # Get Reads Mapped to Too Many Loci Percent
  grep "% of reads mapped to too many loci" $filename | cut -f2 >> tmp_too_many_loci_percent

  # Get Percent Reads Unmapped: Too Many Mismatches
  grep " % of reads unmapped: too many mismatches" $filename | cut -f2 >> tmp_too_many_mismatch_percent

  # Get Percent Reads Unmapped: Too Short
  grep "% of reads unmapped: too short" $filename | cut -f2 >> tmp_too_short_percent

  # Get Percent Reads Unmapped: Other
  grep "% of reads unmapped: other" $filename | cut -f2 >> tmp_unmapped_other_percent;
done

paste -d "\t" tmp_sample_name tmp_input_reads tmp_unq_number tmp_unq_percent tmp_multi_number tmp_multi_percent tmp_too_many_loci_number tmp_too_many_loci_percent tmp_too_many_mismatch_percent tmp_too_short_percent tmp_unmapped_other_percent > tmp_out.txt

# sort on the first column (sample name)
sort -k1,1 tmp_out.txt > $outfile

sed -i '1 i\SampleName\tInputReads\tUnqMap\tUnqMapPerc\tMultiMap\tMultiMapPerc\tTooManyLoci\tTooManyLociPerc\tUnmappedTooManyMismatch\tUnmappedTooShort\tUnmappedOther' $outfile

rm tmp_*
