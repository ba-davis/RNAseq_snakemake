#!/bin/bash

fqc=$1
trim=$2
aln=$3
outfile=$4

join -t $'\t' -1 1 -2 1 <(sort -t$'\t' -k1,1 $fqc) <(sort -t$'\t' -k1,1 $trim) | join -t $'\t' -1 1 -2 1 - <(sort -t$'\t' -k1,1 $aln) > $outfile

# add a header
sed -i '1s/^/Sample_Name\tR1_raw_reads\tR2_raw_reads\tR1_read_length\tR2_read_length\tR1_GC%\tR2_GC%\tR1_dup%\tR2_dup%\tInputReadPairs\tBothSurviving\tBothSurvivingPerc\tForwardSurviving\tForwardSurvivingPerc\tReverseSurviving\tReverseSurvivingPerc\tReadPairsDropped\tReadPairsDroppedPerc\tInputReads\tUnqMap\tUnqMapPerc\tMultiMap\tMultiMapPerc\tTooManyLoci\tTooManyLociPerc\tUnmappedTooManyMismatch\tUnmappedTooShort\tUnmappedOther\n/' $outfile
