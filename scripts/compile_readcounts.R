#!/usr/bin/env Rscript

source("scripts/downstream_functions.R")

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
strand <- as.numeric(args[2])
suffix <- args[3]
outfile <- args[4]

# execute compile readcounts function
compile_readcounts(dir = indir,
                   strand = strand,
                   suffix = suffix,
                   outfile = outfile)
