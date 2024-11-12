#!/usr/bin/Rscript

source("scripts/ide_fxns.R")


#------------------------------------------------------#
# LOAD LIBRARIES

library(edgeR)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(genefilter)
library(stringr)
library(readxl)

options(scipen = 999)

#------------------------------------------------------#
# Get Command Line Arguments

args <- commandArgs(trailingOnly = TRUE)

print(args[1])
print(args[2])
print(args[3])
print(args[4])
print(args[5])
print(args[6])
print(args[7])
print(args[8])
print(args[9])
print(args[10])
print(args[11])

# raw counts file
raw_counts_file <- args[1]

# path to metadata.txt file
metadata_file <- args[2]

# path to gene info file
# var is set as "FALSE" if not provided
gene_info_file <- args[3]
if (gene_info_file == "FALSE") {
    print("No gene info file provided.")
    gene_info_file <- as.logical(gene_info_file)
    print(paste0("gene info file var is now: ", gene_info_file))
    gene_desc_col <- args[7]
} else {
    print(paste0("Using gene info file: ", gene_info_file))
    gene_id_col <- args[4]
    gene_name_col <- args[5]
    description <- as.logical(args[6])
    gene_desc_col <- args[7]
}
#gene_id_col <- args[4]
#gene_name_col <- args[5]
#description <- args[6]
#gene_desc_col <- args[7]

# output directory
outdir <- args[8]

# var of interest
myvar <- args[9]
# pca color var
pca_color_var <- args[10]
# pca label var
pca_label_var <- args[11]

# width of total counts barplot
tcounts_width <- 7

# type of plot files to save
plot_type <- "pdf"
outfile_prefix <- "./"

# column of metadata file to use as a group for low count filtering
# Various functions are kind of set up for this to be "Treatmnent"
#   - variance stabilizing transformation
#   - PCA plots
#lcf_col <- ""

#------------------------------------------------------#
# Read in data

print("reading in raw counts.")
# read in raw counts
raw_counts <- read.delim(raw_counts_file, header=T)
# set column 1 (genes) to rownames and delete column 1
rownames(raw_counts) <- raw_counts[, 1]
raw_counts[, 1] <- NULL

print("reading in metadata.")
# read in metadata
coldata <- read.delim(metadata_file, header=T)

print("reading in gene info file.")
# if provided, read in gene info and clean colnames
if (gene_info_file) {
  print("attempting to read gene_info_file")
  ref <- clean_gene_info(gene_info_file,
    gene_id_col = gene_id_col,
    gene_name_col = gene_name_col,
    description = description,
    gene_desc = gene_desc_col)
} else {
    ref <- NULL
}

# check that column names of raw_counts are in same order as 
# sample rows in metadata file

#------------------------------------------------------#
# Set working directory as outdir
print("creating outdir.")
dir.create(file.path(outdir), showWarnings = FALSE)
setwd(file.path(outdir))

# Plot Total Raw Counts Barplot
print("attempting raw counts plot.")
myplot <- plot_total_counts(raw_counts)
ggsave(filename = paste0("total_counts_barplot.", plot_type),
  plot = myplot,
  width = tcounts_width)

# Plot top 20 most abundant genes based on raw counts
print("attempting abundant genes plot")
myplot <- plot_top_count_genes(raw_counts, n = 20, gene_info = ref)
ggsave(filename=paste0("max_counts_barplot.", plot_type),
  plot = myplot)

#------------------------------------------------------#
# LOW-COUNT FILTERING

# Strategy: utilize edgeR's function "filterByExpr"

#counts.keep <- lcf_edger(raw_counts, group = coldata$Treatment)
counts.keep <- lcf_edger(raw_counts, group = coldata[[myvar]])
# export filtered counts
counts.keep.export <- counts.keep
counts.keep.export$GeneID <- rownames(counts.keep.export)
counts.keep.export <- counts.keep.export[ ,c(ncol(counts.keep.export),1:ncol(counts.keep.export)-1)]
write.table(counts.keep.export,
  "filtered_counts.txt",
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

# summarize low count filtering in ide log
sink("ide_complete.txt")
cat("LOW COUNT FILTERING\n")
cat("------------------------------\n")
print(paste0("Number of genes in raw counts table: ", nrow(raw_counts)))
print(paste0("Variable used to define groups for low count filtering: ", myvar))
print(paste0("Number of genes in filtered counts table: ", nrow(counts.keep)))
sink()

#------------------------------------------------------#
# PERFORM VARIANCE STABILIZING TRANSFORMATION

# Options of vst, rlog, or logcpm
# for now, do vst
sink("ide_complete.txt", append = TRUE)
cat("\nVARIANCE STABILIZING TRANSFORMATION\n")
cat("------------------------------\n")
vsd <- var_stable_transform(counts = counts.keep,
  coldata = coldata,
  method = "vst",
  formula_vars = c(0, myvar),
  blind = FALSE)
sink()

#------------------------------------------------------#
# SAMPLE CORRELATION HEATMAP

sample_dist_heatmap(norm_counts = t(assay(vsd)),
  outfile = "eucl_dist_heatmap",
  plot_type = plot_type)

#------------------------------------------------------#
# PCA

# use all features/genes
sink("ide_complete.txt", append = TRUE)
cat("\nPCA\n")
cat("------------------------------\n")
mypca <- perform_manual_pca(vsd = vsd, n=nrow(counts.keep))
sink()

# produce PCA biplots
myplot <- plot_pca(mypca = mypca,
  coldata = coldata,
  pc_a = "PC1",
  pc_b = "PC2",
  color_var = pca_color_var,
  label_var = pca_label_var)
ggsave(filename=paste0(outfile_prefix, "pca_biplot.pc1_pc2.", plot_type), plot=myplot)

myplot <- plot_pca(mypca = mypca,
  coldata = coldata,
  pc_a = "PC1",
  pc_b = "PC3",
  color_var = pca_color_var,
  label_var = pca_label_var)
ggsave(filename=paste0(outfile_prefix, "pca_biplot.pc1_pc3.", plot_type), plot=myplot)

myplot <- plot_pca(mypca = mypca,
  coldata = coldata,
  pc_a = "PC2",
  pc_b = "PC3",
  color_var = pca_color_var,
  label_var = pca_label_var)
ggsave(filename=paste0(outfile_prefix, "pca_biplot.pc2_pc3.", plot_type), plot=myplot)

sink("ide_complete.txt", append = TRUE)
cat("------------------------------\n")
cat("IDE is complete.\n")
sink()