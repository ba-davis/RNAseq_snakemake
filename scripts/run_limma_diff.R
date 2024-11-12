

# use limma to fit a mixed effects model, treating donor as random effects
# instead of fixed effects (as they would be if included in the design matrix)

# conda activate RNAseq_v2020

#------------------------------------------------------#

# LOAD LIBRARIES
library(edgeR)
library(readxl)
library(statmod)
library(writexl)

#------------------------------------------------------#

# USER VARIABLES

# define paths to filtered counts file and metadata excel file
filtered_counts <- "../data/counts/filtered_counts.txt"
my_metadata <- "../metadata/rnaseq_metadata.xlsx"
gene_info <- "/home/groups/hoolock2/u0/genomes/ensembl/homo_sapiens/primary_assembly/annotation/GRCh38.103.gene_info.txt"
outfile_prefix <- "../data/diff/"

# read in data
counts.keep <- read.delim(filtered_counts, header=T)
rownames(counts.keep) <- counts.keep[ ,1]
counts.keep[ ,1] <- NULL

# read in coldata
coldata <- as.data.frame(read_xlsx(my_metadata))

# Clean coldata for better use
# change seq_run to categorical
coldata$seq_run <- as.character(coldata$seq_run)

ref <- read.delim(gene_info, header=TRUE)
# change colnames
colnames(ref) <- c("GeneID","gene.name","gene.type","gene.description")
# remove any duplicate GeneID rows
ref <- ref[!duplicated(ref$GeneID), ]

# Create design/model formula
my_design <- model.matrix(~ 0 + Group + seq_run, data=coldata)
colnames(my_design) <- gsub("Group", "", colnames(my_design))

# create outdir if doesn't exist
dir.create(outfile_prefix)

#------------------------------------------------------#

# create DGEList object
y <- DGEList(counts=counts.keep, group=coldata$Group)
# Normalize (TMM normalization for RNA composition)
y <- calcNormFactors(y, method = "TMM")
# Get logCPM counts
logCPM <- as.data.frame(cpm(y, log=TRUE))
logCPM$GeneID <- rownames(logCPM)
logCPM2 <- logCPM[ ,c(ncol(logCPM),1:ncol(logCPM)-1)]
write.table(logCPM2, "filtered_logCPM_counts.txt", sep="\t", col.names=T, row.names=F, quote=F)

#------------------------------------------------------#

# Define contrast matrix
contr.matrix <- makeContrasts(
  c1 = M6_rest_oud_hcv - M0_rest_oud_hcv,
  c2 = M3_rest_oud_hcv - M0_rest_oud_hcv,
  c3 = M6_LPS_oud_hcv - M0_LPS_oud_hcv,
  c4 = M3_LPS_oud_hcv - M0_LPS_oud_hcv,
  c5 = M6_rest_oud - M0_rest_oud,
  c6 = M3_rest_oud - M0_rest_oud,
  c7 = M6_LPS_oud - M0_LPS_oud,
  c8 = M3_LPS_oud - M0_LPS_oud,
  c9 = M0_rest_oud_hcv - M0_rest,
  c10 = M0_LPS_oud_hcv - M0_LPS,
  c11 = M0_rest_oud - M0_rest,
  c12 = M0_LPS_oud - M0_LPS,
  c13 = M0_rest_oud_hcv - M0_rest_oud,
  c14 = M0_LPS_oud_hcv - M0_LPS_oud,
  levels=my_design)

#------------------------------------------------------#

# Removing heteroscedascity with voom
pdf("voomplot.pdf")
v <- voom(y, my_design, plot=TRUE)
dev.off()

#------------------------------------------------------#

# Account for paired data (donors)
Donor <- factor(coldata$Donor)  # to be used as block-variable
table(Donor)

corfit <- duplicateCorrelation(v, my_design, block=Donor)
corfit$consensus  ## = 0.2816011, excpected to be positive in study with repeated sampling of the same patient at different time points
vfit <- lmFit(v, my_design, block=Donor, correlation=corfit$consensus) 
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)

summary(decideTests(efit))
head(efit$coefficients, n = 2)

#List all genes and corrsponding Log2FC, P-values, BH adjusted P-values to investigate further
res1 <- topTable(efit, coef=1, n=Inf, sort.by = "P")
res2 <- topTable(efit, coef=2, n=Inf, sort.by = "P")
res3 <- topTable(efit, coef=3, n=Inf, sort.by = "P")
res4 <- topTable(efit, coef=4, n=Inf, sort.by = "P")
res5 <- topTable(efit, coef=5, n=Inf, sort.by = "P")
res6 <- topTable(efit, coef=6, n=Inf, sort.by = "P")
res7 <- topTable(efit, coef=7, n=Inf, sort.by = "P")
res8 <- topTable(efit, coef=8, n=Inf, sort.by = "P")
res9 <- topTable(efit, coef=9, n=Inf, sort.by = "P")
res10 <- topTable(efit, coef=10, n=Inf, sort.by = "P")
res11 <- topTable(efit, coef=11, n=Inf, sort.by = "P")
res12 <- topTable(efit, coef=12, n=Inf, sort.by = "P")
res13 <- topTable(efit, coef=13, n=Inf, sort.by = "P")
res14 <- topTable(efit, coef=14, n=Inf, sort.by = "P")

# make a list of results dfs
dfs <- list(res1, res2, res3, res4, res5, res6, res7,
  res8, res9, res10, res11, res12, res13, res14)

# clean dfs
for (i in 1:length(dfs)) {
  dfs[[i]]$GeneID <- rownames(dfs[[i]])
  dfs[[i]] <- dfs[[i]][ ,c(ncol(dfs[[i]]), 1:ncol(dfs[[i]])-1)]
  rownames(dfs[[i]]) <- NULL
}

# add names to dfs
names(dfs) <- c("c1", "c2", "c3", "c4", "c5", "c6",
  "c7", "c8", "c9", "c10", "c11", "c12", "c13", "c14")

# merge with gene info
for (i in 1:length(dfs)) {
  dfs[[i]] <- merge(dfs[[i]], ref, by="GeneID")
}

# export excel files of diff results
write_xlsx(dfs, "ECP64.diff.results.xlsx")

# export txt files of results
for (i in 1:length(dfs)) {
  write.table(dfs[[i]], paste0(names(dfs)[i], ".diff.results.txt"), sep="\t", col.names=T, row.names=F, quote=F)
}
