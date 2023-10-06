

############--------------------------------------------#
# FUNCTION # compile_readcounts: to create counts table #
############--------------------------------------------#
# INPUTS: dir: path to directory containing STAR ReadsPerGene.out.tab files
#         strand: number signifying which column of counts to keep
#                 compare N_noFeature line in ReadsPerGene files
#                 choose 4 if col4 has much lower number than col3
#                 if col3 and col4 are similar, choose 2
#                 often either 2 or 4
#         outfile: name of output counts table text file
#         suffix: ending of ReadsPerGene files to grab (and to remove)
# OUTPUTS: 1. a raw counts table tab-separated text file
#          2. returns a raw counts data.frame object
compile_readcounts <- function(dir,
                               strand,
                               suffix = ".ReadsPerGene.out.tab",
                               outfile = "raw_counts.txt") {
  # create a list that contains the path to all input files
  counts_list <- list.files(dir, paste0(suffix, "$"), full = TRUE)
  names(counts_list) <- counts_list
  names(counts_list) <- gsub(suffix, "", names(counts_list))
  names(counts_list) <- gsub(dir, "", names(counts_list))
  names(counts_list) <- gsub("/", "", names(counts_list))

  # read in each counts file as a data frame and store in a list of df's
  counts_dfs <- parallel::mclapply(counts_list,
                                   function(x) read.table(x), mc.cores = 8)

  # convert the correct counts column to numeric
  for (i in 1:length(counts_dfs)) {
    counts_dfs[[i]][, strand] <- as.numeric(counts_dfs[[i]][, strand])
  }

  # change the name of the correct column to be the sample name
  for (i in 1:length(counts_dfs)) {
    colnames(counts_dfs[[i]])[strand] <- names(counts_dfs)[i]
  }

  # change the name of column 1 to "GeneID"
  for (i in 1:length(counts_dfs)) {
    colnames(counts_dfs[[i]])[1] <- "GeneID"
  }

  # drop unnecessary counts columns
  for (i in 1:length(counts_dfs)) {
    counts_dfs[[i]] <- subset(counts_dfs[[i]],
                              select = colnames(counts_dfs[[i]][c(1, strand)]))
  }

  # remove rows 1,2,3,4
  for (i in 1:length(counts_dfs)) {
    counts_dfs[[i]] <- counts_dfs[[i]][-c(1, 2, 3, 4), ]
  }

  # merge all data frames on "GeneID"
  my_table <- Reduce(function(df1, df2) merge(df1,
                                              df2,
                                              by = "GeneID",
                                              all = TRUE), counts_dfs)

  # export the final table
  write.table(my_table,
    outfile,
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE,
  )

  #return(my_table)
}
