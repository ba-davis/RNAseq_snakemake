
# Snakefile to analyze bulk RNAseq PE data

configfile:"proj_config.yaml"
#project_id = config["project_id"]


SAMPLES, = glob_wildcards("data/fastq/{sample}_R1.fastq.gz")
#COMPARISONS = config["contrasts"]

localrules: compile_readcounts, collect_fqc_metrics, collect_trimmomatic_metrics, collect_star_metrics, join_metrics, ide

rule all:
    input:
        expand("data/fastqc/raw/{sample}_{dir}_fastqc.zip", sample = SAMPLES, dir = ["R1", "R2"]),
        expand("data/trimming/{sample}.paired_{dir}.fq.gz", sample = SAMPLES, dir = ["R1", "R2"]),
        #expand("data/star/{sample}.Aligned.sortedByCoord.out.bam", sample = SAMPLES),
        #expand("data/star/{sample}.ReadsPerGene.out.tab", sample = SAMPLES),
        "data/counts/raw_counts.txt",
        "data/fastqc/raw/fqc_stats_table.txt",
        "data/trimming/trimmomatic_stats_table.txt",
        #"data/star/star_stats_table.txt",
        #"data/preprocessing_metrics/metrics.txt",
        "data/ide/ide_complete.txt"

rule fastqc_raw:
    input:
        fwd = "data/fastq/{sample}_R1.fastq.gz",
        rev = "data/fastq/{sample}_R2.fastq.gz"
    output:
        fwd = "data/fastqc/raw/{sample}_R1_fastqc.zip",
        rev = "data/fastqc/raw/{sample}_R2_fastqc.zip"
    conda:
        "envs/fastqc.yaml"
    params:
        outdir = "data/fastqc/raw"
    shell:
        "fastqc -o {params.outdir} {input.fwd} {input.rev}"

# perform multiqc on fastqc results

rule trimmomatic:
    input:
        fwd = "data/fastq/{sample}_R1.fastq.gz",
	    rev = "data/fastq/{sample}_R2.fastq.gz"
    output:
        fwd = "data/trimming/{sample}.paired_R1.fq.gz",
	    rev = "data/trimming/{sample}.paired_R2.fq.gz",
	    fwd_unpaired = "data/trimming/{sample}.unpaired_R1.fq.gz",
	    rev_unpaired = "data/trimming/{sample}.unpaired_R2.fq.gz"
    conda:
        "envs/trimmomatic.yaml"
    params:
        trimmer = ["LEADING:3", "TRAILING:3", "SLIDINGWINDOW:4:15", "MINLEN:36"],
	    adapters = config["adapters"]
    shell:
        "trimmomatic PE -phred33 {input.fwd} {input.rev} {output.fwd} {output.fwd_unpaired} {output.rev} {output.rev_unpaired} {params.adapters} {params.trimmer}"

#rule star:
#    input:
#        fwd = "data/trimming/{sample}.paired_R1.fq.gz",
#        rev = "data/trimming/{sample}.paired_R2.fq.gz"
#    output:
#        bam_file = "data/star/{sample}.Aligned.sortedByCoord.out.bam",
#        counts = "data/star/{sample}.ReadsPerGene.out.tab",
#	logs = "data/star/{sample}.Log.final.out"
#    conda:
#        "envs/star.yaml"
#    params:
#        genome = config["star_genome"],
#        gtf = config["gtf"],
#        out_prefix = "data/star/{sample}."
#    shell:
#        """
#        STAR --runThreadN 8 --genomeDir {params.genome} --readFilesIn {input.fwd} {input.rev} --readFilesCommand zcat --sjdbGTFfile {params.gtf} --outFileNamePrefix {params.out_prefix} --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outFilterIntronMotifs RemoveNoncanonicalUnannotated --quantMode GeneCounts --twopassMode Basic --outReadsUnmapped Fastx
#        """

rule compile_readcounts:
    input:
        expand("data/star/{sample}.ReadsPerGene.out.tab", sample = SAMPLES)
    output:
        counts_table = "data/counts/raw_counts.txt"
    conda:
        "envs/bioinfo_r.yaml"
    params:
        indir = "data/star",
        strand = config["strand"],
        suffix = ".ReadsPerGene.out.tab"
    shell:
        """
        Rscript scripts/compile_readcounts.R {params.indir} {params.strand} {params.suffix} {output.counts_table}
        """

rule collect_fqc_metrics:
    input:
        expand("data/fastqc/raw/{sample}_{dir}_fastqc.zip", sample = SAMPLES, dir = ["R1", "R2"])
    output:
        "data/fastqc/raw/fqc_stats_table.txt"
    params:
        inpath = "data/fastqc/raw"
    shell:
        "scripts/collect_fastqc_metrics_PE.sh {params.inpath}"

rule collect_trimmomatic_metrics:
    input:
        expand("data/trimming/{sample}.paired_{dir}.fq.gz", sample = SAMPLES, dir = ["R1", "R2"])
    output:
        "data/trimming/trimmomatic_stats_table.txt"
    params:
        inpath = "logs",
        raw_file_path = "data/fastq",
        suffix = "_R1.fastq.gz",
        outfile = "data/trimming/trimmomatic_stats_table.txt"
    shell:
        "scripts/collect_trimmomatic_stats_PE.sh {params.inpath} {params.raw_file_path} {params.suffix} {params.outfile}"

rule collect_star_metrics:
    input:
        expand("data/star/{sample}.Log.final.out", sample = SAMPLES)
    output:
        "data/star/star_stats_table.txt"
    params:
        inpath = "data/star",
        outfile = "data/star/star_stats_table.txt"
    shell:
        "scripts/collect_star_stats.sh {params.inpath} {params.outfile}"

rule join_metrics:
    input:
        fqc = "data/fastqc/raw/fqc_stats_table.txt",
	    trim = "data/trimming/trimmomatic_stats_table.txt",
	    aln = "data/star/star_stats_table.txt"
    output:
        "data/preprocessing_metrics/metrics.txt"
    params:
        outfile = "data/preprocessing_metrics/metrics.txt"
    shell:
        "scripts/join_metrics.sh {input.fqc} {input.trim} {input.aln} {params.outfile}"

rule ide:
    input:
        counts = "data/counts/raw_counts.txt"
    output:
        "data/ide/ide_complete.txt"
    conda:
        "envs/RNAseq_v2020.yaml"
    params:
        metadata = config["metadata_file"],
        outdir = "data/ide",
        gene_info = config["gene_info_file"],
        gene_id_col = config["gene_id_col"],
        gene_name_col = config["gene_name_col"],
        description = config["description"],
        gene_desc_col = config["gene_desc_col"],
        var_of_int = config["var_of_int"],
        pca_color_var = config["pca_color_var"],
        pca_label_var = config["pca_label_var"]
    shell:
        "Rscript scripts/run_ide.R {input.counts} {params.metadata} {params.gene_info} {params.gene_id_col} {params.gene_name_col} {params.description} {params.gene_desc_col} {params.outdir} {params.var_of_int} {params.pca_color_var} {params.pca_label_var}"
