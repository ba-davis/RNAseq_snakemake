
# Snakefile to analyze bulk RNAseq PE data

configfile:"proj_config.yaml"
#project_id = config["project_id"]


SAMPLES, = glob_wildcards("data/fastq/{sample}_R1.fastq.gz")
#COMPARISONS = config["contrasts"]

localrules: compile_readcounts

rule all:
    input:
        expand("data/fastqc/raw/{sample}_{dir}_fastqc.zip", sample = SAMPLES, dir = ["R1", "R2"]),
        expand("data/trimming/{sample}.paired_{dir}.fq.gz", sample = SAMPLES, dir = ["R1", "R2"]),
        expand("data/star/{sample}.bam", sample = SAMPLES),
        expand("data/star/{sample}.ReadsPerGene.out.tab", sample = SAMPLES),
        "data/counts/raw_counts.txt"

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

rule star:
    input:
        fwd = "data/trimming/{sample}.paired_R1.fq.gz",
        rev = "data/trimming/{sample}.paired_R2.fq.gz"
    output:
        bam_file = "data/star/{sample}.bam",
        counts = "data/star/{sample}.ReadsPerGene.out.tab"
    conda:
        "envs/star.yaml"
    params:
        genome = config["star_genome"],
        gtf = config["gtf"],
        out_prefix = "data/star/{sample}"
    shell:
        """
        STAR --runThreadN 8 --genomeDir {params.genome} --readFilesIn {input.fwd} {input.rev} --readFilesCommand zcat --sjdbGTFfile {params.gtf} --outFileNamePrefix {params.out_prefix} --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outFilterIntronMotifs RemoveNoncanonicalUnannotated --quantMode GeneCounts --twopassMode Basic --outReadsUnmapped Fastx
        """

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
