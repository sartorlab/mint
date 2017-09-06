import os
from snakemake.workflow import expand
import glob
import getpass

################################################################################

# Get things from .yaml and do some massage
PROJECT_DIR = config.get("project_dir")
DATA_DIR = config.get("data_dir")
GENOME_DIR = config.get("genome_dir")
CHROM_LENGTHS = config.get("chrom_lengths")
BOWTIE2_INDEX = GENOME_DIR + "/genome"

GENOME = config.get("genome")

BIS_SAMPLE_DICT = config.get("bisulfite_samples")
PULL_SAMPLE_DICT = config.get("pulldown_samples")

BISULFITE_COMPARISONS = config.get("bisulfite_comparisons")
PULLDOWN_COMPARISONS = config.get("pulldown_comparisons")

################################################################################

# Create the path and switch
if not os.path.exists(PROJECT_DIR):
	os.makedirs(PROJECT_DIR)

os.chdir(PROJECT_DIR)

################################################################################

rule all:
    input:
        expand("bisulfite/01-raw_fastq/{sample}.fastq.gz", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/02-fastqc/{sample}_fastqc.zip", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/03-trim_galore/{sample}_trimmed.fq.gz", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/04-fastqc/{sample}_trimmed_fastqc.zip", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/05-bismark/{sample}_trimmed_bismark_bt2.bam", sample = BIS_SAMPLE_DICT.keys()),
        "bisulfite/06-multiqc/multiqc_report.html",
        expand("pulldown/01-raw_fastq/{sample}.fastq.gz", sample = PULL_SAMPLE_DICT.keys()),
        expand("pulldown/02-fastqc/{sample}_fastqc.zip", sample = PULL_SAMPLE_DICT.keys()),
        expand("pulldown/03-trim_galore/{sample}_trimmed.fq.gz", sample = PULL_SAMPLE_DICT.keys()),
        expand("pulldown/04-fastqc/{sample}_trimmed_fastqc.zip", sample = PULL_SAMPLE_DICT.keys()),
        expand("pulldown/05-bowtie2/{sample}_trimmed.fq.gz_aligned.bam", sample = PULL_SAMPLE_DICT.keys()),
        "pulldown/06-multiqc/multiqc_report.html"

################################################################################

rule bisulfite_setup:
    input:
        lambda wildcards: DATA_DIR + "/" + str(BIS_SAMPLE_DICT[wildcards.sample]) + ".fastq.gz"
    output:
        "bisulfite/01-raw_fastq/{sample}.fastq.gz"
    shell:  """
            ln -s {input} {output}
            """

rule pulldown_setup:
    input:
        lambda wildcards: DATA_DIR + "/" + str(PULL_SAMPLE_DICT[wildcards.sample]) + ".fastq.gz"
    output:
        "pulldown/01-raw_fastq/{sample}.fastq.gz"
    shell:  """
            ln -s {input} {output}
            """

################################################################################

rule bisulfite_raw_fastqc:
    input:
        "bisulfite/01-raw_fastq/{sample}.fastq.gz"
    output:
        "bisulfite/02-fastqc/{sample}_fastqc.zip"
    params:
        out_dir = "bisulfite/02-fastqc"
    shell:  """
            module purge && module load fastqc/0.11.4
            fastqc --format fastq --noextract --outdir {params.out_dir} {input}
            """

rule bisulfite_trimgalore:
    input:
        "bisulfite/01-raw_fastq/{sample}.fastq.gz"
    output:
        "bisulfite/03-trim_galore/{sample}_trimmed.fq.gz"
    params:
        adapter = "TGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGC",
        quality = 20,
        stringency = 6,
        error = 0.2,
        length = 25,
        out_dir = "bisulfite/03-trim_galore"
    shell:  """
            module purge && module load python/2.7.9 cutadapt/1.8.1 trim_galore/0.4.0
            trim_galore --quality {params.quality} --adapter {params.adapter} --stringency {params.stringency} -e {params.error} --gzip --length {params.length} --rrbs --output_dir {params.out_dir} {input}
            """

rule bisulfite_trim_fastqc:
    input:
        "bisulfite/03-trim_galore/{sample}_trimmed.fq.gz"
    output:
        "bisulfite/04-fastqc/{sample}_trimmed_fastqc.zip"
    shell:  """
            module purge && module load fastqc/0.11.4
            fastqc --format fastq --noextract --outdir bisulfite/04-fastqc {input}
            """

rule bisulfite_bismark:
    input:
        "bisulfite/03-trim_galore/{sample}_trimmed.fq.gz"
    output:
        bam = "bisulfite/05-bismark/{sample}_trimmed_bismark_bt2.bam",
        report = "bisulfite/05-bismark/{sample}_trimmed_bismark_bt2_SE_report.txt"
    params:
        genome_dir = GENOME_DIR,
        out_dir = "bisulfite/05-bismark"
    shell:  """
            module purge && module load bowtie2/2.2.4 samtools/1.2 bismark/0.16.3
            bismark --bowtie2 {params.genome_dir} --output_dir {params.out_dir} --temp_dir {params.out_dir} {input}
            samtools sort -f {output.bam} {output.bam}
            samtools index {output.bam}
            """

rule bisulfite_multiqc:
    input:
        expand("bisulfite/02-fastqc/{sample}_fastqc.zip", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/04-fastqc/{sample}_trimmed_fastqc.zip", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/05-bismark/{sample}_trimmed_bismark_bt2_SE_report.txt", sample = BIS_SAMPLE_DICT.keys())
    output:
        "bisulfite/06-multiqc/multiqc_report.html"
    params:
        out_dir = "bisulfite/06-multiqc"
    shell:  """
            multiqc --force ./bisulfite --outdir {params.out_dir}
            """

################################################################################

rule pulldown_raw_fastqc:
    input:
        "pulldown/01-raw_fastq/{sample}.fastq.gz"
    output:
        "pulldown/02-fastqc/{sample}_fastqc.zip"
    params:
        out_dir = "pulldown/02-fastqc"
    shell:  """
            module purge && module load fastqc/0.11.4
            fastqc --format fastq --noextract --outdir {params.out_dir} {input}
            """

rule pulldown_trimgalore:
    input:
        "pulldown/01-raw_fastq/{sample}.fastq.gz"
    output:
        "pulldown/03-trim_galore/{sample}_trimmed.fq.gz"
    params:
        adapter = "TGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGC",
        quality = 20,
        stringency = 6,
        error = 0.2,
        length = 25,
        out_dir = "pulldown/03-trim_galore"
    shell:  """
            module purge && module load python/2.7.9 cutadapt/1.8.1 trim_galore/0.4.0
            trim_galore --quality {params.quality} --illumina --stringency {params.stringency} -e {params.error} --gzip --length {params.length} --output_dir {params.out_dir} {input}
            """

rule pulldown_trim_fastqc:
    input:
        "pulldown/03-trim_galore/{sample}_trimmed.fq.gz"
    output:
        "pulldown/04-fastqc/{sample}_trimmed_fastqc.zip"
    shell:  """
            module purge && module load fastqc/0.11.4
            fastqc --format fastq --noextract --outdir pulldown/04-fastqc {input}
            """

rule pulldown_bowtie2:
    input:
        "pulldown/03-trim_galore/{sample}_trimmed.fq.gz"
    output:
        "pulldown/05-bowtie2/{sample}_trimmed.fq.gz_aligned.bam"
    params:
        align_summary = "pulldown/05-bowtie2/{sample}_bowtie2_summary.txt",
        bowtie2_index = BOWTIE2_INDEX
    shell:  """
            module purge && module load bowtie2/2.2.4 samtools/1.2
            bowtie2 -q --no-unal -x {params.bowtie2_index} -U {input} 2> {params.align_summary} | samtools view -bS - > {output}
            samtools sort -f {output} {output}
            samtools index {output}
            """

rule pulldown_multiqc:
    input:
        expand("pulldown/02-fastqc/{sample}_fastqc.zip", sample = PULL_SAMPLE_DICT.keys()),
        expand("pulldown/04-fastqc/{sample}_trimmed_fastqc.zip", sample = PULL_SAMPLE_DICT.keys()),
        expand("pulldown/05-bowtie2/{sample}_trimmed.fq.gz_aligned.bam", sample = PULL_SAMPLE_DICT.keys())
    output:
        "pulldown/06-multiqc/multiqc_report.html"
    params:
        out_dir = "pulldown/06-multiqc"
    shell:  """
            multiqc --force ./pulldown --outdir {params.out_dir}
            """
