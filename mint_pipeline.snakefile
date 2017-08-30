import os
from snakemake.workflow import expand
import glob
import getpass

# Get things from .yaml and do some massage
PROJECT_DIR = config.get("project_dir")
DATA_DIR = config.get("data_dir")
GENOME_DIR = config.get("genome_dir")
CHROM_LENGTHS = config.get("chrom_lengths")

GENOME = config.get("genome")
BIS_SAMPLE_DICT = config.get("bisulfite_samples")
PULL_SAMPLE_DICT = config.get("pulldown_samples")

BIS_SAMPLE_IDS = list(BIS_SAMPLE_DICT.keys())
PULL_SAMPLE_IDS = list(PULL_SAMPLE_DICT.keys())
ALL_SAMPLE_IDS = BIS_SAMPLE_IDS + PULL_SAMPLE_IDS

BIS_SAMPLE_ORIG = list(BIS_SAMPLE_DICT.values())
PULL_SAMPLE_ORIG = list(PULL_SAMPLE_DICT.values())
ALL_SAMPLE_ORIG = BIS_SAMPLE_ORIG + PULL_SAMPLE_ORIG

BISULFITE_COMPARISONS = config.get("bisulfite_comparisons")
PULLDOWN_COMPARISONS = config.get("pulldown_comparisons")

# print(PROJECT_DIR)
# print(DATA_DIR)
# print(GENOME_DIR)
# print(CHROM_LENGTHS)
# print(GENOME)
#
# print(BIS_SAMPLE_ORIG)
# print(PULL_SAMPLE_ORIG)
# print(ALL_SAMPLE_ORIG)
#
# print(BIS_SAMPLE_IDS)
# print(PULL_SAMPLE_IDS)
# print(ALL_SAMPLE_IDS)
#
# print(expand("{data_dir}/{sample}.fastq.gz", data_dir = DATA_DIR, sample = ALL_SAMPLE_IDS))
# print(expand("01-raw_fastq/{sample}.fastq.gz", sample = ALL_SAMPLE_IDS))
# print(expand("bisulfite/02-fastqc/{sample}_fastqc.zip", sample = BIS_SAMPLE_IDS))
# print(expand("bisulfite/03-trim_galore/{sample}_trimmed.fq.gz", sample = BIS_SAMPLE_IDS))
# print(expand("bisulfite/04-fastqc/{sample}_trimmed_fastqc.zip", sample = BIS_SAMPLE_IDS))
# print(expand("bisulfite/05-bismark/{sample}_trimmed.fq.gz_bismark_bt2.bam", sample = BIS_SAMPLE_IDS))
# print(expand("pulldown/02-fastqc/{sample}_fastqc.zip", sample = PULL_SAMPLE_IDS))
# print(expand("pulldown/03-trim_galore/{sample}_trimmed.fq.gz", sample = PULL_SAMPLE_IDS))
# print(expand("pulldown/04-fastqc/{sample}_trimmed_fastqc.zip", sample = PULL_SAMPLE_IDS))
# print(expand("pulldown/05-bowtie2/{sample}_trimmed.fq.gz_aligned.bam", sample = PULL_SAMPLE_IDS))

Create the path and switch
if not os.path.exists(PROJECT_DIR):
	os.makedirs(PROJECT_DIR)

os.chdir(PROJECT_DIR)

rule all:
    input:
        setup

rule setup:
    input:
        expand("{data_dir}/{sample}.fastq.gz", data_dir = DATA_DIR, sample = ALL_SAMPLE_IDS)
    output:
        expand("01-raw_fastq/{sample}.fastq.gz", sample = ALL_SAMPLE_IDS)
    shell:
        "ln -s {input} {output}"
#
# rule bisulfite_align:
#     output:
#         expand("bisulfite/02-fastqc/{sample}_fastqc.zip", sample = BIS_SAMPLE_IDS),
#         expand("bisulfite/03-trim_galore/{sample}_trimmed.fq.gz", sample = BIS_SAMPLE_IDS),
#         expand("bisulfite/04-fastqc/{sample}_trimmed_fastqc.zip", sample = BIS_SAMPLE_IDS),
#         expand("bisulfite/05-bismark/{sample}_trimmed.fq.gz_bismark_bt2.bam", sample = BIS_SAMPLE_IDS)
#
# rule pulldown_align:
#     output:
#         expand("pulldown/02-fastqc/{sample}_fastqc.zip", sample = PULL_SAMPLE_IDS),
#         expand("pulldown/03-trim_galore/{sample}_trimmed.fq.gz", sample = PULL_SAMPLE_IDS),
#         expand("pulldown/04-fastqc/{sample}_trimmed_fastqc.zip", sample = PULL_SAMPLE_IDS),
#         expand("pulldown/05-bowtie2/{sample}_trimmed.fq.gz_aligned.bam", sample = PULL_SAMPLE_IDS)

# print(GENOME)
# print(BISULFITE_SAMPLES)
# print(len(BISULFITE_SAMPLES))
# if PULLDOWN_SAMPLES is None:
#     print("No pulldown samples")
# else:
#     print(PULLDOWN_SAMPLES)
#
# for sample_id in BISULFITE_SAMPLES.keys():
#     print(sample_id)
#
# for sample_id in BISULFITE_SAMPLES.items():
#     print(sample_id[1])
#
# print(BISULFITE_COMPARISONS)
# print(PULLDOWN_COMPARISONS)
#
# print(PULLDOWN_COMPARISONS['IDH2mut_v_NBM']['exp_input'])
#
# for bis_comp in BISULFITE_COMPARISONS.keys():
#     print(bis_comp)
#     print(BISULFITE_COMPARISONS[bis_comp])
#     print(BISULFITE_COMPARISONS[bis_comp]['model'])
#     print(unpack(BISULFITE_COMPARISONS[bis_comp]))
#
# if JUST_BISULFITE:
#     print("Just Bisulfite!")
# else:
#     print("Hybrid!")
#
# print(expand("{one}/file.{two}", one = ['dir1','dir2'], two = ['ext1','ext2']))
