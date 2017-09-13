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

BIS_COMPARISONS_DICT = config.get("bisulfite_comparisons")
PULL_COMPARISONS_DICT = config.get("pulldown_comparisons")

EXECUTE_DIR = os.getcwd()

print(expand("bisulfite/10-diffMeth/{comparison}_dss_significant.txt", comparison = BIS_COMPARISONS_DICT.keys()))
print(str(BIS_COMPARISONS_DICT['IDH2mut_v_NBM']['exp'] + BIS_COMPARISONS_DICT['IDH2mut_v_NBM']['con']))
print(BIS_COMPARISONS_DICT['IDH2mut_v_NBM']['model'])
print(BIS_COMPARISONS_DICT['IDH2mut_v_NBM']['contrast'])
print(BIS_COMPARISONS_DICT['IDH2mut_v_NBM']['covariates'])
print(BIS_COMPARISONS_DICT['IDH2mut_v_NBM']['numerical_covariates'])
print(BIS_COMPARISONS_DICT['IDH2mut_v_NBM']['groups'])
print(BIS_COMPARISONS_DICT['IDH2mut_v_NBM']['interpretation'])

################################################################################

# Create the path and switch
if not os.path.exists(PROJECT_DIR):
	os.makedirs(PROJECT_DIR)
os.chdir(PROJECT_DIR)

################################################################################
################################################################################

rule all:
    input:
        expand("bisulfite/01-raw_fastq/{sample}.fastq.gz", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/02-fastqc/{sample}_fastqc.zip", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/03-trim_galore/{sample}_trimmed.fq.gz", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/04-fastqc/{sample}_trimmed_fastqc.zip", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/05-bismark/{sample}_trimmed_bismark_bt2.bam", sample = BIS_SAMPLE_DICT.keys()),
        "bisulfite/06-multiqc/multiqc_report.html",
        expand("bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.bedGraph.gz", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.bismark.cov.gz", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.CpG_report.txt.gz", sample = BIS_SAMPLE_DICT.keys()),
        expand("RData/{sample}_trimmed_bismark_annotatr_analysis.RData", sample = BIS_SAMPLE_DICT.keys()),
        expand("trackhub/{sample}_trimmed_bismark_bt2.bw", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/09-simple_classification/{sample}_bismark_simple_classification.bed", sample = BIS_SAMPLE_DICT.keys()),
        expand("RData/{sample}_bismark_simple_classification_annotatr_analysis.RData", sample = BIS_SAMPLE_DICT.keys()),
        expand("trackhub/{sample}_bismark_simple_classification.bb", sample = BIS_SAMPLE_DICT.keys()),
        "bisulfite/08-multiqc/multiqc_report.html",
        expand("pulldown/01-raw_fastq/{sample}.fastq.gz", sample = PULL_SAMPLE_DICT.keys()),
        expand("pulldown/02-fastqc/{sample}_fastqc.zip", sample = PULL_SAMPLE_DICT.keys()),
        expand("pulldown/03-trim_galore/{sample}_trimmed.fq.gz", sample = PULL_SAMPLE_DICT.keys()),
        expand("pulldown/04-fastqc/{sample}_trimmed_fastqc.zip", sample = PULL_SAMPLE_DICT.keys()),
        expand("pulldown/05-bowtie2/{sample}_trimmed.fq.gz_aligned.bam", sample = PULL_SAMPLE_DICT.keys()),
        "pulldown/06-multiqc/multiqc_report.html"

################################################################################
################################################################################

rule bisulfite_align:
    input:
        expand("bisulfite/01-raw_fastq/{sample}.fastq.gz", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/02-fastqc/{sample}_fastqc.zip", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/03-trim_galore/{sample}_trimmed.fq.gz", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/04-fastqc/{sample}_trimmed_fastqc.zip", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/05-bismark/{sample}_trimmed_bismark_bt2.bam", sample = BIS_SAMPLE_DICT.keys()),
        "bisulfite/06-multiqc/multiqc_report.html"

rule bisulfite_sample:
    input:
        expand("bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.bedGraph.gz", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.bismark.cov.gz", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.CpG_report.txt.gz", sample = BIS_SAMPLE_DICT.keys()),
        expand("RData/{sample}_trimmed_bismark_annotatr_analysis.RData", sample = BIS_SAMPLE_DICT.keys()),
        expand("trackhub/{sample}_trimmed_bismark_bt2.bw", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/09-simple_classification/{sample}_bismark_simple_classification.bed", sample = BIS_SAMPLE_DICT.keys()),
        expand("RData/{sample}_bismark_simple_classification_annotatr_analysis.RData", sample = BIS_SAMPLE_DICT.keys()),
        expand("trackhub/{sample}_bismark_simple_classification.bb", sample = BIS_SAMPLE_DICT.keys()),
        "bisulfite/08-multiqc/multiqc_report.html"

rule bisulfite_compare:
    input:
        expand("bisulfite/10-diffMeth/{comparison}_dss_significant.txt", comparison = BIS_COMPARISONS_DICT.keys())

################################################################################
################################################################################

rule pulldown_align:
    input:
        expand("pulldown/01-raw_fastq/{sample}.fastq.gz", sample = PULL_SAMPLE_DICT.keys()),
        expand("pulldown/02-fastqc/{sample}_fastqc.zip", sample = PULL_SAMPLE_DICT.keys()),
        expand("pulldown/03-trim_galore/{sample}_trimmed.fq.gz", sample = PULL_SAMPLE_DICT.keys()),
        expand("pulldown/04-fastqc/{sample}_trimmed_fastqc.zip", sample = PULL_SAMPLE_DICT.keys()),
        expand("pulldown/05-bowtie2/{sample}_trimmed.fq.gz_aligned.bam", sample = PULL_SAMPLE_DICT.keys()),
        "pulldown/06-multiqc/multiqc_report.html"

################################################################################
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

rule bisulfite_align_raw_fastqc:
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

rule bisulfite_align_trimgalore:
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

rule bisulfite_align_trim_fastqc:
    input:
        "bisulfite/03-trim_galore/{sample}_trimmed.fq.gz"
    output:
        "bisulfite/04-fastqc/{sample}_trimmed_fastqc.zip"
    shell:  """
            module purge && module load fastqc/0.11.4
            fastqc --format fastq --noextract --outdir bisulfite/04-fastqc {input}
            """

rule bisulfite_align_bismark:
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

rule bisulfite_align_multiqc:
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

rule bisulfite_sample_extractor:
    input:
        "bisulfite/05-bismark/{sample}_trimmed_bismark_bt2.bam"
    output:
        bdg = "bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.bedGraph.gz",
        cov = "bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.bismark.cov.gz",
        cpg = "bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.CpG_report.txt.gz",
        mbtxt = "bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.M-bias.txt",
        mbplot = "bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.M-bias_R1.png",
        srep = "bisulfite/07-methCall/{sample}_trimmed_bismark_bt2_splitting_report.txt"
    params:
        genome_dir = GENOME_DIR,
        min_cov = 5,
        out_dir = "bisulfite/07-methCall"
    threads: 4
    shell:  """
            module purge && module load bowtie2/2.2.4 samtools/1.2 bismark/0.16.3
            bismark_methylation_extractor --single-end --gzip --bedGraph --cutoff {params.min_cov} --cytosine_report --genome_folder {params.genome_dir} --multicore {threads} --output {params.out_dir} {input}
            """

rule bisulfite_sample_to_annotatr:
    input:
        "bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.bismark.cov.gz"
    output:
        temp("bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.bismark.cov")
    shell: """
            gunzip -c {input} | awk -v OFS="\t" '{{print $1, $2 - 1, $3, ".", $4, ".", $5 + $6}}' > {output}
            """

rule bisulfite_sample_annotatr:
    input:
        "bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.bismark.cov"
    output:
        "RData/{sample}_trimmed_bismark_annotatr_analysis.RData"
    params:
        genome = GENOME,
        exec_dir = EXECUTE_DIR
    shell:  """
            module purge && module load java/1.8.0 gcc/4.9.3 R/3.4.0
            Rscript {params.exec_dir}/scripts/annotatr_annotations.R --file {input} --genome {params.genome} --annot_type bismark --group1 NULL --group0 NULL
            """

rule bisulfite_sample_pretrack:
    input:
        "bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.bedGraph.gz"
    output:
        temp("bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.bedGraph")
    shell:  """
            gunzip -c {input} | awk 'NR > 1 {{print $0}}' | sort -T . -k1,1 -k2,2n > {output}
            """

rule bisulfite_sample_track:
    input:
        "bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.bedGraph"
    output:
        "trackhub/{sample}_trimmed_bismark_bt2.bw"
    params:
        chrom_lengths = CHROM_LENGTHS
    shell:  """
            bedGraphToBigWig {input} {params.chrom_lengths} {output}
            """

rule bisulfite_sample_multiqc:
    input:
        expand("bisulfite/02-fastqc/{sample}_fastqc.zip", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/04-fastqc/{sample}_trimmed_fastqc.zip", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/05-bismark/{sample}_trimmed_bismark_bt2_SE_report.txt", sample = BIS_SAMPLE_DICT.keys()),
        expand("bisulfite/07-methCall/{sample}_trimmed_bismark_bt2_splitting_report.txt", sample = BIS_SAMPLE_DICT.keys())
    output:
        "bisulfite/08-multiqc/multiqc_report.html"
    params:
        out_dir = "bisulfite/08-multiqc"
    shell:  """
            multiqc --force ./bisulfite --outdir {params.out_dir}
            """

################################################################################

rule bisulfite_sample_simple:
    input:
        "bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.bedGraph.gz"
    output:
        "bisulfite/09-simple_classification/{sample}_bismark_simple_classification.bed"
    params:
        exec_dir = EXECUTE_DIR
    shell:  """
            module purge && module load java/1.8.0 gcc/4.9.3 R/3.4.0
            Rscript {params.exec_dir}/scripts/classify_simple.R --inFile {input} --outFile {output} --group1 NULL --group0 NULL
            """

rule bisulfite_sample_simple_annotatr:
    input:
        "bisulfite/09-simple_classification/{sample}_bismark_simple_classification.bed"
    output:
        "RData/{sample}_bismark_simple_classification_annotatr_analysis.RData"
    params:
        genome = GENOME,
        exec_dir = EXECUTE_DIR
    shell:  """
            module purge && module load java/1.8.0 gcc/4.9.3 R/3.4.0
            Rscript {params.exec_dir}/scripts/annotatr_annotations.R --file {input} --genome {params.genome} --annot_type simple_bisulfite_bismark --group1 NULL --group0 NULL
            """

rule bisulfite_sample_simple_track:
    input:
        "bisulfite/09-simple_classification/{sample}_bismark_simple_classification.bed"
    output:
        "trackhub/{sample}_bismark_simple_classification.bb"
    params:
        chrom_lengths = CHROM_LENGTHS
    shell:  """
            bedToBigBed {input} {params.chrom_lengths} {output}
            """

################################################################################

rule bisulfite_compare_predss:
    input:
        "bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.CpG_report.txt.gz"
    output:
        "bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.CpG_report_for_dss.txt"
    params:
        min_cov = 5
    shell:  """
            awk -v OFS='\t' -v MIN_COV={params.min_cov} '$4 + $5 >= MIN_COV {print $1, $2, $3, $4 + $5, $4}' <(gunzip -c {input}) | sort -T . -k1,1 -k2,2n > {output}
            """

rule bisulfite_compare_dss:
    input:
        lambda wildcards: expand("bisulfite/07-methCall/{sample}_trimmed_bismark_bt2.CpG_report_for_dss.txt", sample = str(BIS_COMPARISONS_DICT[{wildcards.comparison}]['exp'] +  BIS_COMPARISONS_DICT[{wildcards.comparison}]['con']))
    output:
        "bisulfite/10-diffMeth/{comparison}_dss_significant.txt"
    params:
        genome = GENOME,
        exec_dir = EXECUTE_DIR,
        names = str(lambda wildcards: str(BIS_COMPARISONS_DICT[{wildcards.comparison}]['exp'] + BIS_COMPARISONS_DICT[{wildcards.comparison}]['con'])),
        model = str(lambda wildcards: BIS_COMPARISONS_DICT[{wildcards.comparison}]['model']),
        contrast = str(lambda wildcards: BIS_COMPARISONS_DICT[{wildcards.comparison}]['contrast']),
        covariates = str(lambda wildcards: BIS_COMPARISONS_DICT[{wildcards.comparison}]['covariates']),
        numerical_covariates = str(lambda wildcards: BIS_COMPARISONS_DICT[{wildcards.comparison}]['numerical_covariates']),
        groups = str(lambda wildcards: BIS_COMPARISONS_DICT[{wildcards.comparison}]['groups']),
        interpretation = str(lambda wildcards: BIS_COMPARISONS_DICT[{wildcards.comparison}]['interpretation']),
        dm_diff = 10,
        dm_fdr = 0.05,
        dm_p = 0.005,
        destrand = 'TRUE',
        tilewidth = 50,
        outprefix = str(lambda wildcards: {wildcards.comparison})
    shell:  """
            module purge && module load java/1.8.0 gcc/4.9.3 R/3.4.0
            Rscript {params.exec_dir}/scripts/dss_run.R --genome {params.genome} --files {input} --samplenames {params.names} --model {params.model} --groups {params.groups} --contrast {params.contrast} --covariates {params.covariates} --covIsNumeric {params.numerical_covariates} --interpretation {params.interpretation} --outprefix {params.outprefix} --destrand {params.destrand} --tilewidth {params.tilewidth} --methdiffthreshold {params.dm_diff} --FDRthreshold {params.dm_fdr}  --pvalthreshold {params.dm_p} --quiet FALSE
            """

################################################################################
################################################################################
################################################################################

rule pulldown_align_raw_fastqc:
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

rule pulldown_align_trimgalore:
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

rule pulldown_align_trim_fastqc:
    input:
        "pulldown/03-trim_galore/{sample}_trimmed.fq.gz"
    output:
        "pulldown/04-fastqc/{sample}_trimmed_fastqc.zip"
    shell:  """
            module purge && module load fastqc/0.11.4
            fastqc --format fastq --noextract --outdir pulldown/04-fastqc {input}
            """

rule pulldown_align_bowtie2:
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

rule pulldown_align_multiqc:
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
