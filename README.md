# mint: A pipeline for analysis, integration, classification, and annotation of genome-wide DNA methylation and hydroxymethylation data

## Overview

DNA methylation is known to occur in a variety of forms. Canonical 5-methylcytosine (5mc) is the best studied, and has been shown to have a variety of roles in differentiation and regulation. More recently, studies of 5-hydroxymethylcytosine (5hmc) imply it is a stable base-modification with biological roles distinct from those of 5mc, including roles in DNA demethylation.

Current bisulfite-conversion + sequencing technologies (e.g. BS-seq and RRBS) are unable to distinguish between 5mc and 5hmc because both marks protect the cytosine from bisulfite-conversion. Newer technologies designed to detect only 5mc or 5hmc (e.g. oxBS-seq and TAB-seq) currently have low-reproducibility and conversion efficiency. There are, however, specific antibodies capable of pulling down either 5mc or 5hmc alone. These pulldown methods, unfortunately, do not achieve base-pair resolution, and are qualitative in nature.

We have developed the mint pipeline to classify regions of the genome into those containing 5mc, 5hmc, both, or neither by integrating information from multiple sequencing experiments.

## Tools and Outputs

For bisulfite-conversion data we use: [trim_galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) to trim adapter sequence from reads and trim based on quality scores, [bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/) to align reads to the reference genome and determine methylation rates at CpG resolution, and [methylSig](https://github.com/sartorlab/methylSig) to determine regions of differential methylation.

For pulldown data we use: [trim_galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) to trim adapter sequence from reads and trim based on quality scores, [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to align the reads to the reference genome, [MACS2](https://github.com/taoliu/MACS/tree/master/MACS2) to determine peaks for individual samples, and [PePr](https://github.com/shawnzhangyx/PePr) to determine regions of differential methylation.

Classification is done by intersecting the resulting regions from the above tools by using [bedtools](https://bedtools.readthedocs.org/en/latest/) and [bedops](https://bedops.readthedocs.org/en/latest/).

Sites and regions of methylation, differential methylation, and methylation classifications are annotated to genomic regions using the [annotatr](https://www.github.com/rcavalcante/annotatr) R package. A variety of summary graphs are automatically generated at the sample and comparison level.

A track hub is created for the [UCSC Genome Browser](https://genome.ucsc.edu) that includes the following tracks:
	* Percent methylation tracks from the bismark methylation extractor.
	* Pileup tracks of pulldown coverage for pulldowns and inputs.
	* Peak tracks of methylated regions from pulldowns.
	* Percent methylation difference between groups from methylSig.
	* Peak tracks of differentially methylated regions from PePr.
	* Classification tracks.

# Usage

The mint pipeline can be used for any combination of the following experimental setups and analysis workflows. At present, only single-end reads are supported.

## Supported Experimental Setups

* **Hybrid** experimental setups include data from both bisulfite-conversion experiments (BS-seq, WGBS, RRBS, etc.) and pulldown experiments (hMeDIP-seq, hMeSeal, etc.).

* **Pulldown** experimental setups include only data from pulldowns with an antibody such as MeDIP-seq, hMeDIP-seq, etc.

* *NOTE*: Purely bisulfite setups (e.g. RRBS + oxBS-seq, RRBS + TAB-seq, etc.) are not currently supported.

## Supported Analysis Workflows

* **Sample-wise** analysis workflows combine information from experiments at the sample level only.

* **Comparison-wise** analysis workflows combine information from experiments at the level of groups to be compared for differential methylation.

## Starting a project

0. Install prerequisites found in [VERSIONS](https://github.com/sartorlab/mint/blob/make-refactor/VERSIONS.md).
1. `git clone https://github.com/sartorlab/mint` in the desired install directory.
2. In `mint/` do `mkdir projects`.
3. Put a tab-delimited annotation file `project_name_annotation.txt` in the `mint/projects` directory. It **must** have 9 columns:
	1. `projectID`: The name of the project.
	2. `sampleID`: An alphanumeric ID (perhaps from SRA, GEO, a sequencing core, etc.). Typically these will be the names of the `.fastq` files.
	3. `humanID`: The human readable ID for the sample.
	4. `pulldown`: A binary value indicating whether the sample is the result of a pulldown experiment (1) or not (0).
	5. `bisulfite`: A binary value indicating whether the sample is the result of a bisulfite-conversion experiment (1) or not (0).
	6. `mc`: A binary value indicating whether the sample represents 5mc methylation.
	7. `hmc`: A binary value indicating whether the sample represents 5hmc methylation.
	8. `input`: A binary value indicating whether the sample represents an input.
	9. `group`: A binary value indicating which samples belong to one of two groups.
4. Know where the data matching the `sampleID` column in the annotation file is located.
5. In `mint/` do `Rscript init.R --project name --genome g --datapath path/to/data/matching/sampleID.fastq.gz`
	* This script does the following:
		* Creates the directory structure for the project based on the project annotation file.
		* Copies the project annotation file to `mint/projects/name/name_annotation.txt`.
		* Copies `template_config.mk` and `template_makefile` into `mint/projects/name`.
		* Dynamically generates adds the appropriate rules to the `makefile` based on the project annotation file.
		* Generates UCSC Genome Browser track hub structure and files based on the project annotation file.
		* Generates PBS job scripts in `mint/projects/name/pbs_jobs` for use on computing clusters.
6. In `mint/projects/name` modify the `config.mk` file to reflect the location of tools, and the desired parameters for analysis.
7. Run the analyses. `mint` automatically generates PBS scripts (which need to be customized) for computing clusters in `mint/projects/name/pbs_jobs`, or the following make commands from `mint/projects/name` can be used. NOTE: Depending on the experimental and workflow setup encoded in the project annotation file, some of the `make` commands might not be available for the project.
	* `make bisulfite_align`
	* `make pulldown_align`
	* `make pulldown_sample`
	* `make bisulfite_compare`
	* `make pulldown_compare`
	* `make sample_classification`
	* `make compare_classification`
	* Add the `-n` flag to `make` to see what commands will be run before committing resources.
	* Add the `-j` flag followed by a number to run commands in parallel depending on the computing architecture used.

### Project annotation file examples

Bisulfite-conversion experiments that represent both mc and hmc should have a 1 in each column. Input pulldowns can be matched to the pulldown (e.g. `mc=1` and `hmc=0` and `input=1`) or shared among pulldowns. In the case of sharing between antibodies, the entry should be doubled, as in the first example.

An example of a pulldown experimental setup with sample-wise analysis. NOTE: The inputs are shared between the mc and hmc pulldown. To represent that, the input block is duplicated, within one `mc=1` and `hmc=0` and within the other `mc=0` and `hmc=1`.

```{bash}
projectID       sampleID        humanID pulldown        bisulfite       mc      hmc     input   group
GSE63743        SRR1686689      preeclamptic_1  1       0       0       1       0       0
GSE63743        SRR1686690      preeclamptic_2  1       0       0       1       0       0
GSE63743        SRR1686693      normal_1        1       0       0       1       0       0
GSE63743        SRR1686694      normal_2        1       0       0       1       0       0
GSE63743        SRR1686697      preeclamptic_1  1       0       1       0       0       0
GSE63743        SRR1686698      preeclamptic_2  1       0       1       0       0       0
GSE63743        SRR1686701      normal_1        1       0       1       0       0       0
GSE63743        SRR1686702      normal_2        1       0       1       0       0       0
GSE63743        SRR1686705      preeclamptic_1  1       0       1       0       1       0
GSE63743        SRR1686706      preeclamptic_2  1       0       1       0       1       0
GSE63743        SRR1686709      normal_1        1       0       1       0       1       0
GSE63743        SRR1686710      normal_2        1       0       1       0       1       0
GSE63743        SRR1686705      preeclamptic_1  1       0       0       1       1       0
GSE63743        SRR1686706      preeclamptic_2  1       0       0       1       1       0
GSE63743        SRR1686709      normal_1        1       0       0       1       1       0
GSE63743        SRR1686710      normal_2        1       0       0       1       1       0
```

An example of a hybrid experimental setup with comparison-wise analysis. Note the added lines for comparisons indicate that pulldown and bisulfite experiments should be compared according to groups 0 and 1.

```{bash}
projectID	sampleID	humanID	pulldown	bisulfite	mc	hmc	input	group
GSE52945	SRR1041959	IDH2mut_1	1	0	0	1	0	1
GSE52945	SRR1041960	IDH2mut_2	1	0	0	1	0	1
GSE52945	SRR1041977	IDH2mut_1	1	0	0	1	1	1
GSE52945	SRR1041978	IDH2mut_2	1	0	0	1	1	1
GSE52945	SRR1041992	IDH2mut_1	0	1	1	1	0	1
GSE52945	SRR1041993	IDH2mut_2	0	1	1	1	0	1
GSE52945	SRR1638715	NBM_1	1	0	0	1	0	0
GSE52945	SRR1638716	NBM_2	1	0	0	1	0	0
GSE52945	SRR1638720	NBM_1	1	0	0	1	1	0
GSE52945	SRR1638721	NBM_2	1	0	0	1	1	0
GSE52945	SRR1638726	NBM_2	0	1	1	1	0	0
GSE52945	SRR1638727	NBM_1	0	1	1	1	0	0
GSE52945	comparison	IDH2mut_v_NBM	1	0	0	1	0	0,1
GSE52945	comparison	IDH2mut_v_NBM	0	1	1	1	0	0,1
```

A more complicated hybrid experimental setup where multiple groups are to be compared in the same analysis. In the first comparisons the groups 0 and 1 are used, whereas in the second comparisons use groups 2 and 3. An arbitrary number of comparisons are supported.

```{bash}
projectID	sampleID	humanID	pulldown	bisulfite	mc	hmc	input	group
hnscc_13	Sample_42741	HPV+2	0	1	1	1	0	1,2
hnscc_13	Sample_42742	HPV+5	0	1	1	1	0	1,2
hnscc_13	Sample_45194	HPV-1	0	1	1	1	0	0
hnscc_13	Sample_45195	HPV+1	0	1	1	1	0	1,2
hnscc_13	Sample_45196	HPV+3	0	1	1	1	0	1,3
hnscc_13	Sample_45197	HPV+4	0	1	1	1	0	1,3
hnscc_13	Sample_45199	HPV+8	0	1	1	1	0	1,3
hnscc_13	Sample_45200	HPV-2	0	1	1	1	0	0
hnscc_13	Sample_45201	HPV-4	0	1	1	1	0	0
hnscc_13	Sample_45202	HPV+7	0	1	1	1	0	1,3
hnscc_13	Sample_45203	HPV-7	0	1	1	1	0	0
hnscc_13	Sample_45204	HPV-8	0	1	1	1	0	0
hnscc_13	Sample_51324	HPV-8	1	0	0	1	1	0
hnscc_13	Sample_51327	HPV+8	1	0	0	1	0	1,3
hnscc_13	Sample_51325	HPV-8	1	0	0	1	0	0
hnscc_13	Sample_51326	HPV+8	1	0	0	1	1	1,3
hnscc_13	Sample_49651	HPV-1	1	0	0	1	1	0
hnscc_13	Sample_49654	HPV+1	1	0	0	1	0	1,2
hnscc_13	Sample_49652	HPV-1	1	0	0	1	0	0
hnscc_13	Sample_49653	HPV+1	1	0	0	1	1	1,2
hnscc_13	Sample_49655	HPV+2	1	0	0	1	1	1
hnscc_13	Sample_49658	HPV-2	1	0	0	1	0	0
hnscc_13	Sample_49656	HPV+2	1	0	0	1	0	1,2
hnscc_13	Sample_49657	HPV-2	1	0	0	1	1	0
hnscc_13	Sample_49659	HPV+3	1	0	0	1	1	1,3
hnscc_13	Sample_49662	HPV+4	1	0	0	1	0	1,3
hnscc_13	Sample_49660	HPV+3	1	0	0	1	0	1,3
hnscc_13	Sample_49661	HPV+4	1	0	0	1	1	1,3
hnscc_13	Sample_49663	HPV+5	1	0	0	1	1	1,2
hnscc_13	Sample_49664	HPV+5	1	0	0	1	0	1,2
hnscc_13	Sample_51316	HPV-4	1	0	0	1	1	0
hnscc_13	Sample_51317	HPV-4	1	0	0	1	0	0
hnscc_13	Sample_51320	HPV+7	1	0	0	1	1	1,3
hnscc_13	Sample_51323	HPV-7	1	0	0	1	0	0
hnscc_13	Sample_51321	HPV+7	1	0	0	1	0	1,3
hnscc_13	Sample_51322	HPV-7	1	0	0	1	1	0
hnscc_13	comparison1	HPV+_v_HPV-	1	0	0	1	0	0,1
hnscc_13	comparison1	HPV+_v_HPV-	0	1	1	1	0	0,1
hnscc_13	comparison2	int+_v_int-	1	0	0	1	0	2,3
hnscc_13	comparison2	int+_v_int-	0	1	1	1	0	2,3
```
