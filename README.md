# mint: a pipeline for methylation integration

# Overview

DNA methylation is known to occur in a variety of forms. Canonical 5-methylcytosine (5mC) is the best studied, and has been shown to have a variety of roles in differentiation and regulation. More recently, studies of 5-hydroxymethylcytosine (5hmC) imply it is a stable base-modification with biological roles distinct from those of 5mC.

Current bisulfite-conversion + sequencing technologies (e.g. BS-seq and RRBS) are unable to distinguish between 5mC and 5hmC because both marks protect the cytosine from bisulfite-conversion. Newer technologies designed to detect only 5hmC (e.g. oxBS-seq and TAB-seq) currently have low-reproducibility and conversion efficiency. There are, however, specific antibodies capable of pulling down either 5mC or 5hmC alone. These pulldown methods, unfortunately, do not achieve base-pair resolution, and are qualitative in nature.

We have developed the mint pipeline to classify regions of the genome into those containing 5mC, 5hmC, both, or neither on the basis of combined information from multiple sequencing experiments.

# Methods

For bisulfite-conversion data we use [cutadapt](https://cutadapt.readthedocs.org/en/stable/) to trim adapter sequence from reads and [bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/) to align reads to the reference genome and determine methylation rates at CpG resolution. We use [methylSig](https://github.com/sartorlab/methylSig) to determine regions of differential methylation.

For pulldown data we use [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to align the reads to the reference genome, [MACS2](https://github.com/taoliu/MACS/tree/master/MACS2) to determine peaks for individual samples, and [PePr](https://github.com/shawnzhangyx/PePr) to determine regions of differential methylation.

Classification is done by intersecting the resulting regions from the above tools using the [GenomicRanges](http://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) package in R, or by using [bedtools](https://bedtools.readthedocs.org/en/latest/).

# Usage

The mint pipeline can be used for any combination of the following two experimental setups and two analysis workflows. At present, we support only single-end sequencing data.

## Supported Experimental Setups

* **Hybrid** experimental setups include data from both bisulfite-conversion experiments (BS-seq, WGBS, RRBS, etc.) and pulldown experiments (hMeDIP-seq, hMeSeal, etc.).

* **Pulldown** experimental setups include only data from pulldowns with an antibody such as MeDIP-seq, hMeDIP-seq, etc.

* Note, purely bisulfite-conversion workflows (e.g. RRBS + oxBS-seq, RRBS + TAB-seq, etc.) are not currently supported.

## Supported Analysis Workflows

* **Sample-wise** analysis workflows combine information from experiments at the sample level only.

* **Comparison-wise** analysis workflows combine information from experiments at the level of groups to be compared for differential methylation.

## Initializing a project

  1. After obtaining mint, and installing all dependencies, users can navigate to the `mint/scripts/` directory and do the following:
  ```{bash}
  sh project_init.sh project_name
  ```
  This initiates a project with a fixed directory structure for organizing files output by the workflow.

  2. Provide a tab-delimited annotation file `project_name_annotation.txt` in the `mint/project_name/data/` directory. It should include 9 columns:
    1. `projectID`: The name given in the call to `project_init.sh`.
    2. `sampleID`: An alphanumeric ID (e.g. from SRA, GEO, sequencing core, etc.). Typically these will be the names of the `.fastq` files.
    3. `humanID`: The corresponding human readable ID.
    4. `pulldown`: A binary value indicating whether the sample is the result of a pulldown experiment (1) or not (0).
    5. `bisulfite`: A binary value indicating whether the sample is the result of a bisulfite-conversion experiment (1) or not (0).
    6. `mc`: A binary value indicating whether the sample represents 5mC methylation.
    7. `hmc`: A binary value indicating whether the sample represents 5hmC methylation.
    8. `input`: A binary value indicating whether the sample represents an input.
    9. `group`: A binary value indicating which samples belong to one of two groups.

    Note that bisulfite-conversion experiments that represent both mC and hmC should have a 1 in each column. Input pulldowns can be matched to the pulldown (e.g. `mc=1` and `hmc=0` and `input=1`) or not (e.g. `mc=0` and `hmc=0` and `input=1`).

    An example of a pulldown experimental setup with sample-wise analysis is:
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
    GSE63743        SRR1686705      preeclamptic_1  1       0       0       0       1       0
    GSE63743        SRR1686706      preeclamptic_2  1       0       0       0       1       0
    GSE63743        SRR1686709      normal_1        1       0       0       0       1       0
    GSE63743        SRR1686710      normal_2        1       0       0       0       1       0
    ```
    An example of a hybrid experimental setup with comparison-wise analysis is below. Note the added lines for comparisons indicating that pulldown and bisulfite experiments should be compared according to groups 0 and 1.
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
	A more complicated example with multiple comparisons to be made.
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
