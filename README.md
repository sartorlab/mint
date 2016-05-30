## mint: Analysis, integration, classification, and annotation of DNA methylation and hydroxymethylation data

## Contents

* [Overview](https://github.com/sartorlab/mint#overview)
* [Installation](https://github.com/sartorlab/mint#installation)
	* [Dependencies](https://github.com/sartorlab/mint#dependencies)
	* [Getting mint](https://github.com/sartorlab/mint#getting-mint)
	* [Testing mint](https://github.com/sartorlab/mint#testing-mint)
* [Usage](https://github.com/sartorlab/mint#usage)
	* [Supported Experiments and Designs](https://github.com/sartorlab/mint#supported-experiments-and-designs)
	* [Setting up a project](https://github.com/sartorlab/mint#setting-up-a-project)
		* [Annotation file](https://github.com/sartorlab/mint#annotation-file)
		* [Instantiating a project](https://github.com/sartorlab/mint#instantiating-a-project)
	* [Configuring a project](https://github.com/sartorlab/mint#configuring-a-project)
	* [Running a project](https://github.com/sartorlab/mint#running-a-project)
* [Outputs](https://github.com/sartorlab/mint#outputs)
	* [FastQC](https://github.com/sartorlab/mint#fastqc)
	* [Trim Galore](https://github.com/sartorlab/mint#trim-galore)
	* [Bisulfite](https://github.com/sartorlab/mint#bisulfite)
		* [bismark](https://github.com/sartorlab/mint#bismark)
		* [methylSig](https://github.com/sartorlab/mint#methylsig)
	* [Pulldown](https://github.com/sartorlab/mint#pulldown)
		* [bowtie2](https://github.com/sartorlab/mint#bowtie2)
		* [Genome coverage](https://github.com/sartorlab/mint#genome-coverage)
		* [macs2](https://github.com/sartorlab/mint#macs2)
		* [PePr](https://github.com/sartorlab/mint#pepr)
	* [Classifications](https://github.com/sartorlab/mint#classifications)
		* [Simple](https://github.com/sartorlab/mint#simple-classification)
		* [Sample](https://github.com/sartorlab/mint#sample-classification)
		* [Comparison](https://github.com/sartorlab/mint#comparison-classification)
	* [Annotations and Visualizations](https://github.com/sartorlab/mint#annotations-and-visualizations)
	* [UCSC Browser track hub](https://github.com/sartorlab/mint#ucsc-browser-track-hub)

## Overview

The `mint` pipeline analyzes single-end reads coming from sequencing assays measuring DNA methylation. The pipeline analyzes reads from both bisulfite-conversion assays such as WGBS and RRBS, and from pulldown assays such as MeDIP-seq, hMeDIP-seq, and hMeSeal. Moreover, with data measuring both 5-methylcytosine (5mc) and 5-hydroxymethylcytosine (5hmc), the `mint` pipeline integrates the two data types to classify genomic regions of 5mc, 5hmc, a mixture, or neither.

The `mint` pipeline is executed with `make` and includes configurable steps for:

* Quality control (`FastQC` and `MultiQC`)
* Adapter and quality trimming (`trim_galore`)
* Alignment (`bismark` and `bowtie2`)
* Sample-wise quantification (`bismark` and `macs2`)
* Group-wise differential methylation detection (`methylSig` and `PePr`)
* Classification
	* of samples into regions of no, low, medium, or high methylation
	* of 5mc + 5hmc sample-wise integration into regions of 5mc, 5hmc, a mixture, or neither
	* of group comparisons into regions of hyper/hypo DMR or hyper/hypo DhMR
	* of 5mc + 5hmc group-wise integration into regions of hyper/hypo DMR, hyper/hypo DhMR, a mixture, or neither
* Genomic annotation and visualization of methylation quantifications and classifications (`annotatr`)
* Visualization in the UCSC Genome Browser

## Installation

### Dependencies

The `mint` pipeline is dependent on several software packages to carry out its analysis, integration, annotation, and visualization. Links to appropriate versions are included in the list below.

* [`bedtools` v2.25.0](https://github.com/arq5x/bedtools2/releases/tag/v2.25.0)
* [`bedops` v2.4.14](https://github.com/bedops/bedops/releases/tag/v2.4.14)
* [`bismark` v0.16.1](https://github.com/FelixKrueger/Bismark/releases/tag/0.16.1)
* [`bowtie2` v2.2.4](https://github.com/BenLangmead/bowtie2/releases/tag/v2.2.4)
* [`macs2` v2.1.0.20140616](https://pypi.python.org/pypi/MACS2/2.1.0.20140616)
* [`PePr` v1.1.5](https://github.com/shawnzhangyx/PePr/releases/tag/1.1.5)
* [`R` >= v3.2.5](https://cran.r-project.org)
	* [`annotatr` v0.7.3](https://github.com/rcavalcante/annotatr/releases/tag/v0.7.3)
	* [`methylSig` v0.4.3](https://github.com/sartorlab/methylSig/releases/tag/v0.4.3)
* [`samtools` v0.1.19](https://github.com/samtools/samtools/releases/tag/0.1.19)
* [`trim_galore` v0.4.1](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.1.zip)
* [`cutadapt` v1.9.1](https://pypi.python.org/pypi/cutadapt/1.9.1)
* [`FastQC` v0.11.5](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5_source.zip)

### Getting `mint`

The easiest way to install `mint` is to `cd` into the directory you'd like to place the `mint` folder and do:

```{bash}
git clone https://github.com/sartorlab/mint.git
```

Cloning `mint` makes getting updates easy with `git pull` from within `mint/`. Updating `mint` will never erase existing analyses.

### Testing `mint`

After installing the dependencies and cloning `mint`, there are two datasets totalling `996M` available for testing at [LINK COMING](http://sartorlab.ccmb.med.umich.edu/software/). The tests consist of one hybrid experiment (based on [GSE52945](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52945)) and one pulldown experiment (based on [GSE63743](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63743)). The examples in the Usage section below are for these two test datasets.

## Usage

### Supported Experiments and Designs

The `mint` pipeline supports **either**:

* **a hybrid setup**, with data from bisulfite-conversion experiments representing 5mc + 5hmc methylation (RRBS, WGBS, etc.) and from pulldown experiments representing 5hmc methylation (hMeDIP-seq, hMeSeal, etc.), **or**
* **a pulldown setup**, with data representing 5mc methylation (MeDIP-seq, etc.) and 5hmc methylation (hMeDIP-seq, hMeSeal, etc.).

And the `mint` pipeline can analyze:

* **sample-wise**, where no groups are present for comparison, and 5mc / 5hmc integration is done per-sample, **and/or**
* **comparison-wise**, where two or more groups are tested for differential methylation, and 5mc / 5hmc differential methylation integration is done per-group.

### Setting up a project

After setting up the [dependencies](https://github.com/sartorlab/mint#dependencies) and [getting `mint`](https://github.com/sartorlab/mint#getting-mint), you're ready to start a project. In the two examples below, we'll use the [test data](https://github.com/sartorlab/mint#testing-mint) to setup and run a project. To set up a project in `mint`, the pipeline needs the following three things:

1. A project annotation file describing the samples, and any group comparisons to make.
2. The location of the raw sequencing reads corresponding to the samples.
3. The specific parameters desired for the tools in the pipeline.

#### Annotation file

The annotation file is a tab-delimited text file named `projectID_annotation.txt`, and contains the following 9 columns (with headers):

1. `projectID`: The name of the project.
2. `sampleID`: Often an alphanumeric ID (perhaps from SRA, GEO, a sequencing core, etc.). These will be the names of the `.fastq.gz` files containing the raw reads.
3. `humanID`: An ID that makes the `sampleID`s easier to understand. For example, instead of a `sampleID` from the SRA like 'SRA876549', use 'IDH2mut_1' to indicate what kind of mutant and what replicate the sample is. **NOTE:** The `humanID` may not be unique if a particular condition has a corresponding input sample (e.g., pulldown) or a corresponding sample on a different platform (e.g., WGBS). *See examples below*.
4. `pulldown`: A binary indicating whether the sample is the result of a pulldown experiment (1) or not (0).
5. `bisulfite`: A binary indicating whether the sample is the result of a bisulfite-conversion experiment (1) or not (0).
6. `mc`: A binary indicating whether the sample represents 5mc methylation (1) or not (0).
7. `hmc`: A binary indicating whether the sample represents 5hmc methylation (1) or not (0).
8. `input`: A binary indicating whether the sample represents an input (1) or not (0).
9. `group`: A comma-separated string of integers where each integer collects all samples belonging to a particular group. **NOTE:** When testing for differential methylation in `PePr`, the group with the higher number will be `chip1` and the lower number will be `chip2`. When testing for differential methylation in `methylSig`, the group with the higher number will be `group 1` and the group with lower number will be `group 0`. Hypermethylation is then higher methylation in `group 1`, and hypomethylation is then lower methylation in `group 1`.

Additional columns can be included, but their heading cannot be one of the above 9 headings.

In particular, for the hybrid test data we have:

```{bash}
projectID	sampleID	humanID	pulldown	bisulfite	mc	hmc	input	group
test_hybrid	IDH2mut_1_hmc_pulldown	IDH2mut_1	1	0	0	1	0	1
test_hybrid	IDH2mut_2_hmc_pulldown	IDH2mut_2	1	0	0	1	0	1
test_hybrid	IDH2mut_1_hmc_input_pulldown	IDH2mut_1	1	0	0	1	1	1
test_hybrid	IDH2mut_2_hmc_input_pulldown	IDH2mut_2	1	0	0	1	1	1
test_hybrid	IDH2mut_1_mc_hmc_bisulfite	IDH2mut_1	0	1	1	1	0	1
test_hybrid	IDH2mut_2_mc_hmc_bisulfite	IDH2mut_2	0	1	1	1	0	1
test_hybrid	NBM_1_hmc_pulldown	NBM_1	1	0	0	1	0	0
test_hybrid	NBM_2_hmc_pulldown	NBM_2	1	0	0	1	0	0
test_hybrid	NBM_1_hmc_input_pulldown	NBM_1	1	0	0	1	1	0
test_hybrid	NBM_2_hmc_input_pulldown	NBM_2	1	0	0	1	1	0
test_hybrid	NBM_2_mc_hmc_bisulfite	NBM_2	0	1	1	1	0	0
test_hybrid	NBM_1_mc_hmc_bisulfite	NBM_1	0	1	1	1	0	0
test_hybrid	comparison	IDH2mut_v_NBM	1	0	0	1	0	0,1
test_hybrid	comparison	IDH2mut_v_NBM	0	1	1	1	0	0,1
```

And for the pulldown test data we have:

```{bash}
projectID	sampleID	humanID	pulldown	bisulfite	mc	hmc	input	group
test_pulldown	preeclamptic_1_hmc_pulldown	preeclamptic_1	1	0	0	1	0	1
test_pulldown	preeclamptic_2_hmc_pulldown	preeclamptic_2	1	0	0	1	0	1
test_pulldown	normal_1_hmc_pulldown	normal_1	1	0	0	1	0	0
test_pulldown	normal_2_hmc_pulldown	normal_2	1	0	0	1	0	0
test_pulldown	preeclamptic_1_mc_pulldown	preeclamptic_1	1	0	1	0	0	1
test_pulldown	preeclamptic_2_mc_pulldown	preeclamptic_2	1	0	1	0	0	1
test_pulldown	normal_1_mc_pulldown	normal_1	1	0	1	0	0	0
test_pulldown	normal_2_mc_pulldown	normal_2	1	0	1	0	0	0
test_pulldown	preeclamptic_1_input_pulldown	preeclamptic_1	1	0	0	1	1	1
test_pulldown	preeclamptic_2_input_pulldown	preeclamptic_2	1	0	0	1	1	1
test_pulldown	normal_1_input_pulldown	normal_1	1	0	0	1	1	0
test_pulldown	normal_2_input_pulldown	normal_2	1	0	0	1	1	0
test_pulldown	preeclamptic_1_input_pulldown	preeclamptic_1	1	0	1	0	1	1
test_pulldown	preeclamptic_2_input_pulldown	preeclamptic_2	1	0	1	0	1	1
test_pulldown	normal_1_input_pulldown	normal_1	1	0	1	0	1	0
test_pulldown	normal_2_input_pulldown	normal_2	1	0	1	0	1	0
test_pulldown	comparison	preeclamptic_v_normal	1	0	1	0	0	0,1
test_pulldown	comparison	preeclamptic_v_normal	1	0	0	1	0	0,1
```

#### Instantiating a project

After creating an appropriate annotation file for your project, within the `mint/` folder do the following:

1. `mkdir projects`
2. Put the `projectID_annotation.txt` file in `mint/projects/`.
3. In `mint/` do:
```{bash}
Rscript init.R --project projectID --genome genome --datapath /path/to/data
```

The `init.R` script creates an appropriate directory structure in `mint/projects/projectID/`, creates symlinks to the `.fastq.gz` files in `/path/to/data`, and creates the `makefile` and `config.mk` files that control the analysis of your project.

#### Configuring a project

The `mint/projects/projectID/config.mk` file contains options for analysis implied in the project annotation file.

Below is an example of the `config.mk` file for the `test_hybrid` data. In this file, the `PROJECT` and `GENOME` variables have been automatically populated from the `Rscript init.R` call.

At minimum, paths to reference genome information must be provided in the `Genome paths` block. Optionally, specific paths to any of the software used for analysis can be given. Options for software used and cutoffs for classification are considered reasonable defaults, but these can be changed as desired.

**NOTE:** The documentation for each program should be consulted before adding or removing parameters.

```{make}
# Configuration for mint pipeline analyses

################################################################################
# Project and experimental information

# Project name
PROJECT = test_hybrid
# Genome
GENOME = hg19

################################################################################
# Genome paths
# Genomes, like as downloaded from iGenomes, will also contain bisulfite
# converted genomes for
GENOME_PATH := /path/to/genome
# Location of bowtie2 indexes
BOWTIE2_GENOME_PATH := /path/to/bowtie2/indices
# Location of chromosome length file
CHROM_PATH := /path/to/chromosome/length/file

################################################################################
# Path to tools

PATH_TO_FASTQC := $(shell which fastqc)
PATH_TO_TRIMGALORE := $(shell which trim_galore)
PATH_TO_BISMARK := $(shell which bismark)
PATH_TO_EXTRACTOR := $(shell which bismark_methylation_extractor)
PATH_TO_BOWTIE2 := $(shell which bowtie2)
PATH_TO_R := $(shell which Rscript)
PATH_TO_MACS := $(shell which macs2)
PATH_TO_PEPR := $(shell which PePr)
PATH_TO_SAMTOOLS := $(shell which samtools)
PATH_TO_BEDTOOLS := $(shell which bedtools)
PATH_TO_BEDOPS := $(shell which bedops)
PATH_TO_AWK := $(shell which awk)
PATH_TO_BDG2BW := $(shell which bedGraphToBigWig)
PATH_TO_BDG2BB := $(shell which bedToBigBed)

################################################################################
# Command line options for tools

# FastQC
OPTS_FASTQC = --format fastq --noextract
# trim_galore bisulfite
OPTS_TRIMGALORE_BISULFITE = --quality 20 --illumina --stringency 6 -e 0.2 --gzip --length 25 --rrbs
# trim_galore pulldown
OPTS_TRIMGALORE_PULLDOWN = --quality 20 --illumina --stringency 6 -e 0.2 --gzip --length 25
# bismark
OPTS_BISMARK = --bowtie2 $(GENOME_PATH)
# bismark_methylation_extractor
OPTS_EXTRACTOR = --single-end --gzip --bedGraph --cutoff 5 --cytosine_report --genome_folder $(GENOME_PATH) --multicore 1
# bowtie2
OPTS_BOWTIE2 = -q -x $(BOWTIE2_GENOME_PATH) -U
# macs2
OPTS_MACS = -t $bowtie2Bam -c $bowtie2InputBam -f BAM -g hs --outdir ./analysis/macs_peaks -n $macsPrefix

################################################################################
# Command line options for methylSig and compare classifications

# DMC for CpG resolution, and DMR for region resolution (window size parameter
# used in the methylSig options below).
OPT_DM_TYPE = DMR

# Thresholds to use for DMCs or DMRs (above) in methylSig
OPT_MSIG_DM_FDR_THRESHOLD = 0.05
OPT_MSIG_DM_DIFF_THRESHOLD = 10

################################################################################
# Comparison specific options
# See ?methylSig::methylSigReadData and ?methylSig::methylSigCalc after installing methylSig in R for parameter information
OPTS_METHYLSIG_IDH2mut_v_NBM_mc_hmc_bisulfite = --context CpG --resolution base --destranded TRUE --maxcount 500 --mincount 5 --filterSNPs TRUE --dmtype $(OPT_DM_TYPE) --winsize.tile 50 --dispersion both --local.disp FALSE --winsize.disp 200 --local.meth FALSE --winsize.meth 200 --minpergroup 2,2 --T.approx TRUE --ncores 4 --quiet FALSE

OPTS_PEPR_IDH2mut_v_NBM_hmc_pulldown = --file-format=bam --peaktype=sharp --diff --threshold=1e-05 --num-processors=1
```

#### Running a project

Once a project is instantiated and configured, the analysis steps can begin. The `mint` pipeline is controlled by `make` and there are up to 7 simple commands to run the modules, depending on the analyses implied by the project annotation file. In general, the bisulfite and pulldown commands are independent of each other until classification, and the `*_sample` and `*_compare` commands rely on the corresponding `*_align` steps.

For example, the following commands will analyze the test hybrid data:

```{bash}
cd /path/to/mint/projects/test_hybrid
make pulldown_align
make bisulfite_align
make pulldown_sample
make bisulfite_compare
make pulldown_compare
make sample_classification
make compare_classification
```

And the following commands will analyze the test pulldown data:

```{bash}
cd /path/to/mint/projects/test_pulldown
make pulldown_align
make pulldown_sample
make pulldown_compare
make sample_classification
make compare_classification
```

Depending on the computing hardware used, projects can be run with the `make -j n` command where `n` is a positive integer. The `-j` flag specifies how many commands `make` is allowed to run simultaneously. When it is not present, the default is to run commands in serial.

**NOTE:** Some software in the `mint` pipeline have options for the number of processors to use, so some care should be taken not to exceed the computing limitations of the hardware. For example, `bismark` tends to use 5 cores per instantiation, so using `make -j 5` would really use 25 cores.

## Outputs

In the case of the `test_hybrid` project the following directory structure is created by `Rscript init.R` and the outputs of the `make` commands above are organized within:

```{bash}
test_hybrid
test_hybrid/bisulfite
test_hybrid/bisulfite/trim_fastqcs
test_hybrid/bisulfite/methylsig_calls
test_hybrid/bisulfite/raw_fastqcs
test_hybrid/bisulfite/trim_fastqs
test_hybrid/bisulfite/raw_fastqs
test_hybrid/bisulfite/bismark
test_hybrid/config.mk
test_hybrid/test_hybrid_hub
test_hybrid/test_hybrid_hub/test_hybrid_hub.html
test_hybrid/test_hybrid_hub/hub.txt
test_hybrid/test_hybrid_hub/genomes.txt
test_hybrid/test_hybrid_hub/hg19
test_hybrid/RData
test_hybrid/summary
test_hybrid/summary/tables
test_hybrid/summary/figures
test_hybrid/summary/reports
test_hybrid/makefile
test_hybrid/narrowPeak.as
test_hybrid/pulldown
test_hybrid/pulldown/bowtie2_bams
test_hybrid/pulldown/trim_fastqcs
test_hybrid/pulldown/raw_fastqcs
test_hybrid/pulldown/trim_fastqs
test_hybrid/pulldown/raw_fastqs
test_hybrid/pulldown/pepr_peaks
test_hybrid/pulldown/macs2_peaks
test_hybrid/pulldown/pulldown_coverages
test_hybrid/data
test_hybrid/data/test_hybrid_annotation.txt
test_hybrid/data/raw_fastqs
test_hybrid/classifications
test_hybrid/classifications/sample
test_hybrid/classifications/comparison
test_hybrid/classifications/simple
```

### FastQC

The `FastQC` output is by default not extracted, so the `.zip` files and `.html` files for the raw reads are located in `test_hybrid/bisulfite/raw_fastqcs` and `test_hybrid/pulldown/raw_fastqcs`. The output `FastQC` done after trimming are located in `test_hybrid/bisulfite/trim_fastqcs` and `test_hybrid/pulldown/trim_fastqcs`. The output is standard.

### Trim Galore

The adapter and quality trimmed raw reads, as well as trimming reports, are output in `test_hybrid/bisulfite/trim_fastqs` and `test_hybrid/pulldown/trim_fastqs`. The trimming reports are the normal output of `trim_galore`.

### Bisulfite

#### bismark

The results of the `bismark` alignment and the `bismark_methylation_extractor` go in `test_hybrid/bisulfite/bismark`.

#### methylSig

The results of the tests for differential methylation with `methylSig` go in `test_hybrid/bisulfite/methylsig_calls`. The `R` session image for each `methylSig` run is saved as an `.RData` file in `test_hybrid/RData`. This allows for retesting for differential methylation with different parameters without having to read in all the data again.

### Pulldown

#### bowtie2

The results of `bowtie2` alignments go in `test_hybrid/pulldown/bowtie2_bams`. All alignments are sorted and indexed after they are aligned, and the mapping efficiencies are output to a text file in the same folder.

#### Genome coverage

The read pileups using `bedtools genomecov` go in `test_hybrid/pulldown/pulldown_coverages`. This includes both pulldowns with an antibody and input pulldowns. The corresponding coverage `.bedGraph` files are compressed into `.bigWig`s and placed in `test_hybrid/test_hybrid_hub/hg19`. For downstream use, a 'merged' coverage `.bed` file is created that fills coverage gaps of up to 20bp. These files are used to determine where no signal is present for classifications.

#### macs2

The `.narrowPeak` files resulting from `macs2` peak calling go in `test_hybrid/pulldown/macs2_peaks`. Additionally, the `.pdf` images of the model used for peak calling are in the same folder.

#### PePr

The `*_chip1_peaks.bed` and `*_chip2_peaks.bed` files resulting from `PePr`'s test for differentially methylated regions go in `test_hybrid/pulldown/pepr_peaks`.

### Classifications

#### Simple classification

Simple classifications are done on the results of `bismark_methylation_extractor` and `macs2` in order to classify sites/regions in each sample as having no, low, medium, or high methylation. In the case of `bismark_methylation_extractor` we have absolute quantification of methylation levels and the breakdown classification is:

* None: [0%, 5%)
* Low: [5%, 33%)
* Medium: [33%, 66%)
* High : [66%, 100%]

Since `macs2` peaks determines qualitative methylation, we determine if a peak represents low, medium, or high methylation by dividing the fold-changes reported by `macs2` into tertiles based on the range of the fold-changes: The minimum observed fold-change and the minimum observed fold-change among the top 1%-tile of fold changes. Regions with no peaks are considered to have no methylation.

The resulting simple classifications go in `test_hybrid/classifications/simple` as `.bed` files which include colors for visualization in the UCSC Genome Browser. The `.bed` files are then compressed into `bigBed` and are located in `test_hybrid/test_hybrid_hub/hg19`.

#### Sample classification

Sample classifications are done when data is available measuring 5mc + 5hmc (WGBS, RRBS, etc.) and 5hmc (hMeDIP-seq, hMeSeal, etc.), or 5mc (MeDIP-seq, etc.) and 5hmc (hMeDIP-seq, etc.). The inputs into the classifier are sites/regions of 5mc + 5hmc and regions of 5hmc or regions of 5mc and regions of 5hmc. The tables below determine how the data is integrated into classifications.

Hybrid sample classification:

|					| hmc peak 		| No hmc peak	| No signal 		|
| :---:				| :----------:	| :--------:	| :------:			|
| **High mc+hmc**	| hmc			| mc			| hmc or mc			|
| **Low mc+hmc**	| hmc			| mc (low)		| hmc or mc (low)	|
| **No mc+hmc**		| hmc			| No methylation| No methylation	|
| **No signal**		| hmc			| No methylation| Unclassifiable	|

Pulldown sample classification:

|					| hmc peak 		| No hmc peak	| No signal 	|
| :---:				| :----------:	| :--------:	| :------:		|
| **mc peak**		| hmc and mc	| mc			| mc			|
| **No mc peak**	| hmc			| No methylation| No methylation|
| **No signal**		| hmc			| No methylation| Unlassifiable	|

The sample classifications go in `test_hybrid/classifications/sample` as `.bed` files which include colors for visualization in the UCSC Genome Browser. The `.bed` files are then compressed into `bigBed` and are located in `test_hybrid/test_hybrid_hub/hg19`.

#### Comparison classification

Comparison classifications are done when data is available measuring 5mc + 5hmc (WGBS, RRBS, etc.) and 5hmc (hMeDIP-seq, hMeSeal, etc.), or 5mc (MeDIP-seq, etc.) and 5hmc (hMeDIP-seq, etc.), and there are groups of samples compared against each other for differential methylation (DM). The inputs into the classifier are sites/regions of 5mc + 5hmc DM and regions of 5hmc DM or regions of 5mc DM and 5hmc DM. The tables below determine how the data is integrated into classifications.

Hybrid comparison classification:

|					| Hyper hmc 			| Hypo hmc				| No DM 		| No signal		|
| :------:			| :----------:			| :--------:			| :------:		| :-----:		|
| **Hyper mc + hmc**| Hyper mc & Hyper hmc	| Hyper mc & Hypo hmc	| Hyper mc		| Hyper mc		|
| **Hypo mc + hmc**	| Hypo mc & Hyper hmc	| Hypo mc & Hypo hmc	| Hypo mc		| Hypo mc		|
| **No DM**			| Hyper hmc				| Hypo hmc				| No DM			| No DM			|
| **No signal**		| Hyper hmc				| Hypo hmc				| No DM			| Unclassifiable|

Pulldown comparison classification:

|				| Hyper hmc 			| Hypo hmc				| No DM 		| No signal		|
| :------:		| :----------:			| :--------:			| :------:		| :-----:		|
| **Hyper mc**	| Hyper mc & Hyper hmc	| Hyper mc & Hypo hmc	| Hyper mc		| Hyper mc		|
| **Hypo mc**	| Hypo mc & Hyper hmc	| Hypo mc & Hypo hmc	| Hypo mc		| Hypo mc		|
| **No DM**		| Hyper hmc				| Hypo hmc				| No DM			| No DM			|
| **No signal**	| Hyper hmc				| Hypo hmc				| No DM			| Unclassifiable|

The comparison classifications go in `test_hybrid/classifications/comparison` as `.bed` files which include colors for visualization in the UCSC Genome Browser. The `.bed` files are then compressed into `bigBed` and are located in `test_hybrid/test_hybrid_hub/hg19`.

### Annotations and Visualizations

The results of `bismark_methylation_extractor`, `methylSig`, `PePr`, and the classifications are annotated to genomic regions using `annotatr`, a fast and flexible `R` package that annotates genomic regions to genomic annotations. As with `methylSig`, every `annotatr` session is saved as an `.RData` file in `test_hybrid/RData` to enable users to quickly go back to the annotations, investigate further, alter plots, or create new plots.

Only the hg19, hg38, mm9, or mm10 genomes are currently supported for annotation in the `mint` pipeline. We plan to implement the use of custom annotations *within the pipeline* in the future (`annotatr` already supports custom annotations). All files are annotated against CpG features (islands, shores, shelves, and inter CGI) and genic features based on UCSC knownGene transcripts (1-5kb upstream of promoter, promoter (<1kb upstream of TSS), 5'UTR, exons, introns, and 3'UTR). In the case of hg19, the FANTOM5 robust enhancers are also included. Annotations include UCSC transcript IDs, Entrez Gene IDs, and gene symbols.

All annotation sessions output a table of all genomic annotations intersecting the input regions in `test_hybrid/summary/tables`. Also output are summary tables indicating the number of regions annotated to a particular type, along with corresponding numbers for a set of random regions.

Additionally, a variety of plots are output to help interpret the output of `bismark_methylation_extractor`, `methylSig`, `PePr`, and the classifications.

#### Plot: Number of regions per annotation

![Regions per annotation type](https://github.com/sartorlab/mint/blob/master/docs/IDH2mut_v_NBM_hmc_pulldown_PePr_counts.png)

#### Plot: Distribution of % methylation in annotations



#### Plot: Distribution of coverage in annotations



#### Plot: Distribution of peak widths

![Distribution of peak widths](https://github.com/sartorlab/mint/blob/master/docs/IDH2mut_1_hmc_pulldown_simple_class_peak_widths.png)

#### Plot: Volcano plots of DM



#### Plot: Distribution of chip1/chip2 peaks in annotations

![Counts of DM type in annots](https://github.com/sartorlab/mint/blob/master/docs/IDH2mut_v_NBM_hmc_pulldown_PePr_cat_count_genes.png)
![Prop. of DM type in annots](https://github.com/sartorlab/mint/blob/master/docs/IDH2mut_v_NBM_hmc_pulldown_PePr_cat_prop_genes.png)

#### Plot: Distribution of classifications in annotations



For details on the general features of `annotatr`, visit the [GitHub repository](https://github.com/rcavalcante/annotatr).

### UCSC Browser track hub

Each `mint` project has a UCSC Genome Browser track hub autogenerated in `test_hybrid/test_hybrid_hub`. In order to view it on the genome browser, it must be placed in a web-facing folder and the URL should be given under the `My Hubs` tab.

Included in each track hub, where applicable, are:

* Pulldown coverage pileups from `bedtools genomecov`
* Percent methylation tracks from `bismark_methylation_extractor`
* Simple classification tracks
* `macs2` peaks
* Combined peaks from `PePr` with chip1/chip2 labels
* Differential methylation from `methylSig`
* Sample classifications
* Comparison classifications
