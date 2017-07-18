## mint: Analysis, integration, classification, and annotation of DNA methylation and hydroxymethylation data

v0.4.0

## Contents

* [Overview](#overview)
* [Installation](#installation)
	* [`mint` and Dependencies](#mint-and-dependencies)
	* [Checking Dependency Installation](#checking-dependency-installation)
	* [Reference Genomes](#reference-genomes)
	  * [Symlinks](#symlinks)
	  * [Bisulfite-converted Reference Genome](#bisulfite-converted-reference-genome)
* [Testing `mint`](#testing-mint)
* [Upgrading `mint`](#upgrading-mint)
* [Details](#details)
	* [Supported Experiments and Designs](#supported-experiments-and-designs)
	* [Setting up a project](#setting-up-a-project)
		* [Annotation file](#annotation-file)
		* [Group file](#group-file)
		* [Instantiating a project](#instantiating-a-project)
		* [Configuring a project](#configuring-a-project)
	* [Running a project](#running-a-project)
* [Outputs](#outputs)
	* [FastQC](#fastqc)
	* [Trim Galore](#trim-galore)
	* [Bisulfite](#bisulfite)
		* [bismark](#bismark)
		* [DSS](#dss)
	* [Pulldown](#pulldown)
		* [bowtie2](#bowtie2)
		* [Genome coverage](#genome-coverage)
		* [macs2](#macs2)
		* [csaw](#csaw)
	* [Classifications](#classifications)
		* [Simple](#simple-classification)
		* [Sample](#sample-classification)
		* [Comparison](#comparison-classification)
	* [Annotations and Visualizations](#annotations-and-visualizations)
	* [UCSC Browser track hub](#ucsc-browser-track-hub)

## Overview

The `mint` pipeline analyzes single-end reads coming from sequencing assays measuring DNA methylation and hydroxymethylation. The pipeline analyzes reads from both bisulfite-converted assays such as WGBS and RRBS, and from pulldown assays such as MeDIP-seq, hMeDIP-seq, and hMeSeal. Moreover, with data measuring both 5-methylcytosine (5mc) and 5-hydroxymethylcytosine (5hmc), the `mint` pipeline integrates the two data types to classify genomic regions of 5mc, 5hmc, a mixture, or neither.

The `mint` pipeline is executed with `make` and includes configurable steps for:

* Quality control (`FastQC` and `MultiQC`)
* Adapter and quality trimming (`trim_galore`)
* Alignment (`bismark` and `bowtie2`)
* Sample-wise quantification (`bismark` and `macs2`)
* Differential methylation detection with multi-factor models with categorical or continuous covariates (`DSS` and `csaw`)
* Classification
	* of samples into regions of no, low, medium, or high methylation
	* of 5mc + 5hmc sample-wise integration into regions of 5mc, 5hmc, a mixture, or neither
	* of group comparisons into regions of hyper/hypo DMR or hyper/hypo DhMR
	* of 5mc + 5hmc group-wise integration into regions of hyper/hypo DMR, hyper/hypo DhMR, a mixture, or neither
* Genomic annotation and visualization of methylation quantifications and classifications (`annotatr`)
* Visualization in the UCSC Genome Browser

The mint pipeline is also implemented for [Galaxy](https://usegalaxy.org), and the repository for setting up the Galaxy version of mint is <https://github.com/sartorlab/mint_galaxy>.

[Top](#contents)

## Installation

### `mint` and Dependencies

The following steps will install `mint` and its dependencies on a Linux or macOS system. **NOTE:** Any references to `/path/to` need to be modified throughout the code below.

1. `cd` into the directory you'd like `mint` installed.
2. Clone the `mint` repository with:

  ```
  git clone https://github.com/sartorlab/mint.git
  ```
3. `cd` into the `mint` directory and install the dependencies to `mint/apps` with:

  ```
  bash install/install_deps.sh
  ```
  Some NOTES on the dependencies:
  * A complete list of dependencies with links to corresponding versions can be found in [VERSIONS.md](https://github.com/sartorlab/mint/blob/master/VERSIONS.md).
  * All dependencies are installed in `mint/apps` and should not interfere with any preexisting installs of the dependencies on the system.
  * Users may already have some dependencies installed. In which case, the `install/install_deps.sh` script can be modified to suit those situations.
  * The install script checks for a `pip` installation and installs it with the `--user` flag if it is not present.
  * Python packages are installed in a `virtualenv` in `mint/apps/mint_env`, ensuring existing Python packages remain untouched.
  * The `virtualenv` activation binary is symlinked in `/mint/apps/bin` as `mint_venv` for convenience.
  * The install script expects a version of `R >= 3.3.0` to be installed on the system.
4. Add `/path/to/mint/apps/bin` to your `$PATH` variable in `~/.bashrc` or `~/.bash_profile` (whichever exists), by adding the corrected version of this line:

  ```
  export PATH=/path/to/mint/apps/bin:~/.local/bin:$PATH
  ```
  NOTE: Adding the `/path/to/mint/apps/bin` *before* `$PATH` ensures that the versions installed in `mint/apps` are used rather than previous installations.
5. Source the `~/.bashrc` or `~/.bash_profile` to update your path in your current terminal session.

  ```
  source ~/.bashrc
  ```
6. In `/path/to/mint` do:

  ```
  mkdir projects
  ```

[Top](#contents)

### Checking Dependency Installation

At this point using `which` with one of the dependencies should return a path in `mint/apps`:
```
which bowtie2
# Should return /path/to/mint/apps/bin/bowtie2
```

It is also worthwhile to try loading the `virtualenv` and testing `which` on a Python package:
```
# Activate the virtualenv
# NOTE: Running install/install_deps.sh creates a symlink for mint/apps/mint_env/bin/activate
# as mint/apps/mint_venv
source mint_venv

which multiqc
# Should return /path/to/mint/apps/mint_env/bin/multiqc
```

[Top](#contents)

### Reference Genomes

As with any high-throughput sequencing analysis, the correct reference genome is required to align reads. The most convenient resource for reference genomes is the [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html) page.

[Top](#contents)

#### Symlinks

The following symlinks at the root of the reference genome folder are needed for the `config.mk` file `mint` creates for each project.
```
# Go to the reference genome root
cd /path/to/reference/genome/

ln -s ./Sequence/WholeGenomeFasta/Bisulfite_Genome Bisulfite_Genome
ln -s ./Annotation/Genes/genes.gtf genes.gtf
ln -s ./Sequence/Bowtie2Index/genome.1.bt2 genome.1.bt2
ln -s ./Sequence/Bowtie2Index/genome.2.bt2 genome.2.bt2
ln -s ./Sequence/Bowtie2Index/genome.3.bt2 genome.3.bt2
ln -s ./Sequence/Bowtie2Index/genome.4.bt2 genome.4.bt2
ln -s ./Sequence/Bowtie2Index/genome.fa genome.fa
ln -s ./Sequence/Bowtie2Index/genome.fa.fai genome.fa.fai
ln -s ./Sequence/Bowtie2Index/genome.rev.1.bt2 genome.rev.1.bt2
ln -s ./Sequence/Bowtie2Index/genome.rev.2.bt2 genome.rev.2.bt2

# NOTE: mint requires a chromosome length file for many operations this particular
# line of code is specific to the hg19 reference genome from iGenomes
ln -s ./Annotation/Archives/archive-2013-03-06-11-23-03/Genes/ChromInfo.txt chromInfo_hg19.txt
```

[Top](#contents)

#### Bisulfite-converted Reference Genome

If you download a reference genome from iGenomes, it is possible that a bisulfite-converted reference genome is already included. If not, in order to use the [`bismark`](http://www.bioinformatics.babraham.ac.uk/projects/bismark/) aligner on WGBS or RRBS data you must create your own. For more information [see the documentation](https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html#i-running-bismark-genome-preparation). In short the following will build it:
```
bismark_genome_preparation [options] <path_to_genome_folder>
```

[Top](#contents)

## Testing `mint`

After following the [installation instructions](#installation), we recommend testing the pipeline on a small dataset to ensure everything is working properly. We have created a repository with test data at <https://github.com/sartorlab/mint_test> based on data from [GSE52945](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52945). Follow these steps to test the entire pipeline:

1. Clone the test data repository in a folder **outside the `mint` installation**:

  ```
  git clone https://github.com/sartorlab/mint_test.git
  ```

2. Copy the metadata files for the test data into `/path/to/mint/projects`:

  ```
  cp /path/to/mint_test/hybrid/*txt /path/to/mint/projects/
  ```

3. `cd` to `mint`:

  ```
  cd /path/to/mint
  ```

4. Initialize the `test_hybrid_small` project:

  ```
  Rscript init.R --project test_hybrid_small --genome hg19 --genomepath /path/to/hg19 --chrompath /path/to/hg19/chromInfo_hg19.txt --datapath /path/to/mint_test/hybrid/
  ```

5. Run the modules. In the simplest form, the following would run in series with all output going to `stdout`:

  ```
  # Assumes an interactive server session, and NOT a cluster with queueing

  # cd into test_hybrid_small
  cd /path/to/mint/projects/test_hybrid_small

  # Activate the virtualenv
  source mint_venv

  make pulldown_align
  make bisulfite_align
  make pulldown_sample
  make bisulfite_sample
  make pulldown_compare
  make bisulfite_compare
  make compare_classification
  make sample_classification
  ```

  To help the tests run faster, you can use the `-j` flag of `make` to run operations in parallel, and use `&` to run operations in the background:

  ```
  # Assumes an interactive server session where the server has ample resources
  # For example, a server with 20 cores and 128GB RAM will be fine

  # cd into test_hybrid_small
  cd /path/to/mint/projects/test_hybrid_small

  # Activate the virtualenv
  source mint_venv

  # Here we use nohup to write the stdout to a file
  # These can be run simultaneously
  nohup make -j 4 pulldown_align > nohup_pulldown_align.out &
  nohup make -j 4 bisulfite_align > nohup_bisulfite_align.out &

  # After the above two are finished these can be run simultaneously
  nohup make -j 4 pulldown_sample > nohup_pulldown_sample.out &
  nohup make -j 4 bisulfite_sample > nohup_bisulfite_sample.out &

  # After the above two are finished these can be run simultaneously
  nohup make -j 3 pulldown_compare > nohup_pulldown_compare.out &
  nohup make -j 4 bisulfite_compare > nohup_bisulfite_compare.out &

  # After the above two are finished, these should be run serially
  nohup make -j 3 compare_classification > nohup_compare_classification.out &
  nohup make -j 4 sample_classification > nohup_sample_classification.out &
  ```

As documented in <https://github.com/sartorlab/mint_test>, the test data contains reads drawn from the following regions, and have the following resulting features:

* chr1:2209176-2214175 (hyper_mc_hyper_hmc)
* chr2:71942840-72005839 (hypo_mc_hyper_hmc)
* chr7:36074781-36102780 (hyper_hmc)
* chr8:35090385-35095384 (hyper_mc)
* chr11:67384075-67389174 (hypo_mc_hypo_hmc)
* chr16:3046564-3096563 (hyper_mc_hyper_hmc and hyper_mc_hyper_hmc)
* chr19:39702429-39752428 (hyper_mc_hypo_hmc)

[Top](#contents)

## Upgrading `mint`

To upgrade `mint`, all you need to do is `cd /path/to/mint` and do `git pull`. Existing analyses in the `mint/projects/` folder will be unaffected (since this folder isn't tracked by `git`). To use the upgraded version of `mint` on existing projects, you will need to reinitialize a project with `Rscript init.R ...`.

If there are changes to any directories that *are* tracked by `git`, you might need to `git stash` your changes, and then `git pull`.

## Details

### Supported Experiments and Designs

The `mint` pipeline supports **either**:

* **a hybrid setup**, with data from bisulfite-conversion experiments representing 5mc + 5hmc methylation (RRBS, WGBS, etc.) and from pulldown experiments representing 5hmc methylation (hMeDIP-seq, hMeSeal, etc.), **or**
* **a pulldown setup**, with data representing 5mc methylation (MeDIP-seq, etc.) and 5hmc methylation (hMeDIP-seq, hMeSeal, etc.).

And the `mint` pipeline can analyze:

* **sample-wise**, where no groups are present for comparison, and 5mc / 5hmc integration is done per-sample, **and/or**
* **comparison-wise**, where groups are tested for differential methylation using a multi-factor design with categorical and/or continuous covariates, and 5mc / 5hmc differential methylation integration is done per-comparison.

[Top](#contents)

### Setting up a project

Taking the example in [Testing `mint`](#testing-mint) as a guide, setting up a project requires the following:

1. A [samples file](#samples-file) describing the samples and any covariates associated with them. Put this in `/path/to/mint/projects`.
2. A [comparisons file](#comparisons-file) describing the comparisons to be made and models and covariates to use. Put this in `/path/to/mint/projects`.
3. The location of the raw sequencing reads corresponding to the samples, and the relevant reference genome.

[Top](#contents)

#### Samples file

The sample file is a tab-delimited text file named `[projectID]_samples.txt` in `/path/to/mint/projects`, and contains the following 9+ columns (with headers):

1. `projectID`: The name of the project.
2. `sampleID`: Often an alphanumeric ID (perhaps from SRA, GEO, a sequencing core, etc.). These will be the names of the `.fastq.gz` files containing the raw reads.
3. `humanID`: An ID that makes the `sampleID`s easier to understand. For example, instead of a `sampleID` from the SRA like 'SRA876549', use 'IDH2mut_1' to indicate what kind of mutant and what replicate the sample is. **NOTE:** The `humanID` may not be unique if a particular condition has a corresponding input sample (e.g., pulldown) or a corresponding sample on a different platform (e.g., WGBS). *See examples below*.
4. `pulldown`: A binary indicating whether the sample is the result of a pulldown experiment (1) or not (0).
5. `bisulfite`: A binary indicating whether the sample is the result of a bisulfite-conversion experiment (1) or not (0).
6. `mc`: A binary indicating whether the sample represents 5mc methylation (1) or not (0). If a sample was run on WGBS, this column and the `hmc` column would both be 1.
7. `hmc`: A binary indicating whether the sample represents 5hmc methylation (1) or not (0). If a sample was run on WGBS, this column and the `mc` column would both be 1.
8. `input`: A binary indicating whether the sample represents an input (1) or not (0).
9. `group`: A comma-separated string denoting the group numbers the sample belongs to for comparisons.
10. Any additional columns are categorical or continuous covariates, and their column headers should match what is used in the models in `[projectID]_comparisons.txt`.

In particular, for the [test](#testing-mint) data the annotation file looks like:

```
projectID	sampleID	humanID	pulldown	bisulfite	mc	hmc	input	group	subject	age
test_hybrid_small	IDH2mut_1_hmeseal	IDH2mut_1	1	0	0	1	0	1	1	3
test_hybrid_small	IDH2mut_2_hmeseal	IDH2mut_2	1	0	0	1	0	1	2	6
test_hybrid_small	IDH2mut_1_hmeseal_input	IDH2mut_1	1	0	0	1	1	1	1	3
test_hybrid_small	IDH2mut_2_hmeseal_input	IDH2mut_2	1	0	0	1	1	1	2	6
test_hybrid_small	IDH2mut_1_errbs	IDH2mut_1	0	1	1	1	0	1	1	3
test_hybrid_small	IDH2mut_2_errbs	IDH2mut_2	0	1	1	1	0	1	2	6
test_hybrid_small	NBM_1_hmeseal	NBM_1	1	0	0	1	0	0	1	10
test_hybrid_small	NBM_2_hmeseal	NBM_2	1	0	0	1	0	0	2	15
test_hybrid_small	NBM_1_hmeseal_input	NBM_1	1	0	0	1	1	0	1	10
test_hybrid_small	NBM_2_hmeseal_input	NBM_2	1	0	0	1	1	0	2	15
test_hybrid_small	NBM_1_errbs	NBM_1	0	1	1	1	0	0	1	10
test_hybrid_small	NBM_2_errbs	NBM_2	0	1	1	1	0	0	2	15
```

[Top](#contents)

#### Comparisons file

The comparisons file is a tab-delimeted file named `[projectID]_comparisons.txt` in `/path/to/mint/projects`, and contains 13 columns:

1. `projectID`: The name of the project.
2. `comparison`: The prefix for the output files.
4. `pulldown`: A binary indicating whether the comparison is for a pulldown experiment (1) or not (0).
5. `bisulfite`: A binary indicating whether the comparison is for a bisulfite-conversion experiment (1) or not (0).
6. `mc`: A binary indicating whether the comparison is for 5mc methylation (1) or not (0). If for WGBS, this column and the `hmc` column would both be 1.
7. `hmc`: A binary indicating whether the comparison is for 5hmc methylation (1) or not (0). If for WGBS, this column and the `mc` column would both be 1.
8. `input`: A logical indicating whether to use the input data in the test for differential methylation. This only applies to `csaw`.
9. `model`: A string as one would pass to `formula()`. Note, `group` should be in the model, and any covariates should have matching column headings in `[projectID]_samples.txt`.
10. `contrast`: A comma-separated binary string denoting which coefficient in the model to test.
11. `covariates`: A comma-separated string indicating the names of the covariates to use.
12. `covIsNumeric`: A comma-separated binary string denoting which covariates are numeric (1) and categorical (0).
13. `groups`: A comma-separated string denoting which groups are present in this comparison. This is used to get the correct samples from `[projectID_samples.txt]`.
14. `interpretation`: A comma-separated string with descriptions for what it means for a region to have `logFC` or `methdiff` less than 0 (first entry) and greater than 0 (second entry).

In particular, for the [test](#testing-mint) data the comparisons file looks like:

```
projectID	comparison	pulldown	bisulfite	mc	hmc	input	model	contrast	covariates	covIsNumeric	groups	interpretation
test_hybrid_small	IDH2mut_v_NBM	1	0	0	1	TRUE	~1+group	"0,1"	NA	0	"0,1"	"NBM,IDH2mut"
test_hybrid_small	IDH2mut_v_NBM	0	1	1	1	FALSE	~1+group	"0,1"	NA	0	"0,1"	"NBM,IDH2mut"
test_hybrid_small	IDH2mut_v_NBM_paired	1	0	0	1	TRUE	~1+group+subject	"0,1,0"	subject	0	"0,1"	"NBM,IDH2mut"
test_hybrid_small	IDH2mut_v_NBM_paired	0	1	1	1	FALSE	~1+group+subject	"0,1,0"	subject	0	"0,1"	"NBM,IDH2mut"
```

[Top](#contents)

#### Instantiating a project

After creating an appropriate annotation file for your project, within the `mint/` folder do the following:

1. `mkdir projects`
2. Put the `test_hybrid_small_annotation.txt` file in `mint/projects/`.
3. In `mint/` do:

```
Rscript init.R --project projectID --genome genome_build --genomepath /path/to/genome --chrompath /path/to/genome/chromInfo_file.txt --datapath /path/to/data/
```

The `init.R` script creates an appropriate directory structure in `mint/projects/projectID/`, creates symlinks to the `.fastq.gz` files in `/path/to/data`, and creates the `makefile` and `config.mk` files that control the analysis of your project.

NOTE: The genome should use the UCSC notation of `hg19, hg38, mm9`, or `mm10` in order to take advantage of genomic annotations. The `--genome` parameter should be the same genome as the `--genomepath` in spirit, if not in name. e.g. `--genome hg38` and `--genomepath /path/to/GRCh38`.

[Top](#contents)

#### Configuring a project

The `mint/projects/projectID/config.mk` file contains options for analysis implied in the project annotation file. As appropriate, links to documentation for each software tool is provided at the corresponding section fo the `config.mk` file.

**NOTE:** Default parameters may not be correct for your project, so check them carefully!

[Top](#contents)

#### Running a project

Once a project is instantiated and configured, the analysis steps can begin. The `mint` pipeline is controlled by `make` and there are up to 8 simple commands to run the modules, depending on the analyses implied by the project annotation file. In general, the bisulfite and pulldown commands are independent of each other until classification, and the `*_sample` and `*_compare` commands rely on the corresponding `*_align` steps.

For example, the following commands will analyze the test hybrid data:

```
# Go to the project root
cd /path/to/mint/projects/test_hybrid_small

# Activate the python virtualenv
source mint_venv

# Run modules
make pulldown_align
make bisulfite_align
make pulldown_sample
make bisulfite_sample
make bisulfite_compare
make pulldown_compare
make sample_classification
make compare_classification
```

To see what will be run by the pipeline without *actually* running anything, you can `make -n pulldown_align`, etc. for each of the commands. *NOTE:* This is also a way to check that a module completed without problem. If you run `make -n pulldown_align` after the module finishes, `make` will report that nothing is to be done for that rule.

Depending on the computing hardware used, projects can be run with the `make -j n` command where `n` is a positive integer. The `-j` flag specifies how many commands `make` is allowed to run simultaneously. When it is not present, the default is to run commands in serial.

**NOTE:** Some software in the `mint` pipeline have options for the number of processors to use, so some care should be taken not to exceed the computing limitations of the hardware. Tools that have parameters for the user of multiple processors are: `bismark_methylation_extractor`. By default, the number of processors to use is set to 1.

[Top](#contents)

## Outputs

In the case of the `test_hybrid_small` project the following directory structure is created by `Rscript init.R` and the outputs of the `make` commands above are organized within:

```
test_hybrid_small
test_hybrid_small/tmp
test_hybrid_small/data
test_hybrid_small/data/raw_fastqs
test_hybrid_small/data/test_hybrid_small_annotation.txt
test_hybrid_small/summary
test_hybrid_small/summary/figures
test_hybrid_small/summary/tables
test_hybrid_small/summary/reports
test_hybrid_small/test_hybrid_small_hub
test_hybrid_small/test_hybrid_small_hub/hg19
test_hybrid_small/test_hybrid_small_hub/genomes.txt
test_hybrid_small/test_hybrid_small_hub/hub.txt
test_hybrid_small/test_hybrid_small_hub/test_hybrid_small_hub.html
test_hybrid_small/classifications
test_hybrid_small/classifications/simple
test_hybrid_small/classifications/sample
test_hybrid_small/classifications/comparison
test_hybrid_small/RData
test_hybrid_small/bisulfite
test_hybrid_small/bisulfite/raw_fastqs
test_hybrid_small/bisulfite/raw_fastqcs
test_hybrid_small/bisulfite/trim_fastqs
test_hybrid_small/bisulfite/trim_fastqcs
test_hybrid_small/bisulfite/bismark
test_hybrid_small/bisulfite/dss
test_hybrid_small/pulldown
test_hybrid_small/pulldown/raw_fastqs
test_hybrid_small/pulldown/raw_fastqcs
test_hybrid_small/pulldown/trim_fastqs
test_hybrid_small/pulldown/trim_fastqcs
test_hybrid_small/pulldown/bowtie2_bams
test_hybrid_small/pulldown/pulldown_coverages
test_hybrid_small/pulldown/macs2_peaks
test_hybrid_small/pulldown/csaw
test_hybrid_small/makefile
test_hybrid_small/config.mk
test_hybrid_small/narrowPeak.as
```

[Top](#contents)

### FastQC

The `FastQC` output is by default not extracted, so the `.zip` files and `.html` files for the raw reads are located in `test_hybrid_small/bisulfite/raw_fastqcs` and `test_hybrid_small/pulldown/raw_fastqcs`. The output `FastQC` done after trimming are located in `test_hybrid_small/bisulfite/trim_fastqcs` and `test_hybrid_small/pulldown/trim_fastqcs`.

### Trim Galore

The adapter and quality trimmed raw reads, as well as trimming reports, are output in `test_hybrid_small/bisulfite/trim_fastqs` and `test_hybrid_small/pulldown/trim_fastqs`. The trimming reports are the normal output of `trim_galore`.

### Bisulfite

#### bismark

The results of the `bismark` alignment and methylation quantification from `bismark_methylation_extractor` go in `test_hybrid_small/bisulfite/bismark`. These include sorted and indexed alignment `.bam`s, methylation rate `bedGraph`s, coverage files, CpG reports, M-bias plots, and the splitting reports.

#### DSS

The results of the tests for differential methylation with `DSS` go in `test_hybrid_small/bisulfite/dss`.

### Pulldown

#### bowtie2

The results of `bowtie2` alignments go in `test_hybrid_small/pulldown/bowtie2_bams`. All alignments are sorted and indexed after they are aligned, and the mapping efficiencies are output to a text file in the same folder.

#### Genome coverage

The read pileups using `bedtools genomecov` go in `test_hybrid_small/pulldown/pulldown_coverages`. This includes both pulldowns with an antibody and input pulldowns. The corresponding coverage `.bedGraph` files are compressed into `.bigWig`s and placed in `test_hybrid_small/test_hybrid_small_hub/hg19`. For downstream use, a 'merged' coverage `.bed` file is created that fills coverage gaps of up to 20bp. These files are used to determine where no signal is present for classifications.

#### macs2

The `.narrowPeak` files resulting from `macs2` peak calling go in `test_hybrid_small/pulldown/macs2_peaks`. Additionally, the `.pdf` images of the model used for peak calling are in the same folder.

#### csaw

The `*_csaw_significant.txt` file resulting from `csaw`'s test for differentially methylated regions goes in `test_hybrid_small/pulldown/csaw`.

[Top](#contents)

### Classifications

#### Simple classification

Simple classifications are done on the results of `bismark_methylation_extractor` and `macs2` in order to classify sites/regions in each sample as having no, low, medium, or high methylation. In the case of `bismark_methylation_extractor` we have absolute quantification of methylation levels and the breakdown classification is:

* None: [0%, 5%)
* Low: [5%, 33%)
* Medium: [33%, 66%)
* High : [66%, 100%]

Since `macs2` peaks determines qualitative methylation, we determine if a peak represents low, medium, or high methylation by dividing the fold-changes reported by `macs2` into tertiles based on the range of the fold-changes: The minimum observed fold-change and the minimum observed fold-change among the top 1%-tile of fold changes. Regions with no peaks are considered to have no methylation.

The resulting simple classifications go in `test_hybrid_small/classifications/simple` as `.bed` files which include colors for visualization in the UCSC Genome Browser. The `.bed` files are then compressed into `bigBed` and are located in `test_hybrid_small/test_hybrid_small_hub/hg19`.

```
chr2	39812689	39812804	hmc_low	1000	.	39812689	39812804	102,102,255
chr21	15197751	15198018	hmc_low	1000	.	15197751	15198018	102,102,255
chr21	15198928	15199092	hmc_low	1000	.	15198928	15199092	102,102,255
chr21	15200122	15200346	hmc_low	1000	.	15200122	15200346	102,102,255
chr21	15353077	15353463	hmc_med	1000	.	15353077	15353463	0,0,255
```

[Top](#contents)

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
| **No signal**		| hmc			| No methylation| Unclassifiable	|

The sample classifications go in `test_hybrid_small/classifications/sample` as `.bed` files which include colors for visualization in the UCSC Genome Browser. The `.bed` files are then compressed into `bigBed` and are located in `test_hybrid_small/test_hybrid_small_hub/hg19`.

```
chr21   9826885 9826886 mc      1000    .       9826885 9826886 255,0,0
chr21   9826886 9826887 mc_low  1000    .       9826886 9826887 255,165,0
chr21   9826887 9826888 mc      1000    .       9826887 9826888 255,0,0
chr21   9826891 9826892 no_meth 1000    .       9826891 9826892 0,0,0
chr21   9826892 9826893 mc_low  1000    .       9826892 9826893 255,165,0
chr21   9826895 9826896 no_meth 1000    .       9826895 9826896 0,0,0
chr21   9826896 9826897 mc_low  1000    .       9826896 9826897 255,165,0
chr21   9826899 9826901 mc_low  1000    .       9826899 9826901 255,165,0
chr21   9826903 9826905 mc_low  1000    .       9826903 9826905 255,165,0
```

[Top](#contents)

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

The comparison classifications go in `test_hybrid_small/classifications/comparison` as `.bed` files which include colors for visualization in the UCSC Genome Browser. The `.bed` files are then compressed into `bigBed` and are located in `test_hybrid_small/test_hybrid_small_hub/hg19`.

```
chr21   16423450        16423650        hypo_hmc        1000    .       16423450        16423650        102,102,255
chr21   16423650        16424400        no_DM   1000    .       16423650        16424400        0,0,0
chr21   16424400        16424650        hyper_hmc       1000    .       16424400        16424650        0,0,255
chr21   16424650        16427750        no_DM   1000    .       16424650        16427750        0,0,0
chr21   16427750        16427900        hypo_hmc        1000    .       16427750        16427900        102,102,255
chr21   16427900        16429700        no_DM   1000    .       16427900        16429700        0,0,0
chr21   16429700        16429850        hyper_hmc       1000    .       16429700        16429850        0,0,255
```

[Top](#contents)

### Annotations and Visualizations

The genomic regions given by `bismark_methylation_extractor`, [`DSS`](http://www.bioconductor.org/packages/DSS), [`csaw`](http://bioconductor.org/packages/release/bioc/html/csaw.html), and the classifications are annotated to genomic annotations using [`annotatr`](https://github.com/rcavalcante/annotatr), a fast and flexible `R` package designed for exactly this task. As with `DSS`, every `annotatr` session is saved as an `.RData` file in `test_hybrid_small/RData` to enable users to quickly go back to the annotations, investigate further, alter plots, or create new plots.

Only the hg19, hg38, mm9, or mm10 genomes are currently supported for annotation in the `mint` pipeline. We plan to implement the use of custom annotations *within the pipeline* in the future (`annotatr` already supports custom annotations). All files are annotated against CpG features (islands, shores, shelves, and inter CGI) and genic features based on UCSC knownGene transcripts (1-5kb upstream of promoter, promoter (\<1kb upstream of TSS), 5'UTR, exons, introns, and 3'UTR). In the case of hg19, the FANTOM5 robust enhancers are also included. Annotations include UCSC transcript IDs, Entrez Gene IDs, and gene symbols.

All annotation sessions output a table of all genomic annotations intersecting the input regions in `test_hybrid_small/summary/tables`. Also output are summary tables indicating the number of regions annotated to a particular type, along with corresponding numbers for a set of random regions.

For details on the general features of `annotatr`, visit the [GitHub repository](https://github.com/rcavalcante/annotatr).

[Top](#contents)

#### Table: Complete annotations

```
annot_chrom	annot_start	annot_end	annot_strand	annot_type	annot_id	entrez_id	symbol	data_chrom	data_start	data_end	data_strand	name
chr21	9909278	9966321	-	hg19_knownGenes_introns	uc002zka.2,uc021wgx.1	100132288	TEKT4P2	chr21	9909661	9909789	*	hmc_low
chr21	9909515	9909759	*	hg19_cpg_islands	island:17089	NA	NA	chr21	9909661	9909789	*	hmc_low
chr21	9909759	9912040	*	hg19_cpg_shores	shore:30779	NA	NA	chr21	9909661	9909789	*	hmc_low
chr21	9909278	9966321	-	hg19_knownGenes_introns	uc002zka.2,uc021wgx.1	100132288	TEKT4P2	chr21	9944159	9944343	*	hmc_low
```

[Top](#contents)

#### Table: Annotation summaries

```
data_type       annot_type      n
Data    hg19_cpg_inter  2977
Data    hg19_cpg_islands        134
Data    hg19_cpg_shelves        327
Data    hg19_cpg_shores 519
Data    hg19_enhancers_fantom   294
Data    hg19_knownGenes_1to5kb  399
Data    hg19_knownGenes_3UTRs   131
Data    hg19_knownGenes_5UTRs   1533
Data    hg19_knownGenes_exons   408
Data    hg19_knownGenes_intergenic      1528
Data    hg19_knownGenes_introns 3562
Data    hg19_knownGenes_promoters       158
Random Regions  hg19_cpg_inter  3420
Random Regions  hg19_cpg_islands        48
Random Regions  hg19_cpg_shelves        157
Random Regions  hg19_cpg_shores 177
Random Regions  hg19_enhancers_fantom   39
Random Regions  hg19_knownGenes_1to5kb  208
Random Regions  hg19_knownGenes_3UTRs   60
Random Regions  hg19_knownGenes_5UTRs   652
Random Regions  hg19_knownGenes_exons   169
Random Regions  hg19_knownGenes_intergenic      2394
Random Regions  hg19_knownGenes_introns 1938
Random Regions  hg19_knownGenes_promoters       66
```

Additionally, a variety of plots are output to help interpret the output of `bismark_methylation_extractor`, `DSS`, `csaw`, and the classifications.

[Top](#contents)

#### Plot: Number of regions per annotation

![Regions per annotation type](https://github.com/sartorlab/mint/blob/master/docs/IDH2mut_v_NBM_hmc_pulldown_PePr_counts.png)

[Top](#contents)

#### Plot: Distribution of % methylation in annotations

![Distribution of perc. meth. in annotations](https://github.com/sartorlab/mint/blob/master/docs/IDH2mut_1_mc_hmc_bisulfite_bismark_percmeth.png)

[Top](#contents)

#### Plot: Distribution of coverage in annotations

![Distribution of coverage in annotations](https://github.com/sartorlab/mint/blob/master/docs/IDH2mut_1_mc_hmc_bisulfite_bismark_coverage.png)

[Top](#contents)

#### Plot: Distribution of peak widths

![Distribution of peak widths](https://github.com/sartorlab/mint/blob/master/docs/IDH2mut_1_hmc_pulldown_simple_class_peak_widths.png)

[Top](#contents)

#### Plot: Volcano plots of DM

![Volcano plots from DSS results](https://github.com/sartorlab/mint/blob/master/docs/IDH2mut_v_NBM_mc_hmc_bisulfite_DMR_DSS_volcano.png)

[Top](#contents)

#### Plot: Distribution of DSS calls in annotations

![Distribution of DSS calls in annots](https://github.com/sartorlab/mint/blob/master/docs/IDH2mut_v_NBM_mc_hmc_bisulfite_DMR_DSS_DMstatus_prop_genes.png)

[Top](#contents)

#### Plot: Distribution of chip1/chip2 peaks in annotations

![Counts of DM type in annots](https://github.com/sartorlab/mint/blob/master/docs/IDH2mut_v_NBM_hmc_pulldown_PePr_cat_count_genes.png)

![Prop. of DM type in annots](https://github.com/sartorlab/mint/blob/master/docs/IDH2mut_v_NBM_hmc_pulldown_PePr_cat_prop_genes.png)

[Top](#contents)

#### Plot: Distribution of classifications in annotations

![Simple classification pulldown](https://github.com/sartorlab/mint/blob/master/docs/NBM_2_hmc_pulldown_simple_class_cat_prop_cpgs.png)

![Simple classification bisulfite](https://github.com/sartorlab/mint/blob/master/docs/NBM_2_mc_hmc_bisulfite_simple_class_cat_prop_cpgs.png)

![Sample classification](https://github.com/sartorlab/mint/blob/master/docs/IDH2mut_2_sample_class_cat_prop_genes.png)

![Comparison classification](https://github.com/sartorlab/mint/blob/master/docs/IDH2mut_v_NBM_compare_class_cat_prop_genes.png)

[Top](#contents)

### UCSC Browser track hub

Each `mint` project has a UCSC Genome Browser track hub autogenerated in `test_hybrid_small/test_hybrid_small_hub`. In order to view it on the genome browser, it must be placed in a web-facing folder and the URL should be given under the `My Hubs` tab.

Included in each track hub, where applicable, are:

* Pulldown coverage pileups from `bedtools genomecov`
* Percent methylation tracks from `bismark_methylation_extractor`
* Simple classification tracks for `bismark_methylation_extractor`, `macs2`, and `csaw`
* Differential methylation from `csaw` with group labels
* Differential methylation from `DSS`
* Sample classifications
* Comparison classifications

[Top](#contents)
