# mint
## A pipeline for the analysis, integration, and visualization of DNA methylation data
—————

## Motivation
DNA methylation is known to occur in a variety of forms. Canonical 5-methylcytosine (5mC) is the best studied, and has been shown to have a variety of roles in differentiation and regulation. More recently, studies of 5-hydroxymethylcytosine (5hmC) imply it is a stable base-modification with biological roles distinct from those of 5mC. 5‐formylcytosine (5fC) and 5‐carboxycytosine (5caC) are less well known, but are increasingly being explored.

Current bisulfite-conversion + sequencing technologies (e.g. BS-seq and RRBS) are unable to distinguish between 5mC and 5hmC because both marks protect the cytosine from bisulfite-conversion. Newer technologies designed to detect only 5hmC (e.g. oxBS-seq and TAB-seq) currently have low-reproducibility and conversion efficiency. There are, however, specific antibodies capable of pulling down either 5mC or 5hmC alone. These pulldown methods, unfortunately, do not achieve base-pair resolution, and are qualitative in nature.

We have developed a __m__ethylation __int__egration pipeline using established bioinformatics tools to determine regions of 5mC and 5hmC genome-wide, and newly developed classifiers to attempt to pull apart 5mC regions versus 5hmC regions.

## Methods
Bisulfite-conversion methods (e.g. BS-seq and RRBS) quantify 5mC + 5hmC levels while pulldown methods (e.g. MeDIP-seq and hMeDIP-seq) give relative 5mC or 5hmC levels. We have developed a simple classifier that integrates data from these two classes of experiments to label regions 5mC or 5hmC.

### Experimental Setups and Analysis Steps
Our pipeline allows for the integration of bisulfite-converted and pulldown data (hybrid setup) as well as the integration of purely pulldown data (pulldown setup). Integration can occur at the sample/replicate level, or at the differential level where two conditions are tested.

### Hybrid Sample/Replicate-level Analysis
Bisulfite-converted data from BS-seq or RRBS is analyzed with Bismark (currently v0.14.4) to determine methylation rates per sample/replicate. Pulldown data from hMeDIP-seq or hMe-Seal is aligned to the genome using Bowtie2 and regions of methylation are called using MACS2. Our classifier then combines the methylation calls from Bismark with the peaks from MACS2 to determine sample/replicate-wise regions of 5mC and 5hmC. Quality information and statistics are given along with a track hub for visualization in the UCSC Genome Browser.

### Pulldown Sample/Replicate-level Analysis
Bowtie2 + MACS2 pipeline is used for both the 5mC and 5hmC pulldown data. Quality information, statistics, and a UCSC Genome Browser track hub are provided to the user.

### Hybrid Group Comparison Analysis
Bisulfite-converted data is analyzed with Bismark, and pulldown data with Bowtie2. Regions of differential methylation for bisulfite-converted data are determined using methylSig. Regions of differential methylation for pulldown data are determined using PePr. A classifier integrates the differential methylation information. Summary information and UCSC track hubs are output as before.

### Pulldown Group Comparison
Pulldown data for both 5mC and 5hmC is aligned with Bowtie2 and differential methylation is determined using PePr. Peaks from PePr are then integrated with a classifier and summary information is provided.

## Future Work
As oxBS-seq and TAB-seq become more reliable, we intend to implement a base-resolution classifier.