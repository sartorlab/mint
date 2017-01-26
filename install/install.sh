#!/bin/bash
set -e
set -u
set -o pipefail

# Installation directory
INSTALL_DIR=~

# Setup apps folder
mkdir ${INSTALL_DIR}/mint_apps
mkdir ${INSTALL_DIR}/mint_apps/bin
export PATH=${INSTALL_DIR}/mint_apps/bin:$PATH

cd ${INSTALL_DIR}/mint_apps

##############################################
# Download and install samtools
curl -L "https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2" > samtools-1.3.1.tar.bz2
tar -jxvf samtools-1.3.1.tar.bz2
cd samtools-1.3.1
make
make prefix=${INSTALL_DIR}/mint_apps/samtools install
export PATH=${INSTALL_DIR}/mint_apps/samtools/bin:$PATH

cd ${INSTALL_DIR}/mint_apps

# Cleanup
rm samtools-1.3.1.tar.bz2
rm -r samtools-1.3.1

##############################################
# Download and install bedtools
curl -L "https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz" > bedtools-2.26.0.tar.gz
tar -zxvf bedtools-2.26.0.tar.gz
cd bedtools2
make
export PATH=${INSTALL_DIR}/mint_apps/bedtools2/bin:$PATH

cd ${INSTALL_DIR}/mint_apps

# Cleanup
rm bedtools-2.26.0.tar.gz

##############################################
# Download and install bowtie2
curl -L "https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/bowtie2-2.2.9-linux-x86_64.zip" > bowtie2-2.2.9-linux-x86_64.zip
unzip bowtie2-2.2.9-linux-x86_64.zip
export PATH=${INSTALL_DIR}/mint_apps/bowtie2-2.2.9:$PATH

# Cleanup
rm bowtie2-2.2.9-linux-x86_64.zip

##############################################
# Download and install bismark
curl -L "https://github.com/FelixKrueger/Bismark/archive/0.16.1.zip" > Bismark-0.16.1.zip
unzip Bismark-0.16.1.zip
export PATH=${INSTALL_DIR}/mint_apps/Bismark-0.16.1:$PATH

# Cleanup
rm Bismark-0.16.1.zip

##############################################
# Download and install Genome Browser Utilities
curl -L "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig" > bin/bedGraphToBigWig
curl -L "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed" > bin/bedToBigBed
chmod 755 bin/bedGraphToBigWig
chmod 755 bin/bedToBigBed

##############################################
# Download and install bedops
curl -L "https://github.com/bedops/bedops/releases/download/v2.4.14/bedops_linux_x86_64-v2.4.14.tar.bz2" > bedops_linux_x86_64-v2.4.14.tar.bz2
tar -jxvf bedops_linux_x86_64-v2.4.14.tar.bz2

# Cleanup
rm bedops_linux_x86_64-v2.4.14.tar.bz2

##############################################
# Download and install FastQC
curl -L "http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip" > fastqc_v0.11.5.zip
unzip fastqc_v0.11.5.zip
chmod 755 FastQC/fastqc
export PATH=${INSTALL_DIR}/mint_apps/FastQC:$PATH
rm fastqc_v0.11.5.zip

##############################################
# Download and install trim_galore
curl -L "http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.1.zip" > trim_galore_v0.4.1.zip
unzip trim_galore_v0.4.1.zip
mv trim_galore_zip trim_galore
export PATH=${INSTALL_DIR}/mint_apps/trim_galore:$PATH

# Cleanup
rm trim_galore_v0.4.1.zip

##############################################
# Install pip
# NOTE: YOU MAY ALREADY HAVE PIP INSTALLED
curl -L "https://bootstrap.pypa.io/get-pip.py" > get-pip.py
python get-pip.py --user
export PATH=~/.local/bin:$PATH

# Install virtualenv
pip install --user virtualenv

# Setup mint virtual environment
virtualenv mint_env

source mint_env/bin/activate
pip install pysam
pip install multiqc
pip install MACS2==2.1.0.20140616
pip install cutadapt==1.10
pip install PePr

##############################################
# Install R packages
# NOTE: R >= 3.3.0 is required
Rscript r_packages.R
