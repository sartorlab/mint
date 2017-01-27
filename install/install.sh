#!/bin/bash
set -e
set -u
set -o pipefail

# Installation directory
# INSTALL_DIR should be /path/to/mint/
INSTALL_DIR=/path/to/mint
APP_DIR=${INSTALL_DIR}/apps
BIN_DIR=${APP_DIR}/bin

# Setup apps folder
mkdir ${APP_DIR}
mkdir ${BIN_DIR}
export PATH=${BIN_DIR}:$PATH

cd ${APP_DIR}

##############################################
# Download and install samtools
curl -L "https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2" > samtools-1.3.1.tar.bz2
tar -jxvf samtools-1.3.1.tar.bz2
cd samtools-1.3.1
make
make prefix=${APP_DIR}/samtools install
export PATH=${APP_DIR}/samtools/bin:$PATH

cd ${APP_DIR}

# Cleanup
rm samtools-1.3.1.tar.bz2
rm -r samtools-1.3.1

##############################################
# Download and install bedtools
curl -L "https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz" > bedtools-2.26.0.tar.gz
tar -zxvf bedtools-2.26.0.tar.gz
cd bedtools2
make
export PATH=${APP_DIR}/bedtools2/bin:$PATH

cd ${APP_DIR}

# Cleanup
rm bedtools-2.26.0.tar.gz

##############################################
# Download and install bowtie2
curl -L "https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/bowtie2-2.2.9-linux-x86_64.zip" > bowtie2-2.2.9-linux-x86_64.zip
unzip bowtie2-2.2.9-linux-x86_64.zip
mv ${APP_DIR}/bowtie2-2.2.9 ${APP_DIR}/bowtie2
export PATH=${APP_DIR}/bowtie2:$PATH

# Cleanup
rm bowtie2-2.2.9-linux-x86_64.zip

##############################################
# Download and install bismark
curl -L "https://github.com/FelixKrueger/Bismark/archive/0.16.1.zip" > Bismark-0.16.1.zip
unzip Bismark-0.16.1.zip
mv ${APP_DIR}/Bismark-0.16.1 ${APP_DIR}/bismark
export PATH=${APP_DIR}/bismark:$PATH

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
export PATH=${APP_DIR}/FastQC:$PATH

# Cleanup
rm fastqc_v0.11.5.zip

##############################################
# Download and install trim_galore
curl -L "http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.1.zip" > trim_galore_v0.4.1.zip
unzip trim_galore_v0.4.1.zip
mv trim_galore_zip trim_galore
export PATH=${APP_DIR}/trim_galore:$PATH

# Cleanup
rm trim_galore_v0.4.1.zip

##############################################
# Install pip
# Check to see if pip is installed
which pip
if [ $? -ne 0 ]; then
    curl -L "https://bootstrap.pypa.io/get-pip.py" > get-pip.py
    python get-pip.py --user
    export PATH=~/.local/bin:$PATH
fi

# Install virtualenv
# NOTE: Will exit quietly if already installed
pip install --user virtualenv

# Setup mint virtual environment
virtualenv mint_env

# NOTE: On macOS this fails within the script, so virtualenv has to be built interactively
set +u
source ${APP_DIR}/mint_env/bin/activate
set -u
pip install --requirement ../install/pip_packages.txt

deactivate

##############################################
# Install R packages
# NOTE: R >= 3.3.0 is required
Rscript ../install/r_packages.R
