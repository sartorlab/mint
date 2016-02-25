PROJECT=$1
GENOME=$2
DATAPATH=$3

echo Making directories...
mkdir projects/${PROJECT}
mkdir projects/${PROJECT}/data
mkdir projects/${PROJECT}/data/raw_fastqs
mkdir projects/${PROJECT}/{bisulfite,pulldown}
mkdir projects/${PROJECT}/bisulfite/{raw_fastqs,raw_fastqcs,trim_fastqs,trim_fastqcs,bismark,methylsig_calls}
mkdir projects/${PROJECT}/pulldown/{raw_fastqs,raw_fastqcs,trim_fastqs,trim_fastqcs,bowtie2_bams,pulldown_coverages,macs_peaks,pepr_peaks}
mkdir projects/${PROJECT}/classifications
mkdir projects/${PROJECT}/classifications/{simple,sample,comparison}
mkdir projects/${PROJECT}/summary
mkdir projects/${PROJECT}/summary/{figures,tables,reports}
mkdir projects/${PROJECT}/${PROJECT}_hub
mkdir projects/${PROJECT}/${PROJECT}_hub/${GENOME}

echo Copying makefile and config.mk, and making empty variables.mk...
cp template_makefile projects/${PROJECT}/makefile
cp template_config.mk projects/${PROJECT}/config.mk

echo Copying data annotation file to project data folder...
cp projects/${PROJECT}_annotation.txt projects/${PROJECT}/data/

echo Populating file lists...
awk -f scripts/parse_data.awk projects/${PROJECT}/data/${PROJECT}_annotation.txt > projects/${PROJECT}/variables.mk

echo Symlinking data...
ln -s "${DATAPATH}" projects/${PROJECT}/data/raw_fastqs/
