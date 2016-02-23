#!/usr/bin/awk -f
# The pwd should be the mint root folder
BEGIN {
  OFS="\t"
  pwd=ENVIRON["PWD"]
  pulldownFiles = ""
  bisulfiteFiles = ""
}
{
  if( NR > 1 ) {
    # Establish row variables
    projectID=$1
    sampleID=$2
    humanID=$3
    pulldown=$4
    bisulfite=$5
    mc=$6
    hmc=$7
    input=$8
    group=$9

    sampleFile = sprintf("%s/data/raw_fastqs/%s.fastq.gz", pwd, sampleID)

    # Establish platform for row
    if( pulldown == 1 ) {
      platform = "pulldown"
    } else if ( bisulfite == 1 ) {
      platform = "bisulfite"
    }

    # Establish mark for row
    if( mc == 1 && hmc == 1 ) {
      mark = "mc_hmc"
    } else if ( mc == 1 && hmc == 0 ) {
      mark = "mc"
    } else if ( mc == 0 && hmc == 1 ) {
      mark = "hmc"
    }

    # Make symlink on the basis of input
    if ( input == 1 ) {
      readableFile = sprintf("%s/%s/raw_fastqs/%s_%s_input_%s.fastq.gz", pwd, platform, humanID, mark, platform)
    } else {
      readableFile = sprintf("%s/%s/raw_fastqs/%s_%s_%s.fastq.gz", pwd, platform, humanID, mark, platform)
    }

    # Create the symlinks
    system("ln -s "sampleFile" "readableFile)

  }
}
