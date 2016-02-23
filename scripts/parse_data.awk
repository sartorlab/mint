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

    sampleFile = sprintf("%s/%s/data/raw_fastqs/%s.fastq.gz", pwd, projectID, sampleID)

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
      readableFile = sprintf("%s/%s/%s/%s_%s_%s_input.fastq.gz", pwd, projectID, platform, humanID, platform, mark)
    } else {
      readableFile = sprintf("%s/%s/%s/%s_%s_%s.fastq.gz", pwd, projectID, platform, humanID, platform, mark)
    }

    # Create the symlinks
    print "ln -s "sampleFile" "readableFile

    # Add to the appropriate file group variable
    if ( pulldown == 1 ) {
      if ( input == 1 ) {
        pulldownFiles = pulldownFiles" "sprintf("%s_%s_%s_input.fastq.gz", humanID, platform, mark)
      } else {
        pulldownFiles = pulldownFiles" "sprintf("%s_%s_%s.fastq.gz", humanID, platform, mark)
      }
    } else if ( bisulfite == 1 ) {
      bisulfiteFiles = bisulfiteFiles" "sprintf("%s_%s_%s.fastq.gz", humanID, platform, mark)
    }

  }
}
END {
  print "bisulfiteFiles: "bisulfiteFiles
  print "pulldownFiles: "pulldownFiles
}
