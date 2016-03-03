#!/usr/bin/awk -f
# The pwd should be the mint root folder
BEGIN {
	OFS="\t"
	pwd=ENVIRON["PWD"]
	bisAlignFiles = ""
	pullAlignFiles = ""
	pullSampleFiles = ""
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

    # Add to the appropriate file group variable
    if ( pulldown == 1 ) {
      if ( input == 1 ) {
        pullAlignFiles = pullAlignFiles" "sprintf("%s_%s_input_%s.fastq.gz", humanID, mark, platform)
      } else {
        pullAlignFiles = pullAlignFiles" "sprintf("%s_%s_%s.fastq.gz", humanID, mark, platform)
		pullSampleFiles = pullSampleFiles" "sprintf("%s_%s_%s.fastq.gz", humanID, mark, platform)
      }
    } else if ( bisulfite == 1 ) {
      bisAlignFiles = bisAlignFiles" "sprintf("%s_%s_%s.fastq.gz", humanID, mark, platform)
    }

  }
}
END {
	print "BISULFITE_FASTQ_FILES :="bisAlignFiles
	print "PULLDOWN_FASTQ_FILES :="pullAlignFiles
	print "PULLDOWN_SAMPLE_FILES :="pullSampleFiles
}
