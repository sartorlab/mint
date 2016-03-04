#!/usr/bin/awk -f
BEGIN {OFS="\t"}
{
	outFile = FILENAME

	if ( ($4 + $5 >= 5) && ($4 / ($4 + $5) >= 0.5) ) {
		### count >= 5 and meth >= 50%
		sub(/_trimmed\.fq\.gz_bismark_bt2\.CpG_report\.txt/, "_highmeth.txt", outFile)
		print $1, $2, $2 > outFile
	} else if ( ($4 + $5 >= 5) && ($4 / ($4 + $5) > 0.05) && ($4 / ($4 + $5) < 0.5) ) {
		### count >= 5 and 5% < meth < 50%
		sub(/_trimmed\.fq\.gz_bismark_bt2\.CpG_report\.txt/, "_lowmeth.txt", outFile)
		print $1, $2, $2 > outFile
	} else if ( ($4 + $5 >= 5) && ($4 / ($4 + $5) <= 0.05) ) {
		### count >= 5 and meth < 5%
		sub(/_trimmed\.fq\.gz_bismark_bt2\.CpG_report\.txt/, "_nometh_signal.txt", outFile)
		print $1, $2, $2 > outFile
	} else {
		### count <= 4
		sub(/_trimmed\.fq\.gz_bismark_bt2\.CpG_report\.txt/, "_nometh_nosignal.txt", outFile)
		print $1, $2, $2 > outFile
	}
}