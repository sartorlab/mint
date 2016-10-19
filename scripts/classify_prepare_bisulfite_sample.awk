#!/usr/bin/awk -f
BEGIN {OFS="\t"}
{
	outFile = FILENAME

	if ( ($4 + $5 >= MIN_COV) && ($4 / ($4 + $5) >= 0.5) ) {
		### count >= MIN_COV and meth >= 50%
		sub(/_trimmed_bismark_bt2\.CpG_report\.txt/, "_highmeth.txt", outFile)
		print $1, $2 - 1, $2 > outFile
	} else if ( ($4 + $5 >= MIN_COV) && ($4 / ($4 + $5) > 0.05) && ($4 / ($4 + $5) < 0.5) ) {
		### count >= MIN_COV and 5% < meth < 50%
		sub(/_trimmed_bismark_bt2\.CpG_report\.txt/, "_lowmeth.txt", outFile)
		print $1, $2 - 1, $2 > outFile
	} else if ( ($4 + $5 >= MIN_COV) && ($4 / ($4 + $5) <= 0.05) ) {
		### count >= MIN_COV and meth < 5%
		sub(/_trimmed_bismark_bt2\.CpG_report\.txt/, "_nometh_signal.txt", outFile)
		print $1, $2 - 1, $2 > outFile
	} else {
		### count < MIN_COV
		sub(/_trimmed_bismark_bt2\.CpG_report\.txt/, "_nometh_nosignal.txt", outFile)
		print $1, $2 - 1, $2 > outFile
	}
}
