#!/usr/bin/awk -f
BEGIN {OFS="\t"}
{
	if( $6 < FDR && $7 > DIFF && NR > 1 ) {
		print $1, $2 - 1, $3, GROUP1, $5, ".", $7, $11, $12
	}
	if( $6 < FDR && $7 < DIFF*(-1) && NR > 1 ) {
		print $1, $2 - 1, $3, GROUP0, $5, ".", $7, $11, $12
	}
	if( $6 > FDR && NR > 1 ) {
		print $1, $2 - 1, $3, "noDM", $5, ".", $7, $11, $12
	}
}
