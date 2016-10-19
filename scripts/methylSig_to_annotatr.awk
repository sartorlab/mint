#!/usr/bin/awk -f
BEGIN {OFS="\t"}
{
	if( $6 < FDR && $7 > DIFF ) {
		print $1, $2, $3, GROUP1, $5, ".", $7, $11, $12
	}
	if( $6 < FDR && $7 < DIFF*(-1) ) {
		print $1, $2, $3, GROUP0, $5, ".", $7, $11, $12
	}
	if( $6 > FDR ) {
		print $1, $2, $3, "noDM", $5, ".", $7, $11, $12
	}
}
