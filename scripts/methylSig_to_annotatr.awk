#!/usr/bin/awk -f
BEGIN {OFS="\t"}
{
	if( $6 < FDR && $7 > DIFF ) {
		print $1, $2, $3, $5, GROUP1, $7, $11, $12
	}
	if( $6 < FDR && $7 < DIFF*(-1) ) {
		print $1, $2, $3, $5, GROUP0, $7, $11, $12
	}
	if( $6 > FDR ) {
		print $1, $2, $3, $5, "noDM", $7, $11, $12
	}
}
