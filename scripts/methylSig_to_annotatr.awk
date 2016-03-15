#!/usr/bin/awk -f
BEGIN {OFS="\t"}
{
  if ( NR > 1 ) {
	if( $6 < FDR && $7 > DIFF ) {
		print $1, $2, $3, "hyper", $5, $4, $7, $11, $12
	}
	if( $6 < FDR && $7 < DIFF*(-1) ) {
		print $1, $2, $3, "hypo", $5, $4, $7, $11, $12
	}
	if( $6 > FDR ) {
		print $1, $2, $3, "noDM", $5, $4, $7, $11, $12
	}
  }
}
