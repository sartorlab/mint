#!/usr/bin/awk -f
BEGIN {OFS="\t"}
{
	if ( $4 + $5 > 0 ) {
		print $1, $2, $2, $1 "." $2, $4 + $5, $3, ($4 / ($4 + $5))*100
	}
}
