#!/usr/bin/awk -f
BEGIN {OFS="\t"}
{
    if ( $4 + $5 > 0 ) {
        print $1 "." $2, $1, $2, $3, $4 + $5, ($4 / ($4 + $5))*100, ($5 / ($4 + $5))*100
    }
}
