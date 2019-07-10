#!/bin/bash
infile=pmt4.txt
outfile=pmt4.root
mln_number_con=5

scanner_lenght=500
number_parts=100
root_template=root_template.root

./main $infile $root_template $outfile $mln_number_con $scanner_lenght $number_parts




