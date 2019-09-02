#!/bin/bash
infile=ForSTIR_20mln_TRUE_WLS_SiPM_unattenuated.txt
outfile=out.root
mln_number_con=1

scanner_lenght=500
number_parts=100
root_template=root_template.root

./main $infile $root_template $outfile $mln_number_con $scanner_lenght $number_parts


