#!/bin/bash

convertBinary2ASCII()
{
	#Shell script in order to convert binary files from cannon_io.c to ASCII (readable). 
	#The commands should be 'sh bin2ascii.sh <bin_file_name> <Num_of_x> <num_of_y>'

	# Input file (at location 1) <bin-file_name>
	INFILE=$1

	# number of columns (location 2) <Num_of_x>
	X=$2

	# number of rows (location 3) <Num_of_y>
	Y=$3

	#
	YTOT=${Y}

	# Filename of the output ASCII file:
	OUTFILE=${INFILE}_ASCII.txt

	#remove file if already exists:
	rm -f ${OUTFILE}

	#Use hexdump to convert from binary to ascii. Convert headers:
	hexdump -e '1 4 "%d x" 1 4 " %d \n"' -n 8 ${INFILE} >> ${OUTFILE}


	#We have already written 8 bytes (due to the headers)
	OFFSET=8

	#for the new line of matrix
	let YTOT*=8

	#Use hexdump to convert from binary to ascii. Convert matrix:
	for i in `seq 1 ${X}`;
	do
	    hexdump -e ''"${Y}"' 8 "%5.1f " "\n"' -n ${YTOT} -s ${OFFSET} ${INFILE} >> ${OUTFILE}
	    let OFFSET+=${YTOT}
	    #echo $OFFSET
	done
}

echo ============================== Converting binary outputs to ASCII formats...
convertBinary2ASCII c_64x64 64 64
convertBinary2ASCII c_128x128 128 128
convertBinary2ASCII c_256x256 256 256 
convertBinary2ASCII c_512x512 512 512
convertBinary2ASCII c_1024x1024 1024 1024
convertBinary2ASCII c_2048x2048 2048 2048
convertBinary2ASCII c_4096x4096 4096 4096
echo Conversion completed!
mv c_* OutputMatrices
