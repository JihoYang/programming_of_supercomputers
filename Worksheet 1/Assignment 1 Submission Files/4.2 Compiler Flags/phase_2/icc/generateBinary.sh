#!/bin/bash

# Array of flags
array='-march=native -xHost -unroll -ipo'
fullarray="array='$array'"

# Update the array of flags on job.cmd
sed -i '23d' job.cmd
sed -i '23i '"$fullarray"'' job.cmd

# Number of flags
X=4

# Update the number of flags on job.cmd
sed -i '24d' job.cmd
sed -i '24i X='"$X"'' job.cmd

# Create binaries with different flags
for i in $array
do 
	
	FLAGCXX="CXXFLAGS = -O3 -I. -Wall $i"
	FLAGLD="LDFLAGS = -O3 $i"

	# Delete lines and insert flags on the Makefile
	sed -i '29d' Makefile
	sed -i '29i '"$FLAGCXX"'' Makefile

	sed -i '30d' Makefile
	sed -i '30i '"$FLAGLD"'' Makefile
	
	make
	mv lulesh2.0 lulesh_$i
	make clean

done

