#!/bin/bash

# Load gcc-5
module load gcc/5

# Array of flags
array='-march=native -fomit-frame-pointer -floop-block -floop-interchange -floop-strip-mine -funroll-loops -flto'
fullarray="array='$array'"

# Update the array of flags on job.cmd
sed -i '24d' job.cmd
sed -i '24i '"$fullarray"'' job.cmd

# Number of flags
X=7

# Update the number of flags on job.cmd
sed -i '26d' job.cmd
sed -i '26i X='"$X"'' job.cmd

# Create binaries with different flags
for i in $array
do 
	
	FLAGCXX="CXXFLAGS = -g -O3 -I. -Wall $i"
	FLAGLD="LDFLAGS = -g -O3 $i"

	# Delete lines and insert flags on the Makefile
	sed -i '29d' Makefile
	sed -i '29i '"$FLAGCXX"'' Makefile

	sed -i '30d' Makefile
	sed -i '30i '"$FLAGLD"'' Makefile
	
	make
	mv lulesh2.0 lulesh_$i
	make clean

done

