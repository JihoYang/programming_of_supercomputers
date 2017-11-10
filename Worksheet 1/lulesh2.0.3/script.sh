#!/bin/bash

rm output.txt

#Array of flags
array='-flto'

for i in $array
do 
	FLAGCXX="CXXFLAGS = -g -O3 -I. -Wall -march=native -fomit-frame-pointer -floop-strip-mine -funroll-loops -floop-interchange -floop-block $i"
	FLAGLD="LDFLAGS = -g -O3 -march=native -fomit-frame-pointer -floop-strip-mine -funroll-loops -floop-interchange -floop-block $i"

	# Delete lines and insert flags on the Makefile
	sed -i '29d' Makefile
	sed -i '29i '"$FLAGCXX"'' Makefile

	sed -i '30d' Makefile
	sed -i '30i '"$FLAGLD"'' Makefile

	# Loop to make an average
	for N in {1..5}
	do
		make clean
		make 

		echo ========================================================================================================================================= >> output.txt
		echo Current flag - incl -march=native -fomit-frame-pointer -floop-strip-mine -funroll-loops -floop-interchange -floop-block : $i >> output.txt
		echo ========================================================================================================================================= >> output.txt

		./lulesh2.0 >> output.txt

		echo ---------------------------- running lulesh2.0 binary $N times------------------------

	done

done

