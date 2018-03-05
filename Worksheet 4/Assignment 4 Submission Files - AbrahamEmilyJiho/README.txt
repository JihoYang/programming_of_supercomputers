The submission files also include some shell scripts for pre/post-processing and .c code for generating binary files for MPI IO.

Description of the files:

<generateBinaryInput.c> reads the input matrices in .in file format with sequential POSIX IO (exactly the same as the baseline implementation) but exports these with MPI IO as a binary file.
Since this has to be done only once, it was used as a pre-processing step, and executed separately to the cannon_io.c

<preProcessBinary.sh> is the job script used to submit generateBinaryInput.c

<postProcessBinary.sh> is a shell script to convert the resulting binary output matrices (C) from cannon_io.c to ASCII.txt formats

<cannon_matrices_binary> contains the input matrices in binary format generated from generateBinaryInput.c
