#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

// Note we are taking total time as the whole execution time except MPI initialisation (so from start of the MPI Input to before freeing memory)

int main (int argc, char **argv) {
	FILE *fp;
	double **A = NULL, **B = NULL, **C = NULL;
	double *A_local_block = NULL, *B_local_block = NULL, *C_local_block = NULL;
	int A_rows, A_columns, A_local_block_rows, A_local_block_columns, A_local_block_size;
	int B_rows, B_columns, B_local_block_rows, B_local_block_columns, B_local_block_size;
	int rank, size, sqrt_size, matrices_a_b_dimensions[4];
	MPI_Comm cartesian_grid_communicator, row_communicator, column_communicator;
	MPI_Status status; 
	int i, j, k;
	int row, column;
	double compute_time = 0, mpi_time = 0, total_time = 0, input_time = 0, output_time = 0, compute_start, mpi_start,  input_start, output_start;

	// used to manage the cartesian grid
	int dimensions[2], periods[2], coordinates[2], remain_dims[2];

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	/* For square mesh */
	sqrt_size = (int)sqrt((double) size);             
	if(sqrt_size * sqrt_size != size){
		if( rank == 0 ) perror("need to run mpiexec with a perfect square number of processes\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	// create a 2D cartesian grid 
	dimensions[0] = dimensions[1] = sqrt_size;
	periods[0] = periods[1] = 1;    
	MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, periods, 1, &cartesian_grid_communicator);
	MPI_Cart_coords(cartesian_grid_communicator, rank, 2, coordinates);

	// create a row communicator
	remain_dims[0] = 0;            
	remain_dims[1] = 1; 
	MPI_Cart_sub(cartesian_grid_communicator, remain_dims, &row_communicator);

	// create a column communicator
	remain_dims[0] = 1;
	remain_dims[1] = 0;
	MPI_Cart_sub(cartesian_grid_communicator, remain_dims, &column_communicator);

	//------------------------------------------- Start of MPI Input -------------------------------------------//

	// Start measuring MPI Input time - this is also used for total time
    input_start = MPI_Wtime();

	// Define MPI file handlers
	MPI_File fh_A, fh_B;
	// Define file names
	char filenameA[300], filenameB[300];
	// Error output
	int ierr;
	// Set the file names
	sprintf(filenameA, "%s\0", argv[1]);
	sprintf(filenameB, "%s\0", argv[2]);

	// Open matrix A file (filenameA)
	ierr = MPI_File_open(MPI_COMM_WORLD, filenameA, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_A);
	// Check if the MPI open was successful
	if (ierr!=MPI_SUCCESS){
		printf("MPI_FILE_OPEN failed for matrix A\n");
	}

	// Read the A matrix dimensions 
	ierr = MPI_File_read(fh_A, &matrices_a_b_dimensions[0], 1, MPI_INT, &status);
	if (ierr!=MPI_SUCCESS){
		printf("MPI_FILE_READ failed for A matrix dimensions (row)\n");
	}
	ierr = MPI_File_read(fh_A, &matrices_a_b_dimensions[1], 1, MPI_INT, &status);
	if (ierr!=MPI_SUCCESS){
		printf("MPI_FILE_READ failed for A matrix dimensions (column)\n");
	}

	// Set the displacements to ignore the headers for reading the matrices
	MPI_Offset disp = 2*sizeof(int);

	// Update variables related to A matrix dimensions and local metadata for A
	A_rows = matrices_a_b_dimensions[0];
	A_columns = matrices_a_b_dimensions[1];
	A_local_block_rows = A_rows / sqrt_size;
	A_local_block_columns = A_columns / sqrt_size;
	A_local_block_size = A_local_block_rows * A_local_block_columns;
	A_local_block = (double *) malloc (A_local_block_size * sizeof(double));

	// Parameters for sub-array of A. 
	int array_of_sizesA[2];
	int array_of_subsizesA[2];
	int array_of_startsA[2];
	array_of_sizesA[0] = A_rows;
	array_of_sizesA[1] = A_columns;
	array_of_subsizesA[0] = A_local_block_rows;
	array_of_subsizesA[1] = A_local_block_columns;
	array_of_startsA[0] = coordinates[0] * array_of_subsizesA[0];	// Gives the x-coordinate of the matrix where the sub-array begins. Boundaries are defined by sub-array sizes
	array_of_startsA[1] = coordinates[1] * array_of_subsizesA[1]; 	// Gives the y-coordinate of the matrix where the sub-array begins. Boundaries are defined by sub-array sizes

	// Create sub-array of A
	MPI_Datatype sub_array_A;
	ierr = MPI_Type_create_subarray(2, array_of_sizesA, array_of_subsizesA, array_of_startsA, MPI_ORDER_C, MPI_DOUBLE, &sub_array_A);
	MPI_Type_commit(&sub_array_A);

	// Set file view
	ierr = MPI_File_set_view(fh_A, disp, MPI_DOUBLE, sub_array_A, "native", MPI_INFO_NULL);
	if(ierr!=MPI_SUCCESS){
		printf("MPI_FILE_SET_VIEW failed for matrix A\n");
	}
	// Read the sub-array of matrix A
	ierr = MPI_File_read_all(fh_A, A_local_block, A_local_block_size, MPI_DOUBLE, &status);
	if(ierr!=MPI_SUCCESS){
		printf("MPI_FILE_READ_ALL failed for matrix A\n");
	}
	// Close file A
	MPI_File_close(&fh_A);

	// Open matrix B file (filenameB)
	ierr = MPI_File_open(MPI_COMM_WORLD, filenameB, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_B);
	// Check if the MPI open was successful
	if (ierr!=MPI_SUCCESS){
		printf("Cannon failed to open matrix B\n");
	}

	// Read the B matrix dimensions
	MPI_File_read(fh_B, &matrices_a_b_dimensions[2], 1, MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_read(fh_B, &matrices_a_b_dimensions[3], 1, MPI_INT, MPI_STATUS_IGNORE);
	// Update variables related to B matrix dimensions and local metadata for B
	B_rows = matrices_a_b_dimensions[2];
	B_columns = matrices_a_b_dimensions[3];
	B_local_block_rows = B_rows / sqrt_size;
	B_local_block_columns = B_columns / sqrt_size;
	B_local_block_size = B_local_block_rows * B_local_block_columns;
	B_local_block = (double *) malloc (B_local_block_size * sizeof(double));

	// Parameters for sub-array of B. 
	int array_of_sizesB[2];
	int array_of_subsizesB[2];
	int array_of_startsB[2];
	array_of_sizesB[0] = B_rows;
	array_of_sizesB[1] = B_columns;
	array_of_subsizesB[0] = B_local_block_rows;
	array_of_subsizesB[1] = B_local_block_columns;
	array_of_startsB[0] = coordinates[0] * array_of_subsizesB[0];	// Gives the x-coordinate of the matrix where the sub-array begins. Boundaries are defined by sub-array sizes
	array_of_startsB[1] = coordinates[1] * array_of_subsizesB[1]; 	// Gives the y-coordinate of the matrix where the sub-array begins. Boundaries are defined by sub-array sizes
	// Create sub-array of B
	MPI_Datatype sub_array_B;
	MPI_Type_create_subarray(2, array_of_sizesB, array_of_subsizesB, array_of_startsB, MPI_ORDER_C, MPI_DOUBLE, &sub_array_B);
	MPI_Type_commit(&sub_array_B);

	// Set file view
	MPI_File_set_view(fh_B, disp, MPI_DOUBLE, sub_array_B, "native", MPI_INFO_NULL);
	// Read the sub-array of matrix B
	MPI_File_read_all(fh_B, B_local_block, B_local_block_size, MPI_DOUBLE, MPI_STATUS_IGNORE);

	// Close file B
	MPI_File_close(&fh_B);

	// Local meta data for C
	C_local_block = (double *) malloc (A_local_block_rows * B_local_block_columns * sizeof(double));
        // C needs to be initialized at 0 (accumulates partial dot-products)
        for(i=0; i < A_local_block_rows * B_local_block_columns; i++){
                C_local_block[i] = 0;
        }

	// End of MPI Input time
        input_time = MPI_Wtime() - input_start;

	//------------------------------------------- End of MPI Input -------------------------------------------//

	// cannon's algorithm
	int cannon_block_cycle;
	int C_index, A_row, A_column, B_column;

	for(cannon_block_cycle = 0; cannon_block_cycle < sqrt_size; cannon_block_cycle++){
		// compute partial result for this block cycle
		// Start measruing compute time
		compute_start = MPI_Wtime();
		for(C_index = 0, A_row = 0; A_row < A_local_block_rows; A_row++){
			for(B_column = 0; B_column < B_local_block_columns; B_column++, C_index++){
				for(A_column = 0; A_column < A_local_block_columns; A_column++){
					C_local_block[C_index] += A_local_block[A_row * A_local_block_columns + A_column] *
						B_local_block[A_column * B_local_block_columns + B_column];
				}
			}
		}
		// End measuring compute time
		compute_time += MPI_Wtime() - compute_start;
		// Start measuring mpi time
		mpi_start = MPI_Wtime();
		// rotate blocks horizontally
		MPI_Sendrecv_replace(A_local_block, A_local_block_size, MPI_DOUBLE, 
				(coordinates[1] + sqrt_size - 1) % sqrt_size, 0, 
				(coordinates[1] + 1) % sqrt_size, 0, row_communicator, &status);
		// rotate blocks vertically
		MPI_Sendrecv_replace(B_local_block, B_local_block_size, MPI_DOUBLE, 
				(coordinates[0] + sqrt_size - 1) % sqrt_size, 0, 
				(coordinates[0] + 1) % sqrt_size, 0, column_communicator, &status);
		// End measuring mpi time
		mpi_time += MPI_Wtime() - mpi_start;
	}

	//------------------------------------------- Start of MPI Output -------------------------------------------//

	    // Start measuring MPI Output time
	    output_start = MPI_Wtime();
	    
	    //Step 1. Create an mpi file handle for the output file of matrix C:
	    MPI_File fh_c;
	    
	    //Step 2. Create suitable name for output file. The name is created using sprintf:
	    char size_name_c[256];
	    sprintf(size_name_c, "%dx%d\0", A_rows, B_columns);
	    char prefix[3] = "c_\0";
	    char *fn_c = malloc(strlen(prefix)+strlen(size_name_c)+1);
	    strcpy(fn_c, prefix);
	    strcat(fn_c, size_name_c);
	    
	    //step 3. Print out the file name
	    //if (rank == 0){
	    //        printf("Output file C: %s\n", fn_c);
	    //}
	    
	    //Step 4. Create and open the output file for c.
	    //MPI_MODE_CREATE is used so that the file is created if the file does not currently exist (should be the case).
	    //MPI_MODE_WRONLY is used so that the file can only be written to. There is no reason to read from the output file.
	    MPI_File_open(MPI_COMM_WORLD, fn_c, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_c);
	    
	    //Step 5. Write the header of the output file. (eg. 64 64):
	    //Only rank 0 has write access in this step:
	    if (rank == 0) {
		MPI_File_write(fh_c, &A_rows, 1, MPI_INT, MPI_STATUS_IGNORE);
		MPI_File_write(fh_c, &B_columns, 1, MPI_INT, MPI_STATUS_IGNORE);
	    }
	    
	    //Since the header has been output, there needs to be an offset so that the processes will write the matrix AFTER the header (and not overwrite the header).
	    MPI_Offset disp_header = 2*sizeof(int);
	    
	    //Step 6. Create an MPI subarray for the C matrix. the subarray is how the data is written in to the file:
	    // the sizesC array holds the total size of the subarray, accounting for all the blocks together:
	    int sizesC[2];
	    sizesC[0] = A_rows;
	    sizesC[1] = B_columns;
	    
	    //the subsizesC array holds the local block sizes for each process.
	    int subsizesC[2];
	    subsizesC[0] = A_local_block_rows;
	    subsizesC[1] = B_local_block_columns;
	    int C_local_block_size = A_local_block_rows * B_local_block_columns;
	    
	    // the startsC holds the position where each of the local process blocks can be written to.
	    int startsC[2];
	    startsC[0] = coordinates[0] * subsizesC[0];
	    startsC[1] = coordinates[1] * subsizesC[1];
	    
	    // Create and commit the subarray. The subarray is used by MPI_File_set_view
	    MPI_Datatype c_subarray;
	    MPI_Type_create_subarray(2, sizesC, subsizesC, startsC, MPI_ORDER_C, MPI_DOUBLE, &c_subarray);
	    MPI_Type_commit(&c_subarray);
	    
	    //Step 7. Set the file view for the C matrix. This sets where each process can view in the file:
	    MPI_File_set_view(fh_c, disp_header, MPI_DOUBLE, c_subarray, "native", MPI_INFO_NULL);
	    
	    //Step 8. Write the C matrix to the output file:
	    //MPI_File_write_all is a collective version of MPI_File_write. All the ranks can write concurrently. still blocking though.
	    MPI_File_write_all(fh_c, C_local_block, C_local_block_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
	    
	    //Step 9. Close the outputfile.
	    MPI_File_close(&fh_c);

	    // End measuring MPI Output time
	    output_time = MPI_Wtime() - output_start;
	    // End measuring total time
	    total_time = MPI_Wtime() - input_start;

	//------------------------------------------- End of MPI Output -------------------------------------------//

	if (rank == 0) {
                printf("(%d,%d)x(%d,%d)=(%d,%d)\n", A_rows, A_columns, B_rows, B_columns, A_rows, B_columns);
                printf("Computation time: %lf\n", compute_time);
                printf("MPI time:         %lf\n", mpi_time);
                printf("Total time:       %lf\n", total_time);
                printf("MPI Input time:	  %lf\n", input_time);
                printf("MPI Output time:  %lf\n", output_time);
	}
	// Free memory
	free(A_local_block);
	free(B_local_block);
	free(C_local_block);
	free(fn_c);
	// finalize MPI
	MPI_Finalize();
}

