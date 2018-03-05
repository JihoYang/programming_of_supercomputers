#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

// Converts provided matrix.in files to binary files
// This is served as a pre-processing tool and is independent of the MPI IO implementation

int main (int argc, char **argv) {
	FILE *fp;
	double *A_local_block = NULL, *B_local_block = NULL;
	int A_rows, A_columns, A_local_block_rows, A_local_block_columns, A_local_block_size;
	int B_rows, B_columns, B_local_block_rows, B_local_block_columns, B_local_block_size;
	int rank, size, sqrt_size, matrices_a_b_dimensions[4], i;
	MPI_Comm cartesian_grid_communicator, row_communicator, column_communicator;
	MPI_Status status;
	double input_time = 0, output_time = 0, start;

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

	// ====================================== Start of file conversion to binary files ======================================
	// We follow the same methodology taken in the provided cannon.c
	// This approach is not scalable, but consider this as a pre-processing step where it needs to be executed only once
	
	// Define arrays for the matrices
	double **A = NULL, **B = NULL, **C = NULL, *A_array = NULL, *B_array = NULL, *C_array = NULL;

	// getting matrices from files at rank 0 only
	// example: mpiexec -n 64 ./cannon matrix1 matrix2 [test]
	if (rank == 0){

		int row, column;
		if ((fp = fopen (argv[1], "r")) != NULL){
			fscanf(fp, "%d %d\n", &matrices_a_b_dimensions[0], &matrices_a_b_dimensions[1]);
			A = (double **) malloc (matrices_a_b_dimensions[0] * sizeof(double *));
			for (row = 0; row < matrices_a_b_dimensions[0]; row++){
				A[row] = (double *) malloc(matrices_a_b_dimensions[1] * sizeof(double));
				for (column = 0; column < matrices_a_b_dimensions[1]; column++)
					fscanf(fp, "%lf", &A[row][column]);
			}
			fclose(fp);
		} else {
			if(rank == 0) fprintf(stderr, "error opening file for matrix A (%s)\n", argv[1]);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
		if((fp = fopen (argv[2], "r")) != NULL){
			fscanf(fp, "%d %d\n", &matrices_a_b_dimensions[2], &matrices_a_b_dimensions[3]);
			B = (double **) malloc (matrices_a_b_dimensions[2] * sizeof(double *));
			for(row = 0; row < matrices_a_b_dimensions[2]; row++){
				B[row] = (double *) malloc(matrices_a_b_dimensions[3] * sizeof(double *));
				for(column = 0; column < matrices_a_b_dimensions[3]; column++)
					fscanf(fp, "%lf", &B[row][column]);
			}
			fclose(fp);
		} else {
			if(rank == 0) fprintf(stderr, "error opening file for matrix B (%s)\n", argv[2]);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		// need to check that the multiplication is possible given dimensions
		// matrices_a_b_dimensions[0] = row size of A
		// matrices_a_b_dimensions[1] = column size of A
		// matrices_a_b_dimensions[2] = row size of B
		// matrices_a_b_dimensions[3] = column size of B
		if(matrices_a_b_dimensions[1] != matrices_a_b_dimensions[2]){
			if(rank == 0) fprintf(stderr, "A's column size (%d) must match B's row size (%d)\n",
					matrices_a_b_dimensions[1], matrices_a_b_dimensions[2]);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		// this implementation is limited to cases where thematrices can be partitioned perfectly
		if( matrices_a_b_dimensions[0] % sqrt_size != 0
				|| matrices_a_b_dimensions[1] % sqrt_size != 0
				|| matrices_a_b_dimensions[2] % sqrt_size != 0
				|| matrices_a_b_dimensions[3] % sqrt_size != 0 ){
			if(rank == 0) fprintf(stderr, "cannot distribute work evenly among processe\n"
					"all dimensions (A: r:%d c:%d; B: r:%d c:%d) need to be divisible by %d\n",
					matrices_a_b_dimensions[0],matrices_a_b_dimensions[1],
					matrices_a_b_dimensions[2],matrices_a_b_dimensions[3], sqrt_size );
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
	}

	// send dimensions to all peers
	if(rank == 0) {
		for(i = 1; i < size; i++){
			MPI_Send(matrices_a_b_dimensions, 4, MPI_INT, i, 0, cartesian_grid_communicator);
		}
	} else {
		MPI_Recv(matrices_a_b_dimensions, 4, MPI_INT, 0, 0, cartesian_grid_communicator, &status);
	}

	A_rows = matrices_a_b_dimensions[0];
	A_columns = matrices_a_b_dimensions[1];
	B_rows = matrices_a_b_dimensions[2];
	B_columns = matrices_a_b_dimensions[3];

	// local metadata for A
	A_local_block_rows = A_rows / sqrt_size;
	A_local_block_columns = A_columns / sqrt_size;
	A_local_block_size = A_local_block_rows * A_local_block_columns;
	A_local_block = (double *) malloc (A_local_block_size * sizeof(double));

	// local metadata for B
	B_local_block_rows = B_rows / sqrt_size;
	B_local_block_columns = B_columns / sqrt_size;
	B_local_block_size = B_local_block_rows * B_local_block_columns;
	B_local_block = (double *) malloc (B_local_block_size * sizeof(double));

	// full arrays only needed at root
	if(rank == 0){
		A_array = (double *) malloc(sizeof(double) * A_rows * A_columns);
		B_array = (double *) malloc(sizeof(double) * B_rows * B_columns);
		C_array = (double *) malloc(sizeof(double) * A_rows * B_columns);
		// generate the 1D arrays of the matrices at root
		int row, column, i, j;
		for (i = 0; i < sqrt_size; i++){
			for (j = 0; j < sqrt_size; j++){
				for (row = 0; row < A_local_block_rows; row++){
					for (column = 0; column < A_local_block_columns; column++){
						A_array[((i * sqrt_size + j) * A_local_block_size) + (row * A_local_block_columns) + column]
							= A[i * A_local_block_rows + row][j * A_local_block_columns + column];
					}
				}
				for (row = 0; row < B_local_block_rows; row++){
					for (column = 0; column < B_local_block_columns; column++){
						B_array[((i * sqrt_size + j) * B_local_block_size) + (row * B_local_block_columns) + column]
							= B[i * B_local_block_rows + row][j * B_local_block_columns + column];
					}
				}
			}
		}
	}

	// send a block to each process
	if(rank == 0) {
		for(i = 1; i < size; i++){
			MPI_Send((A_array + (i * A_local_block_size)), A_local_block_size, MPI_DOUBLE, i, 0, cartesian_grid_communicator);
			MPI_Send((B_array + (i * B_local_block_size)), B_local_block_size, MPI_DOUBLE, i, 0, cartesian_grid_communicator);
		}
		for(i = 0; i < A_local_block_size; i++){
			A_local_block[i] = A_array[i];
		}
		for(i = 0; i < B_local_block_size; i++){
			B_local_block[i] = B_array[i];
		}
	} else {
		MPI_Recv(A_local_block, A_local_block_size, MPI_DOUBLE, 0, 0, cartesian_grid_communicator, &status);
		MPI_Recv(B_local_block, B_local_block_size, MPI_DOUBLE, 0, 0, cartesian_grid_communicator, &status);
	}

	// Set the filenames for the binary files
	MPI_File fh_A, fh_B;
	// Define char for matrix sizes
	char size_name_a[256];
	sprintf(size_name_a, "%ux%u\0", A_rows, A_columns);	
	char prefix_a[3] = "a_\0";
	char *fn_a = malloc(strlen(prefix_a)+strlen(size_name_a)+1);
	strcpy(fn_a, prefix_a);
	strcat(fn_a, size_name_a);

	char size_name_b[256];
	sprintf(size_name_b, "%ux%u\0", B_rows, B_columns);	
	char prefix_b[3] = "b_\0";
	char *fn_b = malloc(strlen(prefix_b)+strlen(size_name_b)+1);
	strcpy(fn_b, prefix_b);
	strcat(fn_b, size_name_b);

	// Open matrix A file
	MPI_File_open(MPI_COMM_WORLD, fn_a, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_A);	
	// Write the headers of matrix A
	if (rank == 0){
		MPI_File_write(fh_A, &A_rows, 1, MPI_INT, MPI_STATUS_IGNORE);
		MPI_File_write(fh_A, &A_columns, 1, MPI_INT, MPI_STATUS_IGNORE);
	}
	// Set displacements caused by headers
	MPI_Offset disp_header = 2*sizeof(int);

	// Set the parameters for the 2D subarray for A
	int sizesA[2];
	sizesA[0] = A_rows;
	sizesA[1] = A_columns;
	// Set sizes for creating sub-arrays
	int subsizesA[2];
	subsizesA[0] = A_local_block_rows;
	subsizesA[1] = A_local_block_columns;
	int startsA[2];
	startsA[0] = coordinates[0] * subsizesA[0];
	startsA[1] = coordinates[1] * subsizesA[1];
	// Create and commit the sub-array for A
	MPI_Datatype a_subarray;
	MPI_Type_create_subarray(2, sizesA, subsizesA, startsA, MPI_ORDER_C, MPI_DOUBLE, &a_subarray);
	MPI_Type_commit(&a_subarray);
	// Set file view for A
	MPI_File_set_view(fh_A, disp_header, MPI_DOUBLE, a_subarray, "native", MPI_INFO_NULL);
	// Write the A matrix to the output file
	MPI_File_write_all(fh_A, A_local_block, A_local_block_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
	// Close the file for matrix A
	MPI_File_close(&fh_A);

	// Open matrix B file
	MPI_File_open(MPI_COMM_WORLD, fn_b, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_B);	
	// Write the headers of matrix B
	if (rank == 0){
		MPI_File_write(fh_B, &B_rows, 1, MPI_INT, MPI_STATUS_IGNORE);
		MPI_File_write(fh_B, &B_columns, 1, MPI_INT, MPI_STATUS_IGNORE);
	}
	// Set the parameters for the 2D subarray for B
	int sizesB[2];
	sizesB[0] = B_rows;
	sizesB[1] = B_columns;
	// Set sizes for creating sub-arrays
	int subsizesB[2];
	subsizesB[0] = B_local_block_rows;
	subsizesB[1] = B_local_block_columns;
	int startsB[2];
	startsB[0] = coordinates[0] * subsizesB[0];
	startsB[1] = coordinates[1] * subsizesB[1];
	// Create and commit the sub-array for B
	MPI_Datatype b_subarray;
	MPI_Type_create_subarray(2, sizesB, subsizesB, startsB, MPI_ORDER_C, MPI_DOUBLE, &b_subarray);
	MPI_Type_commit(&b_subarray);
	// Set file view for B
	MPI_File_set_view(fh_B, disp_header, MPI_DOUBLE, b_subarray, "native", MPI_INFO_NULL);
	// Write the B matrix to the output file
	MPI_File_write_all(fh_B, B_local_block, B_local_block_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
	// Close the file for matrix B
	MPI_File_close(&fh_B);

	if (rank == 0){
		printf("Generated binary files: %s and %s\n", fn_a, fn_b);
		printf("======================= Input matrix files converted to binary files for A (%dx%d) and B (%dx%d)\n", A_rows, A_columns, B_rows, B_columns);
		printf("\n");
	}
	// ====================================== End of file conversion to binary files ======================================

	// free all memory
	if(rank == 0){
		for(i = 0; i < A_rows; i++){
			free(A[i]);
		}
		for(i = 0; i < B_rows; i++){
			free(B[i]);
		}
		free(A);
		free(B);
		free(A_array);
		free(B_array);
	}

	free(A_local_block);
	free(B_local_block);
	free(fn_a);
	free(fn_b);
	// finalize MPI
	MPI_Finalize();
}
