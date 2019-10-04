// triangles_vs.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <mkl_spblas.h>
#include <malloc.h>
#include <time.h>


MKL_INT m = 5, nnz = 13;



int main()
{
	int NNZ = 6629222;
	MKL_INT sizeOfMatrix = 448695;
	MKL_INT* row = (MKL_INT*)malloc(NNZ * sizeof(MKL_INT));
	MKL_INT* col = (MKL_INT*)malloc(NNZ * sizeof(MKL_INT));
	float* val = (float*)malloc(NNZ * sizeof(float));
	sparse_matrix_t A_COO, A, A2;
	sparse_status_t status;
	clock_t start, end;

	FILE* fp;
	//char buff[255];
	int buff_int=0;
	fp = fopen("auto_int.txt", "r");
	for (int i = 0; i < NNZ; i++)
	{
		fscanf(fp, "%d", &row[i]);
		fscanf(fp, "%d", &col[i]);
		fscanf(fp, "%f", &val[i]);
	}
	
	start = clock();
	
	// Creating the Sparse Matrix A in COO format
	status = mkl_sparse_s_create_coo(&A_COO, SPARSE_INDEX_BASE_ONE, sizeOfMatrix, sizeOfMatrix, NNZ, row, col, val);
	if (status == SPARSE_STATUS_SUCCESS)
		printf("Matrix A created with SUCCESS.\n");

	// Convert Sparse Matrix A to CSR format
	status = mkl_sparse_convert_csr(A_COO, SPARSE_OPERATION_NON_TRANSPOSE, &A);
	if (status == SPARSE_STATUS_SUCCESS)
		printf("Matrix A converted to CSR with SUCCESS.\n\n");


	struct matrix_descr generalDesc;
	generalDesc.type = SPARSE_MATRIX_TYPE_GENERAL;
	status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, A, A,&A2);
	if (status == SPARSE_STATUS_SUCCESS)
		printf("A^2 multiplication done with SUCCESS.\n");
	
	
	end = clock();
	// Calculating total time taken by the program. 
	double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
	printf("Multiplication time = %f\n\n", time_taken);

	MKL_INT no_rows_A2, no_cols_A2, * rows_start_A2, * rows_end_A2, * col_index_A2;
	MKL_INT no_rows_A, no_cols_A, * rows_start_A, * rows_end_A, * col_index_A;
	float *values_A2, *values_A;
	sparse_index_base_t sparse_index_A2, sparse_index_A;

	status = mkl_sparse_s_export_csr(A2, &sparse_index_A2, &no_rows_A2, &no_cols_A2, &rows_start_A2, &rows_end_A2, &col_index_A2, &values_A2);
	if (status == SPARSE_STATUS_SUCCESS)
		printf("A^2 export done with SUCCESS.\n");
	//printf("%d\t%d\t%d\t%p\t%p\t%p\n", sparse_index_A2, no_rows_A2, no_cols_A2, &rows_start_A2[1], rows_end_A2, col_index_A2);

	status = mkl_sparse_s_export_csr(A, &sparse_index_A, &no_rows_A, &no_cols_A, &rows_start_A, &rows_end_A, &col_index_A, &values_A);
	if (status == SPARSE_STATUS_SUCCESS)
		printf("A export done with SUCCESS.\n");


	printf("\nNNZ of A^2 = %d\n", rows_start_A2[no_rows_A2]-sparse_index_A2);

	printf("Values: \t");
	for (int i = 0; i <10; i++)
	{
		printf("%.0f\t", values_A2[i]);
	}
	printf("\n");

	printf("Columns: \t");
	for (int i = 0; i < 10; i++)
	{
		printf("%d\t", col_index_A2[i]);
	}
	printf("\n");

	printf("Rows_start: \t");
	for (int i = 0; i < 10; i++)
	{
		printf("%d\t", rows_start_A2[i]);
	}
	printf("\n");

	printf("Rows_end: \t");
	for (int i = 0; i < 10; i++)
	{
		printf("%d\t", rows_end_A2[i]);
	}
	printf("\n");

	float sum = 0;
	for (int r_index = 0; r_index < no_rows_A; r_index++)
	{
		int A_lower_bound = rows_start_A[r_index] - 1;
		int A_upper_bound = rows_end_A[r_index] - 1;
		int A2_lower_bound = rows_start_A2[r_index] - 1;
		int A2_upper_bound = rows_end_A2[r_index] - 1;
		
		int A_c_index = A_lower_bound;
		int A2_c_index = A2_lower_bound;

		while (A_c_index>=A_lower_bound && A_c_index<A_upper_bound && A2_c_index >= A2_lower_bound && A2_c_index < A2_upper_bound) {
			if (col_index_A[A_c_index] == col_index_A2[A2_c_index]) {
				sum += values_A[A_c_index] * values_A2[A2_c_index];
				A_c_index++;
				A2_c_index++;
			}
			else if (col_index_A[A_c_index] < col_index_A2[A2_c_index]) {
				A_c_index++;
			}
			else {
				A2_c_index++;
			}
		}

	}

	printf("\n\nsum = %.2f", sum);

	float nT = sum / 6;
	printf("\nnT = %.2f", nT);


	mkl_sparse_destroy(A_COO);
	mkl_sparse_destroy(A);
	mkl_sparse_destroy(A2);
	free(row);
	free(col);
	free(val);

}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file