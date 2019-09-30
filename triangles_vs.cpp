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
	sparse_matrix_t A, A_CSR, A2;
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
	status = mkl_sparse_s_create_coo(&A, SPARSE_INDEX_BASE_ONE, sizeOfMatrix, sizeOfMatrix, NNZ, row, col, val);
	if (status == SPARSE_STATUS_SUCCESS)
		printf("Matrix A created with SUCCESS.\n");

	// Convert Sparse Matrix A to CSR format
	status = mkl_sparse_convert_csr(A, SPARSE_OPERATION_NON_TRANSPOSE, &A_CSR);
	if (status == SPARSE_STATUS_SUCCESS)
		printf("Matrix A converted to CSR with SUCCESS.\n");


	struct matrix_descr generalDesc;
	generalDesc.type = SPARSE_MATRIX_TYPE_GENERAL;
	status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, A_CSR, A_CSR,&A2);
	if (status == SPARSE_STATUS_SUCCESS)
		printf("A^2 multiplication done with SUCCESS.\n");
	
	
	end = clock();
	// Calculating total time taken by the program. 
	double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
	printf("Multiplication time = %f\n", time_taken);



	

	mkl_sparse_destroy(A);
	mkl_sparse_destroy(A_CSR);
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
