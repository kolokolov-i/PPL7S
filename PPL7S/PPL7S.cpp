#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>

int nomatrix;
int procCount;
int rank;
int* MatrixChunk;
int mSize;
int* displs;
int* sendcounts;
typedef struct { int v1; int v2; } Edge;
int* MST;
int minWeight;
FILE* f_matrix;
FILE* f_time;
FILE* f_result;


int main(int argc, char* argv[])
{

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procCount);
	if (rank == 0) {
		f_matrix = fopen("Matrix.txt", "r");
		if (f_matrix) {
			fscanf(f_matrix, "%d\n", &mSize);
		}
		if (feof(f_matrix))
			nomatrix = 1;
		else
			nomatrix = 0;
	}
	MPI_Bcast(&mSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int i, j = 0, k;
	displs = (int*)malloc(sizeof(int) * procCount);
	sendcounts = (int*)malloc(sizeof(int) * procCount);
	displs[0] = 0;
	sendcounts[0] = mSize / procCount;
	int remains = procCount - (mSize % procCount);
	for (i = 1; i < remains; ++i) {
		sendcounts[i] = sendcounts[0];
		displs[i] = displs[i - 1] + sendcounts[i - 1];
	}
	for (i = remains; i < procCount; ++i) {
		sendcounts[i] = sendcounts[0] + 1;
		displs[i] = displs[i - 1] + sendcounts[i - 1];
	}
	int* matrix;
	matrix = (int*)malloc(mSize * mSize * sizeof(int));
	if (rank == 0) {
		for (i = 0; i < mSize; ++i) {
			matrix[mSize * i + j] = 0;
			for (j = i + 1; j < mSize; ++j) {
				int buf;
				if (!nomatrix) {
					fscanf(f_matrix, "%d\n", &buf);
					matrix[mSize * i + j] = matrix[mSize * i + j] = buf;
				}
				//printf("matrix(%d,%d) is %d\n", i, j, matrix[mSize * i + j]);
			}
		}
		fclose(f_matrix);
	}
	MatrixChunk = (int*)malloc(sendcounts[rank] * mSize * sizeof(int));
	MPI_Datatype matrixString;
	MPI_Type_contiguous(mSize, MPI_INT, &matrixString);
	MPI_Type_commit(&matrixString);
	MPI_Scatterv(matrix, sendcounts, displs, matrixString, MatrixChunk, sendcounts[rank], matrixString, 0, MPI_COMM_WORLD);
	//if (rank == 0)
	//{
	//	for (i = 0; i < mSize; ++i)
	//	{
	//		for (j = 0; j < mSize; ++j)
	//		{
	//			printf("%d ", MatrixChunk[mSize * i + j]);
	//		}
	//		printf("\n");
	//	}
	//}
	MST = (int*)malloc(sizeof(int) * mSize);
	for (i = 0; i < mSize; ++i) {
		MST[i] = -1;
	}
	double start;
	start = MPI_Wtime();
	MST[0] = 0;
	minWeight = 0;
	int min;
	int v1 = 0;
	int v2 = 0;
	struct { int min; int rank; } minRow, row;
	Edge edge;
	for (k = 0; k < mSize - 1; ++k) {
		min = INT_MAX;
		for (i = 0; i < sendcounts[rank]; ++i) {
			if (MST[i + displs[rank]] != -1) {
				for (j = 0; j < mSize; ++j) {
					if (MST[j] == -1) {
						if (MatrixChunk[mSize * i + j] < min && MatrixChunk[mSize * i + j] != 0) {
							min = MatrixChunk[mSize * i + j];
							v2 = j; // change the current edge
							v1 = i;
						}
					}
				}
			}
		}
		row.min = min;
		row.rank = rank;
		MPI_Allreduce(&row, &minRow, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
		edge.v1 = v1 + displs[rank];
		edge.v2 = v2;
		MPI_Bcast(&edge, 1, MPI_2INT, minRow.rank, MPI_COMM_WORLD);
		MST[edge.v2] = edge.v1;
		minWeight += minRow.min;
	}
	double finish, calc_time;
	finish = MPI_Wtime();
	calc_time = finish - start;
	if (rank == 0) {
		f_result = fopen("result.txt", "w");
		fprintf(f_result, "The min minWeight is %d\n ", minWeight);
		for (i = 0; i < mSize; ++i)
		{
			fprintf(f_result, "Edge %d %d\n", i, MST[i]);
		}
		fclose(f_result);
		f_time = fopen("time.txt", "a+");
		fprintf(f_time, "\n Number of processors: %d\n Number of vertices: %d\n Time of execution: %f\n Total Weight: %d\n\n", procCount, mSize, calc_time, minWeight);
		fclose(f_time);
		printf("\n Number of processors: %d\n Number of vertices: %d\n Time of execution: %f\n Total Weight: %d\n\n", procCount, mSize, calc_time, minWeight);
	}
	MPI_Finalize();
	return 0;

}