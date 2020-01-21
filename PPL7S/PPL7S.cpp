#include <mpi.h>
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
	int ierr;
	ierr = MPI_Init(&argc, &argv);
	if (ierr != MPI_SUCCESS) {
		return ierr;
	}
	cout << "Hello World!\n";
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	cout << "Rank: " << rank << endl;
	MPI_Finalize();
}
