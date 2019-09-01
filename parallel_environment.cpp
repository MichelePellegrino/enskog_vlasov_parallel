#include "parallel_environment.hpp"
#include "display.hpp"

#include <cassert>

ParallelEnvironment::ParallelEnvironment
(DSMC* dsmc, MPI_Comm comm):
Motherbase(dsmc), communicator(comm)
{
  // Check if it's initialized
  int init;
  MPI_Initialized( &init );
  initialized = (bool)init;
  assert(initialized && "MPI environment not initialized");
  // Get rank and size
  MPI_Comm_rank(communicator, &rank);
  MPI_Comm_size(communicator, &size);
  size2 = size*size;
}

ParallelEnvironment::ParallelEnvironment
(DSMC* dsmc):
ParallelEnvironment(dsmc, MPI_COMM_WORLD)
{

}

int
ParallelEnvironment::seed_rank
(int seed, int offset)
{
  return seed + offset*rank;
}

void
ParallelEnvironment::barrier
(void)
{
  MPI_Barrier(communicator);
}

// DEBUG
void
ParallelEnvironment::test_output(void)
{
  if (rank == MPI_MASTER)
    std::cout << "RANK " << rank << " WELCOMES YOU" << std::endl;
}
