#ifndef EV_PARALLEL_HPP
#define EV_PARALLEL_HPP

#include <mpi.h>

#include "motherbase.hpp"
#include "types.hpp"

#ifndef RANK_RNG_OFFSET
#define RANK_RNG_OFFSET 1000
#endif

#ifndef MPI_MASTER
#define MPI_MASTER 0
#endif

#ifndef MPI_NO_RANK
#define MPI_NO_RANK -1
#endif

class DSMC;

class ParallelEnvironment : protected Motherbase
{

private:
  
  MPI_Comm communicator;
  bool initialized;
  int rank;
  int size;
  int size2;

public:

  ParallelEnvironment(DSMC*, MPI_Comm);
  ParallelEnvironment(DSMC*);
  ~ParallelEnvironment() = default;

  inline int get_rank() const { return rank; };
  inline int get_size() const { return size; };
  inline int get_size2() const { return size2; }

  inline bool is_root() const { return rank==MPI_MASTER; }

  int seed_rank(int seed, int offset = RANK_RNG_OFFSET);

  template<class data_type>
  void broadcast(data_type& buffer, int count = 1, int from = MPI_MASTER)
  {
    MPI_Bcast(&buffer, count, ev_mpi::MPI_Data_convert<data_type>(), from, communicator);
  }

  template<class data_type>
  void all_gather_inp(data_type& buffer, int count = 1)
  {
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &buffer, count,
      ev_mpi::MPI_Data_convert<data_type>(), communicator);
  }

  template<class data_type>
  void send(data_type& buffer, int rank_receiving, int count = 1, int tag = 0)
  {
    MPI_Send(&buffer, count, ev_mpi::MPI_Data_convert<data_type>(),
      rank_receiving, tag, communicator);
  }

  template<class data_type>
  void receive(data_type& buffer, int rank_sending, int count = 1, int tag = 0, MPI_Status* status = MPI_STATUS_IGNORE)
  {
    MPI_Recv(&buffer, count, ev_mpi::MPI_Data_convert<data_type>(),
      rank_sending, tag, communicator, status);
  }

  void barrier(void);

  // DEBUG
  void test_output(void);

};

template<>
inline void ParallelEnvironment::broadcast<bool>(bool& buffer, int count, int from)
{
  int temp;
  MPI_Bcast(&temp, count, MPI_INT, from, communicator);
  buffer = (bool)temp;
}

#endif /* EV_PARALLEL_HPP */
