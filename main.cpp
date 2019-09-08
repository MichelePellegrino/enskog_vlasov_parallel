#include "dsmc.hpp"

int main()
{

MPI_Init(nullptr, nullptr);

try
{
  // DSMC my_dsmc_instanziation("input_files/uniform02.dat.txt");
  DSMC my_dsmc_instanziation("input_files/drop_par.dat.txt");
}
catch(const char* ex)
{
  std::cout << ex << std::endl;
}

MPI_Finalize();

return 0;

}
