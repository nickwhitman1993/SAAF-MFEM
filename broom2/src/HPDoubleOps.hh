#ifndef HPDOUBLEOPS_HH__
#define HPDOUBLEOPS_HH__

#include <mpi.h>

/// Register and unregister functions with MPI so that it knows
/// about HPDouble.
class HPDoubleOps
{
  public:
   HPDoubleOps();
   ~HPDoubleOps();

   MPI_Datatype& getType() { return mType; }
   MPI_Op& getSum() { return mSum; }
  private:
   MPI_Datatype mType;
   MPI_Op mSum;
};

#endif  // HPDOUBLEOPS_HH__
