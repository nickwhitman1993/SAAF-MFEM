

#include "HPDouble.hh"
#include "HPDoubleOps.hh"
#include <mpi.h>

namespace
{
void hpsum(void *vin, void *vinout, int *len, MPI_Datatype *dptr)
{
   HPDouble *in = (HPDouble *)vin;
   HPDouble *inout = (HPDouble *)vinout;

   for (int i = 0; i < *len; ++i)
   {
      *inout += *in;
      in++;
      inout++;
   }
}
}

HPDoubleOps::HPDoubleOps()
{
   MPI_Type_contiguous(2, MPI_DOUBLE, &mType);
   MPI_Type_commit(&mType);

   int commute = 1;
   MPI_Op_create(hpsum, commute, &mSum);
}

HPDoubleOps::~HPDoubleOps()
{
   MPI_Type_free(&mType);
   MPI_Op_free(&mSum);
}
