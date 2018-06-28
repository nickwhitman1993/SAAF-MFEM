
#include "mfem.hpp"

#ifndef DSASOLVE_HPP
#define DSASOLVE_HPP

class DSASolve {

public:
   DSASolve(int dim,
   double (*sigma_s_function)(mfem::Vector &x),
   double (*sigma_a_function)(mfem::Vector &x),
   double (*diff_function)(mfem::Vector &x),
   mfem::FiniteElementSpace &fes,
   mfem::GridFunction &phi_old,
   mfem::GridFunction &lDsaPhi,
   mfem::GridFunction &gphi,
   bool directSolve,
   int myRank);
   
private:
   mfem::GridFunction error;
   mfem::GridFunction sigma_s_vis;
   mfem::GridFunction sigma_a_vis;
   mfem::GridFunction diff;
   mfem::GridFunction DiffCoeff;
   mfem::GridFunction scatCorr;
   mfem::GridFunction diff_vis;
};

#endif
