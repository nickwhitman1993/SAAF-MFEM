
#include "mfem.hpp"

#ifndef DIFFSOLVE_HPP
#define DIFFSOLVE_HPP

class diffSolve {

public:
   diffSolve(int dim,
   double (*sigma_a_function)(mfem::Vector &x),
   double (*diff_Q_function)(mfem::Vector &x),
   double (*diff_function)(mfem::Vector &x),
   mfem::FiniteElementSpace &fes,
   mfem::GridFunction &diffPhi,
   bool directSolve,
   int myRank);
   
private:
   mfem::GridFunction sigma_s_vis;
   mfem::GridFunction sigma_t_vis;
   mfem::GridFunction sigma_a_vis;
   mfem::GridFunction diff;
   mfem::GridFunction DiffCoeff;
};

#endif
