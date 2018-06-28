
#include "mfem.hpp"
#include "ProblemSpecifications.hpp"
#include <vector>
#include <memory>

#ifndef SOLVE_HPP
#define SOLVE_HPP

class Solve {

public:
  Solve(int dim, int &qloop, int &q,
        void (*omega_function)(const mfem::Vector &x, mfem::Vector &v),
        double (*inflow_function)(mfem::Vector &x),
        double (*psi0_function)(mfem::Vector &x),
        double (*sigma_t_function)(mfem::Vector &x),
        double (*sigma_s_function)(mfem::Vector &x),
        double (*S0_function)(mfem::Vector &x),
        mfem::FiniteElementSpace &fes,
        mfem::GridFunction &phi_old, mfem::GridFunction &phi, int myRank,
        int numAngles, int numProcs, bool directSolve,
        std::unique_ptr<ProblemSpecification> &problemSpec);

  // mfem::GridFunction &getPhi() { return phi; }
  const mfem::GridFunction &getSigt() const { return sigma_t_vis; }
  const mfem::GridFunction &getSigs() const { return sigma_s_vis; }
  const mfem::GridFunction &getScatSrc() const { return scat_src_vis; }
  const mfem::GridFunction &getQext() const { return qext_vis; }
  const mfem::GridFunction &getPsi0() const { return psi0_vis; }
  const mfem::GridFunction &getInflow() const { return inflow_vis; }
  const mfem::GridFunction &getPsi() const { return psi; }

private:
  // mfem::GridFunction phi;
  mfem::GridFunction psi;
  mfem::GridFunction sigma_t_vis;
  mfem::GridFunction sigma_s_vis;
  mfem::GridFunction scat_src;
  mfem::GridFunction scat_src_vis;
  mfem::GridFunction qext_vis;
  mfem::GridFunction psi0_vis;
  mfem::GridFunction inflow_vis;
};

#endif
