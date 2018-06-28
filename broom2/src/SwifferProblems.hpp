#ifndef __SWIFFERPROBLEM_HH__
#define __SWIFFERPROBLEM_HH__

/// \file SwifferProblems.hpp
/// \brief Interface for all the problem dependent material properties,
/// solutions, etc.

#include <memory>
#include <string>
#include <vector>
#include "mfem.hpp"

//-----------------------------------------------------------------------

/// \brief Generic interface for Swiffer problem specifications.
///
/// This pure virtual class provides an interface
/// that lets us select different problems at run time.  While it
/// add a layer of indirection, a layer like this exists in
/// production codes so that you can run different problems.
///
/// For the square and cube meshes, the boundary condition
/// ID numbers indirection to:
///    - 1 = minimum x
///    - 2 = maximum x
///    - 3 = minimum y
///    - 4 = maximum y
///    - 5 = minimum z
///    - 6 = maximum z
///
class SwifferProblemSpec
{
  public:
   virtual ~SwifferProblemSpec() = 0;

   /// The external source in the radiation equation.
   virtual double computeS0(mfem::Vector& x) = 0;
   /// The total cross section
   virtual double computeSigmaT(mfem::Vector &x) = 0;
   /// The absorption cross section
   virtual double computeSigmaA(mfem::Vector &x) = 0;
   /// The scattering cross section
   virtual double computeSigmaS(mfem::Vector &x) = 0;
   /// The diffusion coefficient
   virtual double computeD(mfem::Vector &x) = 0;
   /// The diffusion equation isotropic source
   virtual double computeDiffQ(mfem::Vector &x) = 0;
   /// The inflow boundary condition
   virtual double computeInflow(mfem::Vector &x) = 0;
   /// The analytic solution
   virtual double computeAnalyticSolution(mfem::Vector &x) = 0;
   /// The initial condition
   virtual double computeInitialCond(mfem::Vector &x) = 0;

   /// For some problems, we have an analytic answer, and so we should compute an
   /// error metric.
   virtual bool haveAnalyticSolution() const = 0;

   /// Name for printing, and other informational stuff.
   virtual const std::string getName() const = 0;

   /// We can run in three different coordinate systems.  By default, we are
   /// Cartesian, and either XY or XYZ.  We can also run in cylindrical RZ
   /// coordinates for some problems.
   virtual bool isCylindrical() const { return false; }
  protected:
   std::vector<int> mDirichletIds;
   std::vector<int> mReflectingIds;
   std::vector<int> mIncomingFluxIds;
};

//-----------------------------------------------------------------------

/// Set up the various problems, and return the problem parameters via the
/// arguments.
void setProblemSpecification(const std::string& problemName,
                             const double epsThick,
                             std::unique_ptr<SwifferProblemSpec>& problemSpec);

//-----------------------------------------------------------------------

#endif  // __SWIFFERPROBLEM_HH__
