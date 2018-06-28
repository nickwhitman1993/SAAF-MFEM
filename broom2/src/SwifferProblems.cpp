
/// \file SwifferProblemSpec.cc
/// \brief All the problem dependent material properties, solutions, etc.
///        The problem names are:
///            BoundaryCondition
///            uniformInfMed

#include "SwifferProblems.hpp"
#include "physics.hpp"
#include <cassert>
#include <cmath>

///////////////////////////////////////////////////////////////////////////

/// The required destructor for the pure virtual class
SwifferProblemSpec::~SwifferProblemSpec() {}
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------

/// \brief "Adams" problem to determine the diffusion limit.
///
/// Optical thickness is determined by choice of probEpsilon at runtime.
///      sigma_t = 1 / probEpsilon
///
class BoundaryCondition : public SwifferProblemSpec
{
  public:
   /// \brief Set the parameters of the problem.
   ///
   /// \param sigmaA is the absorption opacity \f$ \sigma_a \f$
   /// \param sigmaT is the total (absorption + scattering) opacity \f$ \sigma_t
   /// \f$
   /// \param epsThickness is the optical thickness scaling parameter
   BoundaryCondition(const double epsThick);
   ~BoundaryCondition() {}
   
   double computeS0(mfem::Vector& x) override;
   double computeSigmaT(mfem::Vector &x) override;
   double computeSigmaA(mfem::Vector &x) override;
   double computeSigmaS(mfem::Vector &x) override;
   double computeD(mfem::Vector &x) override;
   double computeDiffQ(mfem::Vector &x) override;
   double computeInflow(mfem::Vector &x) override;
   double computeAnalyticSolution(mfem::Vector &x) override;
   double computeInitialCond(mfem::Vector &x) override;
   const std::string getName() const override { return "BoundaryCondition"; }
   // Actually, we only have the analytic solution in 1-D
   bool haveAnalyticSolution() const override { return true; }

  private:
   const double mEpsThick;
};

/// \brief A uniform infinite medium test problem.
/// 
///  This problem is a homogeneous uniform infinite medium problem. The BC's
///  are fixed incident angular fluxes.
///
class uniformInfMed : public SwifferProblemSpec
{
  public:
   uniformInfMed(const double epsThick);
   ~uniformInfMed() {}
   
   double computeS0(mfem::Vector& x) override;
   double computeSigmaT(mfem::Vector &x) override;
   double computeSigmaA(mfem::Vector &x) override;
   double computeSigmaS(mfem::Vector &x) override;
   double computeD(mfem::Vector &x) override;
   double computeDiffQ(mfem::Vector &x) override;
   double computeInflow(mfem::Vector &x) override;
   double computeAnalyticSolution(mfem::Vector &x) override;
   double computeInitialCond(mfem::Vector &x) override;
   const std::string getName() const override { return "uniformInfMed"; }
   bool haveAnalyticSolution() const override { return true; }
   
  private:
   const double mEpsThick;
};

///////////////////////////////////////////////////////////////////////////

BoundaryCondition::BoundaryCondition(const double epsThick)
    : SwifferProblemSpec(),
      mEpsThick(epsThick)
{
   // min x BC
   mReflectingIds.push_back(1);

   // mDirichletIds.push_back(2);
   // mDirichletIds.push_back(3);
   mDirichletIds.push_back(4);

   // mIncomingFluxIds.push_back(1);
   mIncomingFluxIds.push_back(2);
   mIncomingFluxIds.push_back(3);
   // mIncomingFluxIds.push_back(4);

   // These only exist in 3D
   mReflectingIds.push_back(5);
   mDirichletIds.push_back(6);
}

//-----------------------------------------------------------------------
double BoundaryCondition::computeS0(mfem::Vector& x)
{
   //     [ S_0 = S_0 / (4*pi)]
   return mEpsThick / (4 * M_PI);
}
//-----------------------------------------------------------------------
double BoundaryCondition::computeSigmaT(mfem::Vector &x)
{
   return 1.0 / mEpsThick;
}
//-----------------------------------------------------------------------
double BoundaryCondition::computeSigmaA(mfem::Vector &x)
{
   return mEpsThick;
}
//-----------------------------------------------------------------------
double BoundaryCondition::computeSigmaS(mfem::Vector &x)
{
   // [ sigma_s = sigma_t - sigma_a ]
   return this->computeSigmaT(x) - this->computeSigmaA(x);
}
//-----------------------------------------------------------------------
double BoundaryCondition::computeD(mfem::Vector &x)
{
   // [ D = 1 / (3 * sigma_t) ]
   return 1 / (3 * computeSigmaT(x));
}
//-----------------------------------------------------------------------
double BoundaryCondition::computeDiffQ(mfem::Vector &x)
{
   // [ Q = S_0 ]
   return mEpsThick;
}
//-----------------------------------------------------------------------
double BoundaryCondition::computeInflow(mfem::Vector &x)
{
   return 0.0;
}
//-----------------------------------------------------------------------
double BoundaryCondition::computeAnalyticSolution(mfem::Vector &x)
{
   // 1-D analytic solution
   double exSol(0);
   double sigma_a(computeSigmaA(x)/1);
   double D(1/(3*computeSigmaT(x)*1));
   double L(sqrt(D / sigma_a));
   double Q(computeDiffQ(x)/1);
  
   double c1(1/4.*Q/sigma_a*((1/4. + 1/2.*D/L)*exp(1/L) - (1/4. - 1/2.*D/L))/ ((1/4. - 1/2.*D/L)*(1/4. - 1/2.*D/L) - (1/4. + 1/2.*D/L)*(1/4. + 1/2.*D/L)*exp(2/L)));
  
   double c2((-c1*(1/4. + 1/2.*D/L)*exp(1/L) - 1/4.*Q/sigma_a)/(1/4. - 1/2.*D/L)*exp(1/L));
  
   exSol = c1*exp(x[0]/L) + c2*exp(-x[0]/L) + Q/sigma_a;
   return exSol;
}
//-----------------------------------------------------------------------
double BoundaryCondition::computeInitialCond(mfem::Vector &x)
{
   return 0.0;
}

//-----------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////

uniformInfMed::uniformInfMed(const double epsThick)
    : SwifferProblemSpec(),
      mEpsThick(epsThick)
{
   //std::cout << "mEpsThick:\t" << mEpsThick << "\n";
   // min x BC
   mReflectingIds.push_back(1);

   // mDirichletIds.push_back(2);
   // mDirichletIds.push_back(3);
   mDirichletIds.push_back(4);

   // mIncomingFluxIds.push_back(1);
   mIncomingFluxIds.push_back(2);
   mIncomingFluxIds.push_back(3);
   // mIncomingFluxIds.push_back(4);

   // These only exist in 3D
   mReflectingIds.push_back(5);
   mDirichletIds.push_back(6);
}
//-----------------------------------------------------------------------
double uniformInfMed::computeS0(mfem::Vector& x)
{
   // S0 = Q / (4*pi)
   return this->computeDiffQ(x) / (4*M_PI);
}
//-----------------------------------------------------------------------
double uniformInfMed::computeSigmaT(mfem::Vector &x)
{
   return 1 / mEpsThick;
}
//-----------------------------------------------------------------------
double uniformInfMed::computeSigmaA(mfem::Vector &x)
{
   return mEpsThick;
}
//-----------------------------------------------------------------------
double uniformInfMed::computeSigmaS(mfem::Vector &x)
{
   return this->computeSigmaT(x) - this->computeSigmaA(x);
}
//-----------------------------------------------------------------------
double uniformInfMed::computeD(mfem::Vector &x)
{
   return 1/(this->computeSigmaT(x)*3);
}
//-----------------------------------------------------------------------
double uniformInfMed::computeDiffQ(mfem::Vector &x)
{
   // Q = S0
   return mEpsThick;
}
//-----------------------------------------------------------------------
double uniformInfMed::computeInflow(mfem::Vector &x)
{
   // psi_inc = phi / sigma_t = s / (4 * pi * sigma_a)
   return computeDiffQ(x) / (4*M_PI* computeSigmaA(x));
}
//-----------------------------------------------------------------------
double uniformInfMed::computeAnalyticSolution(mfem::Vector &x)
{
   // phi_exact = S0 / sigma_a
   return computeDiffQ(x) / computeSigmaA(x);
}
//-----------------------------------------------------------------------
double uniformInfMed::computeInitialCond(mfem::Vector &x)
{
   return 0.0;
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void setProblemSpecification(const std::string& problemName,
                             const double epsThick,
                             std::unique_ptr<SwifferProblemSpec>& problemSpec)
{
   if (problemName == "BoundaryCondition")
   {
      std::cout << "The problem name is: " << problemName <<"\n";
      problemSpec.reset(new BoundaryCondition(epsThick));
   }
   else if (problemName == "uniformInfMed")
   {
      std::cout << "The problem name is: " << problemName << "\n";
      problemSpec.reset(new uniformInfMed(epsThick));
   }
   else
   {
      MFEM_ABORT("Unsupported problem specification: " << problemName);
   }
}

///////////////////////////////////////////////////////////////////////////
