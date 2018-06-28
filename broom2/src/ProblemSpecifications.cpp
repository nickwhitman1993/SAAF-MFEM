
/// \file ProblemSpecification.cc
/// \brief All the problem dependent material properties, solutions, etc.
///        The problem names are:
///            AdamsProblem
///            AlternatingIncident (previously called TP2)
///            MMSGleicher
///            MMSLinear
///            MMSWarsa
///            MMSWoods
///            MultiMaterial (previously called TP3)
///            uniformInfMed
///            WangRagusa

#include "ProblemSpecifications.hpp"
#include "physics.hpp"
#include <cassert>
#include <cmath>
//#include "MfemBlackBody.hpp"
//#include "PhysicalConstants.hpp"

///////////////////////////////////////////////////////////////////////////

/// The required destructor for the pure virtual class
ProblemSpecification::~ProblemSpecification() {}
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
class AdamsProblem : public ProblemSpecification
{
  public:
   /// \brief Set the parameters of the problem.
   ///
   /// \param sigmaA is the absorption opacity \f$ \sigma_a \f$
   /// \param sigmaT is the total (absorption + scattering) opacity \f$ \sigma_t
   /// \f$
   /// \param epsThickness is the optical thickness scaling parameter
   AdamsProblem(std::vector<std::vector<double> > vecDir,
                    const double epsThick);
   ~AdamsProblem() {}
   
   double computeS0(mfem::Vector& x) override;
   double computeSigmaT(mfem::Vector &x) override;
   double computeSigmaA(mfem::Vector &x) override;
   double computeSigmaS(mfem::Vector &x) override;
   double computeD(mfem::Vector &x) override;
   double computeDiffQ(mfem::Vector &x) override;
   double computeInflow(mfem::Vector &x) override;
   double computeAnalyticSolution(mfem::Vector &x) override;
   double computeInitialCond(mfem::Vector &x) override;
   const std::string getName() const override { return "Adams"; }
   // Actually we only have the analytic solution in 1-D
   bool haveAnalyticSolution() const override { return true; }

  private:
   const double mEpsThick;
   std::vector<std::vector<double> > mVecDir;
};

/// \brief A test problem from Palmer and Adams.
/// 
///  This problem is homogeneous, optically thick, highly scattering, zero source,
///      with alternating incident flux on part of two boundaries.
///
class AlternatingIncident : public ProblemSpecification
{
  public:
   AlternatingIncident(std::vector<std::vector<double> > vecDir,
                       const double epsThick);
   ~AlternatingIncident() {}
   
   double computeS0(mfem::Vector& x) override;
   double computeSigmaT(mfem::Vector &x) override;
   double computeSigmaA(mfem::Vector &x) override;
   double computeSigmaS(mfem::Vector &x) override;
   double computeD(mfem::Vector &x) override;
   double computeDiffQ(mfem::Vector &x) override;
   double computeInflow(mfem::Vector &x) override;
   double computeAnalyticSolution(mfem::Vector &x) override;
   double computeInitialCond(mfem::Vector &x) override;
   const std::string getName() const override { return "AlternatingIncident"; }
   bool haveAnalyticSolution() const override { return false; }
   
  private:
   std::vector<std::vector<double> > mVecDir;
   const double mEpsThick;
};

/// \brief A MMS problem from Gleicher et al. (2012).
/// 
/// A SLIGHTLY oscillating manufactured solution to investigate issues
///   that may be found my varying the magnitudes of the transport eq. terms.
///
///  The manufactured solution is
///      psi = (1-mu^2)*(1-eta^2)*sin(1*pi x)*cos(1*pi/2 y)
class MMSGleicher : public ProblemSpecification
{
  public:
   MMSGleicher(std::vector<std::vector<double> > vecDir,
               const double epsThick, double alpha, double beta);
   ~MMSGleicher() {}
   
   double computeS0(mfem::Vector& x) override;
   double computeSigmaT(mfem::Vector &x) override;
   double computeSigmaA(mfem::Vector &x) override;
   double computeSigmaS(mfem::Vector &x) override;
   double computeD(mfem::Vector &x) override;
   double computeDiffQ(mfem::Vector &x) override;
   double computeInflow(mfem::Vector &x) override;
   double computeAnalyticSolution(mfem::Vector &x) override;
   double computeInitialCond(mfem::Vector &x) override;
   const std::string getName() const override { return "MMSGleicher"; }
   bool haveAnalyticSolution() const override { return true; }
   
  private:
   const double mEpsThick;
   //const std::vector<double> mVecDir;
   std::vector<std::vector<double> > mVecDir;
   const double mAlpha;
   const double mBeta;
};

/// \brief A spatially (NOT direction) dependent linear MMS problem.
/// 
/// A SLIGHTLY oscillating manufactured solution to investigate issues
///   that may be found my varying the magnitudes of the transport eq. terms.
///
///  The manufactured solution is
///      psi = ax + by + g
class MMSLinear : public ProblemSpecification
{
  public:
   MMSLinear(std::vector<std::vector<double> > vecDir,
               const double epsThick,
               const double alpha,
               const double beta,
               const double gamma);
   ~MMSLinear() {}
   
   double computeS0(mfem::Vector& x) override;
   double computeSigmaT(mfem::Vector &x) override;
   double computeSigmaA(mfem::Vector &x) override;
   double computeSigmaS(mfem::Vector &x) override;
   double computeD(mfem::Vector &x) override;
   double computeDiffQ(mfem::Vector &x) override;
   double computeInflow(mfem::Vector &x) override;
   double computeAnalyticSolution(mfem::Vector &x) override;
   double computeInitialCond(mfem::Vector &x) override;
   const std::string getName() const override { return "MMSLinear"; }
   bool haveAnalyticSolution() const override { return true; }
   
  private:
   const double mEpsThick;
   //const std::vector<double> mVecDir;
   std::vector<std::vector<double> > mVecDir;
   double mAlpha;
   double mBeta;
   double mGamma;
};

/// \brief An exponential MMS problem by Warsa (2008) and Owens (2016).
/// 
/// A SLIGHTLY oscillating manufactured solution to investigate issues
///   that may be found my varying the magnitudes of the transport eq. terms.
///
///  The manufactured solution is
///      psi = x^2 * y^2 * (alpha - x^2)*(beta - y^2)
class MMSWarsa : public ProblemSpecification
{
  public:
   MMSWarsa(std::vector<std::vector<double> > vecDir,
               const double epsThick,
               const double alpha,
               const double beta);
   ~MMSWarsa() {}
   
   double computeS0(mfem::Vector& x) override;
   double computeSigmaT(mfem::Vector &x) override;
   double computeSigmaA(mfem::Vector &x) override;
   double computeSigmaS(mfem::Vector &x) override;
   double computeD(mfem::Vector &x) override;
   double computeDiffQ(mfem::Vector &x) override;
   double computeInflow(mfem::Vector &x) override;
   double computeAnalyticSolution(mfem::Vector &x) override;
   double computeInitialCond(mfem::Vector &x) override;
   const std::string getName() const override { return "MMSWarsa"; }
   bool haveAnalyticSolution() const override { return true; }
   
  private:
   const double mEpsThick;
   //const std::vector<double> mVecDir;
   std::vector<std::vector<double> > mVecDir;
   double mAlpha;
   double mBeta;
};

/// \brief An exponential MMS problem by Warsa (2008) and Owens (2016).
/// 
/// A SLIGHTLY oscillating manufactured solution to investigate issues
///   that may be found my varying the magnitudes of the transport eq. terms.
///
///  The manufactured solution is
///      psi = x^2 * y^2 * (alpha - x^2)*(beta - y^2)
class MMSWoodsMS : public ProblemSpecification
{
  public:
   MMSWoodsMS(std::vector<std::vector<double> > vecDir,
               const double epsThick,
               const double alpha,
               const double beta,
               const double gamma,
               const double delta);
   ~MMSWoodsMS() {}
   
   double computeS0(mfem::Vector& x) override;
   double computeSigmaT(mfem::Vector &x) override;
   double computeSigmaA(mfem::Vector &x) override;
   double computeSigmaS(mfem::Vector &x) override;
   double computeD(mfem::Vector &x) override;
   double computeDiffQ(mfem::Vector &x) override;
   double computeInflow(mfem::Vector &x) override;
   double computeAnalyticSolution(mfem::Vector &x) override;
   double computeInitialCond(mfem::Vector &x) override;
   const std::string getName() const override { return "MMSWoodsMS"; }
   bool haveAnalyticSolution() const override { return true; }
   
  private:
   const double mEpsThick;
   //const std::vector<double> mVecDir;
   std::vector<std::vector<double> > mVecDir;
   double mAlpha;
   double mBeta;
   double mGamma;
   double mDelta;
};

/// \brief This problem is from Todd Palmer's dissertation. It is very
///        heterogeneous with optically thick and thin regions.
///
class MultiMaterial : public ProblemSpecification
{
  public:
   MultiMaterial(std::vector<std::vector<double> > vecDir);
   ~MultiMaterial() {}
   
   double computeS0(mfem::Vector& x) override;
   double computeSigmaT(mfem::Vector &x) override;
   double computeSigmaA(mfem::Vector &x) override;
   double computeSigmaS(mfem::Vector &x) override;
   double computeD(mfem::Vector &x) override;
   double computeDiffQ(mfem::Vector &x) override;
   double computeInflow(mfem::Vector &x) override;
   double computeAnalyticSolution(mfem::Vector &x) override;
   double computeInitialCond(mfem::Vector &x) override;
   const std::string getName() const override { return "MultiMaterial"; }
   bool haveAnalyticSolution() const override { return false; }
   
  private:
   std::vector<std::vector<double> > mVecDir;
};

/// \brief This problem is from Wang and Ragusa (2010). They primarily varied
///        the scattering ratio. We use this with ../Meshes/WangRagusa2010.mesh
///        or ../Meshes/WangRagusaUnit.mesh for cuved meshes.
///
class WangRagusa : public ProblemSpecification
{
  public:
   WangRagusa(std::vector<std::vector<double> > vecDir,
                  const double epsThick);
   ~WangRagusa() {}
   
   double computeS0(mfem::Vector& x) override;
   double computeSigmaT(mfem::Vector &x) override;
   double computeSigmaA(mfem::Vector &x) override;
   double computeSigmaS(mfem::Vector &x) override;
   double computeD(mfem::Vector &x) override;
   double computeDiffQ(mfem::Vector &x) override;
   double computeInflow(mfem::Vector &x) override;
   double computeAnalyticSolution(mfem::Vector &x) override;
   double computeInitialCond(mfem::Vector &x) override;
   const std::string getName() const override { return "WangRagusa"; }
   bool haveAnalyticSolution() const override { return true; }
   
  private:
   std::vector<std::vector<double> > mVecDir;
   const double mEpsThick;
};

/// \brief A uniform infinite medium test problem.
/// 
///  This problem is a homogeneous uniform infinite medium problem. The BC's
///  are fixed incident angular fluxes.
///
class uniformInfMed : public ProblemSpecification
{
  public:
   uniformInfMed(std::vector<std::vector<double> > vecDir,
                  const double epsThick);
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
   std::vector<std::vector<double> > mVecDir;
   const double mEpsThick;
};

///////////////////////////////////////////////////////////////////////////

AdamsProblem::AdamsProblem(std::vector<std::vector<double> > vecDir,
                           const double epsThick)
    : ProblemSpecification(),
      mEpsThick(epsThick),
      mVecDir(vecDir)
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
double AdamsProblem::computeS0(mfem::Vector& x)
{
   //     [ S_0 = S_0 / (4*pi)]
   return mEpsThick / (4 * M_PI);
}
//-----------------------------------------------------------------------
double AdamsProblem::computeSigmaT(mfem::Vector &x)
{
   return 1.0 / mEpsThick;
}
//-----------------------------------------------------------------------
double AdamsProblem::computeSigmaA(mfem::Vector &x)
{
   return mEpsThick;
}
//-----------------------------------------------------------------------
double AdamsProblem::computeSigmaS(mfem::Vector &x)
{
   // [ sigma_s = sigma_t - sigma_a ]
   return this->computeSigmaT(x) - this->computeSigmaA(x);
}
//-----------------------------------------------------------------------
double AdamsProblem::computeD(mfem::Vector &x)
{
   // [ D = 1 / (3 * sigma_t) ]
   double sigmaT(this->computeSigmaT(x));
   return 1 / (3 * sigmaT);
}
//-----------------------------------------------------------------------
double AdamsProblem::computeDiffQ(mfem::Vector &x)
{
   // [ Q = S_0 ]
   return this->computeS0(x) * (4 * M_PI);
}
//-----------------------------------------------------------------------
double AdamsProblem::computeInflow(mfem::Vector &x)
{
   return 0.0;
}
//-----------------------------------------------------------------------
double AdamsProblem::computeAnalyticSolution(mfem::Vector &x)
{
   // 1-D analytic solution
   double exSol(0);
   double sigmaA(this->computeSigmaA(x));
   double D(this->computeD(x));
   double L(sqrt(D / sigmaA));
   exSol = ((exp(-1/L)-1)/(exp(1/L)-exp(-1/L)))*exp(x[0]/L) - ((exp(-1/L)-1)/(exp(1/L)-exp(-1/L))+1)*exp(-x[0]/L) + 1;
   return exSol;
}
//-----------------------------------------------------------------------
double AdamsProblem::computeInitialCond(mfem::Vector &x)
{
   return 0.0;
}

//-----------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////

AlternatingIncident::AlternatingIncident(std::vector<std::vector<double> > vecDir,
                                          const double epsThick)
    : ProblemSpecification(),
      mVecDir(vecDir),
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
double AlternatingIncident::computeS0(mfem::Vector& x)
{
   return this->computeDiffQ(x) / (4*M_PI);
}
//-----------------------------------------------------------------------
double AlternatingIncident::computeSigmaT(mfem::Vector &x)
{
   return 1000.0;
}
//-----------------------------------------------------------------------
double AlternatingIncident::computeSigmaA(mfem::Vector &x)
{
   return this->computeSigmaT(x) - this->computeSigmaS(x);
}
//-----------------------------------------------------------------------
double AlternatingIncident::computeSigmaS(mfem::Vector &x)
{
   return this->computeSigmaT(x)*0.999;
}
//-----------------------------------------------------------------------
double AlternatingIncident::computeD(mfem::Vector &x)
{
   return 1/(this->computeSigmaT(x)*3);
}
//-----------------------------------------------------------------------
double AlternatingIncident::computeDiffQ(mfem::Vector &x)
{
   return 0;
}
//-----------------------------------------------------------------------
double AlternatingIncident::computeInflow(mfem::Vector &x)
{
  if (x[0] > 0.5) {
    if (x[1] > 0.5) {
      if (x[0] < 0.52) {
        return 1/(1*M_PI);}
      else if (x[0] < 0.54) {
        return 0;}
      else if (x[0] < 0.56) {
        return 1/(1*M_PI);}
      else if (x[0] < 0.58) {
        return 0;}
      else if (x[0] < 0.6) {
        return 1/(1*M_PI);}
      if (x[1] < 0.6) {
        return 1/(1*M_PI);}
      else if (x[1] < 0.7) {
        return 0;}
      else if (x[1] < 0.8) {
        return 1/(1*M_PI);}
      else if (x[1] < 0.9) {
        return 0;}
      else if (x[1] < 1.0) {
        return 1/(1*M_PI);}
      }
    else {
      return 0;}
    }
  else {
    return 0;}
}
//-----------------------------------------------------------------------
double AlternatingIncident::computeAnalyticSolution(mfem::Vector &x)
{
   return 0.0;
}
//-----------------------------------------------------------------------
double AlternatingIncident::computeInitialCond(mfem::Vector &x)
{
   return 0.0;
}

//-----------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////

MMSGleicher::MMSGleicher(std::vector<std::vector<double> > vecDir,
                           const double epsThick,
                           double alpha,
                           double beta)
    : ProblemSpecification(),
      mEpsThick(epsThick),
      mVecDir(vecDir),
      mAlpha(alpha),
      mBeta(beta)
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
double MMSGleicher::computeS0(mfem::Vector& x)
{
   return this->computeDiffQ(x) / (4*M_PI);
}
//-----------------------------------------------------------------------
double MMSGleicher::computeSigmaT(mfem::Vector &x)
{
   return 1/mEpsThick;
}
//-----------------------------------------------------------------------
double MMSGleicher::computeSigmaA(mfem::Vector &x)
{
   //return mEpsThick;
   return this->computeSigmaT(x) - this->computeSigmaS(x);
}
//-----------------------------------------------------------------------
double MMSGleicher::computeSigmaS(mfem::Vector &x)
{
   //return this->computeSigmaT(x) - this->computeSigmaA(x);
   return this->computeSigmaT(x)*0.999;
}
//-----------------------------------------------------------------------
double MMSGleicher::computeD(mfem::Vector &x)
{
   return 1/(this->computeSigmaT(x)*3);
}
//-----------------------------------------------------------------------
double MMSGleicher::computeDiffQ(mfem::Vector &x)
{
   return 4.*M_PI*(mVecDir[0][q] * mAlpha*M_PI * cos(mAlpha*M_PI*x[1]) * (1.-mVecDir[1][q]*mVecDir[1][q]) * (1.-mVecDir[0][q]*mVecDir[0][q]) * cos(mBeta*M_PI*x[0])
      - mVecDir[1][q] * mBeta*M_PI * sin(mBeta*M_PI*x[0]) * (1.-mVecDir[1][q]*mVecDir[1][q]) * (1.-mVecDir[0][q]*mVecDir[0][q]) * sin(mAlpha*M_PI*x[1])
      + this->computeSigmaT(x)*this->computeInflow(x)
      - 1./(4.*M_PI)*this->computeSigmaS(x)*this->computeAnalyticSolution(x));
}
//-----------------------------------------------------------------------
double MMSGleicher::computeInflow(mfem::Vector &x)
{
   return (1-mVecDir[1][q]*mVecDir[1][q])*(1-mVecDir[0][q]*mVecDir[0][q])*sin(mAlpha*M_PI*x[1])*cos(mBeta*M_PI*x[0]);
}
//-----------------------------------------------------------------------
double MMSGleicher::computeAnalyticSolution(mfem::Vector &x)
{
   // psi = (1-mu^2)*(1-eta^2)*sin(2*pi y)*cos(3*pi/2 x)
   return 8*M_PI/5*sin(mAlpha*M_PI*x[1])*cos(mBeta*M_PI*x[0]);
}
//-----------------------------------------------------------------------
double MMSGleicher::computeInitialCond(mfem::Vector &x)
{
   return 0.0;
}

//-----------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////

MMSLinear::MMSLinear(std::vector<std::vector<double> > vecDir,
                           const double epsThick,
                           const double alpha,
                           const double beta,
                           const double gamma)
    : ProblemSpecification(),
      mEpsThick(epsThick),
      mVecDir(vecDir),
      mAlpha(alpha),
      mBeta(beta),
      mGamma(gamma)
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
double MMSLinear::computeS0(mfem::Vector& x)
{
   return this->computeDiffQ(x) / (4*M_PI);
}
//-----------------------------------------------------------------------
double MMSLinear::computeSigmaT(mfem::Vector &x)
{
   //std::cout << "sigmaT:\t" << 1/mEpsThick << "\n";
   return 1/mEpsThick;
}
//-----------------------------------------------------------------------
double MMSLinear::computeSigmaA(mfem::Vector &x)
{
   //return mEpsThick;
   return this->computeSigmaT(x) - this->computeSigmaS(x);
}
//-----------------------------------------------------------------------
double MMSLinear::computeSigmaS(mfem::Vector &x)
{
   //return this->computeSigmaT(x) - this->computeSigmaA(x);
   return this->computeSigmaT(x)*0.999;
}
//-----------------------------------------------------------------------
double MMSLinear::computeD(mfem::Vector &x)
{
   return 1/(this->computeSigmaT(x)*3);
}
//-----------------------------------------------------------------------
double MMSLinear::computeDiffQ(mfem::Vector &x)
{
   return (4*M_PI)*(mVecDir[0][q]*mAlpha + mVecDir[1][q]*mBeta
            + computeSigmaT(x)*computeInflow(x)
            - computeSigmaS(x)*computeInflow(x));
}
//-----------------------------------------------------------------------
double MMSLinear::computeInflow(mfem::Vector &x)
{
   return mAlpha*x[0] + mBeta*x[1] + mGamma;
}
//-----------------------------------------------------------------------
double MMSLinear::computeAnalyticSolution(mfem::Vector &x)
{
   // phi = 4*M_PI * (ax + by + g)
   return 4*M_PI*computeInflow(x);
}
//-----------------------------------------------------------------------
double MMSLinear::computeInitialCond(mfem::Vector &x)
{
   return 0.0;
}

//-----------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////

MMSWarsa::MMSWarsa(std::vector<std::vector<double> > vecDir,
                           const double epsThick,
                           double alpha,
                           double beta)
    : ProblemSpecification(),
      mEpsThick(epsThick),
      mVecDir(vecDir),
      mAlpha(alpha),
      mBeta(beta)
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
double MMSWarsa::computeS0(mfem::Vector& x)
{
   return this->computeDiffQ(x) / (4*M_PI);
}
//-----------------------------------------------------------------------
double MMSWarsa::computeSigmaT(mfem::Vector &x)
{
   return 1/mEpsThick;
}
//-----------------------------------------------------------------------
double MMSWarsa::computeSigmaA(mfem::Vector &x)
{
   //return mEpsThick;
   return this->computeSigmaT(x) - this->computeSigmaS(x);
}
//-----------------------------------------------------------------------
double MMSWarsa::computeSigmaS(mfem::Vector &x)
{
   //return this->computeSigmaT(x) - this->computeSigmaA(x);
   return this->computeSigmaT(x)*0.9999;
}
//-----------------------------------------------------------------------
double MMSWarsa::computeD(mfem::Vector &x)
{
   return 1/(this->computeSigmaT(x)*3);
}
//-----------------------------------------------------------------------
double MMSWarsa::computeDiffQ(mfem::Vector &x)
{
   return 4.*M_PI * (mVecDir[1][q] * (2*mAlpha*x[0]-4*x[0]*x[0]*x[0]) * (mBeta*x[1]*x[1]-x[1]*x[1]*x[1]*x[1]) * (1+mVecDir[1][q]*mVecDir[1][q]+mVecDir[0][q]*mVecDir[0][q])
      + mVecDir[0][q] * (mAlpha*x[0]*x[0]-x[0]*x[0]*x[0]*x[0]) * (2*mBeta*x[1]-4*x[1]*x[1]*x[1]) * (1+mVecDir[1][q]*mVecDir[1][q]+mVecDir[0][q]*mVecDir[0][q])
      + this->computeSigmaT(x)*this->computeInflow(x)
      - 1./(4.*M_PI)*this->computeSigmaS(x)*this->computeAnalyticSolution(x));
}
//-----------------------------------------------------------------------
double MMSWarsa::computeInflow(mfem::Vector &x)
{
   // psi = x^2 * y^2 * (alpha - x^2)*(beta - y^2) * (1 + mu^2 + eta^2)
   return x[0]*x[0] * x[1]*x[1] * (mAlpha-x[0]*x[0]) * (mBeta-x[1]*x[1]) * (1+mVecDir[1][q]*mVecDir[1][q]+mVecDir[0][q]*mVecDir[0][q]);
}
//-----------------------------------------------------------------------
double MMSWarsa::computeAnalyticSolution(mfem::Vector &x)
{
   // phi = 20*pi/3 * x^2 * y^2 * (alpha - x^2) * (beta - y^2)
   return 20*M_PI/3 * x[0]*x[0] * x[1]*x[1] * (mAlpha-x[0]*x[0]) * (mBeta-x[1]*x[1]);
}
//-----------------------------------------------------------------------
double MMSWarsa::computeInitialCond(mfem::Vector &x)
{
   return 0.0;
}

//-----------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////

MMSWoodsMS::MMSWoodsMS(std::vector<std::vector<double> > vecDir,
                           const double epsThick,
                           double alpha,
                           double beta,
                           double gamma,
                           double delta)
    : ProblemSpecification(),
      mEpsThick(epsThick),
      mVecDir(vecDir),
      mAlpha(alpha),
      mBeta(beta),
      mGamma(gamma),
      mDelta(delta)
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
double MMSWoodsMS::computeS0(mfem::Vector& x)
{
   return this->computeDiffQ(x) / (4*M_PI);
}
//-----------------------------------------------------------------------
double MMSWoodsMS::computeSigmaT(mfem::Vector &x)
{
   return 10.0;
   return 1/mEpsThick;
}
//-----------------------------------------------------------------------
double MMSWoodsMS::computeSigmaA(mfem::Vector &x)
{
   //return mEpsThick;
   return this->computeSigmaT(x) - this->computeSigmaS(x);
}
//-----------------------------------------------------------------------
double MMSWoodsMS::computeSigmaS(mfem::Vector &x)
{
   //return this->computeSigmaT(x) - this->computeSigmaA(x);
   return this->computeSigmaT(x)*0.3;
}
//-----------------------------------------------------------------------
double MMSWoodsMS::computeD(mfem::Vector &x)
{
   return 1/(this->computeSigmaT(x)*3);
}
//-----------------------------------------------------------------------
double MMSWoodsMS::computeDiffQ(mfem::Vector &x)
{
   return 4.*M_PI * (mVecDir[1][q] * 3*M_PI * mDelta * cos(3.*M_PI*x[0])*cos(4.*M_PI*x[1])
      - mVecDir[0][q] * 4*M_PI * mDelta * sin(3.*M_PI*x[0])*sin(4.*M_PI*x[1])
      + this->computeSigmaT(x)*this->computeInflow(x)
      - 1./(4.*M_PI)*this->computeSigmaS(x)*this->computeAnalyticSolution(x));
}
//-----------------------------------------------------------------------
double MMSWoodsMS::computeInflow(mfem::Vector &x)
{
   // psi = a + b*mu + g*eta + d*sin(3*pi*x)*cos(4*pi*y)
   return mAlpha + mBeta*mVecDir[1][q] + mGamma*mVecDir[0][q] + mDelta*sin(3.*M_PI*x[0])*cos(4.*M_PI*x[1]);
}
//-----------------------------------------------------------------------
double MMSWoodsMS::computeAnalyticSolution(mfem::Vector &x)
{
   // phi = 4*pi * (a + d*sin(3*pi*x)*cos(4*pi*y))
   return 4*M_PI * (mAlpha + mDelta*sin(3*M_PI*x[0])*cos(4*M_PI*x[1]));
}
//-----------------------------------------------------------------------
double MMSWoodsMS::computeInitialCond(mfem::Vector &x)
{
   return 0.0;
}

//-----------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////

MultiMaterial::MultiMaterial(std::vector<std::vector<double> > vecDir)
    : ProblemSpecification(),
      mVecDir(vecDir)
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
double MultiMaterial::computeS0(mfem::Vector& x)
{
   // S0 = Q / (4*pi)
   return this->computeDiffQ(x) / (4*M_PI);
}
//-----------------------------------------------------------------------
double MultiMaterial::computeSigmaT(mfem::Vector &x)
{
   if (x[0] < 0.2)
   {
      return 1.0;
   }
   else if (x[0] < 0.4 && x[1] > 0.2)
   {
      return 10.0;
   }
   else if (x[0] < 0.6 && x[1] < 0.2)
   {
      return 0.0001;
   }
   else if (x[1] > 0.2)
   {
      return 1000.0;
   }
   else if (x[1] < 0.2)
   {
      return 100.0;
   }
}
//-----------------------------------------------------------------------
double MultiMaterial::computeSigmaA(mfem::Vector &x)
{
   return this->computeSigmaT(x) - this->computeSigmaS(x);
}
//-----------------------------------------------------------------------
double MultiMaterial::computeSigmaS(mfem::Vector &x)
{
   if (x[0] < 0.2)
   {
      return 1.0;
   }
      else if (x[0] < 0.4 && x[1] > 0.2)
   {
      return 0.0;
   }
   else if (x[0] < 0.6 && x[1] < 0.2)
   {
      return 0.0;
   }
   else if (x[1] > 0.2)
   {
      return 1000.0;
   }
   else if (x[1] < 0.2)
   {
      return 0.0;
   }
}
//-----------------------------------------------------------------------
double MultiMaterial::computeD(mfem::Vector &x)
{
   return 1/(this->computeSigmaT(x)*3);
}
//-----------------------------------------------------------------------
double MultiMaterial::computeDiffQ(mfem::Vector &x)
{
   // Q = S0
   if (x[0] < 0.2)
   {
      return 1.0 / (4*M_PI);
   }
   else
   {
      return 0;
   }
}
//-----------------------------------------------------------------------
double MultiMaterial::computeInflow(mfem::Vector &x)
{
   if (x[0] == 0)
   {
      return 1.0 / (4*M_PI);
   }
   else
   {
      return 0;
   }
}
//-----------------------------------------------------------------------
double MultiMaterial::computeAnalyticSolution(mfem::Vector &x)
{
   return 0.0;
}
//-----------------------------------------------------------------------
double MultiMaterial::computeInitialCond(mfem::Vector &x)
{
   return 0.0;
}

//-----------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////

WangRagusa::WangRagusa(std::vector<std::vector<double> > vecDir,
                           const double epsThick)
    : ProblemSpecification(),
      mEpsThick(epsThick),
      mVecDir(vecDir)
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
double WangRagusa::computeS0(mfem::Vector& x)
{
   // S0 = Q / (4*pi)
   return this->computeDiffQ(x) / (4*M_PI);
}
//-----------------------------------------------------------------------
double WangRagusa::computeSigmaT(mfem::Vector &x)
{
   return 1 / mEpsThick;
}
//-----------------------------------------------------------------------
double WangRagusa::computeSigmaA(mfem::Vector &x)
{
   return this->computeSigmaT(x) - this->computeSigmaS(x);
}
//-----------------------------------------------------------------------
double WangRagusa::computeSigmaS(mfem::Vector &x)
{
   return computeSigmaT(x)*0.9999;
}
//-----------------------------------------------------------------------
double WangRagusa::computeD(mfem::Vector &x)
{
   return 1/(this->computeSigmaT(x)*3);
}
//-----------------------------------------------------------------------
double WangRagusa::computeDiffQ(mfem::Vector &x)
{
   // Q = S0
   return 1.0;
}
//-----------------------------------------------------------------------
double WangRagusa::computeInflow(mfem::Vector &x)
{
   // psi_inc = phi / sigma_t = s / (4 * pi * sigma_a)
   return computeDiffQ(x) / (4*M_PI* computeSigmaA(x));
}
//-----------------------------------------------------------------------
double WangRagusa::computeAnalyticSolution(mfem::Vector &x)
{
   return 0.0;
   return computeDiffQ(x) / computeSigmaA(x);
}
//-----------------------------------------------------------------------
double WangRagusa::computeInitialCond(mfem::Vector &x)
{
   return 0.0;
}

//-----------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////

uniformInfMed::uniformInfMed(std::vector<std::vector<double> > vecDir,
                           const double epsThick)
    : ProblemSpecification(),
      mEpsThick(epsThick),
      mVecDir(vecDir)
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
   return 1.0;
   //return 1/mEpsThick;
}
//-----------------------------------------------------------------------
double uniformInfMed::computeSigmaA(mfem::Vector &x)
{
   return this->computeSigmaT(x) - this->computeSigmaS(x);
}
//-----------------------------------------------------------------------
double uniformInfMed::computeSigmaS(mfem::Vector &x)
{
   return 0.3;
   //return this->computeSigmaT(x) - this->computeSigmaA(x);
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
   return 0.7;
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
                             //const std::vector<double> vecDir,
                             std::vector<std::vector<double> > vecDir,
                             std::unique_ptr<ProblemSpecification>& problemSpec)
{
   if (problemName == "Adams")
   {
      //std::cout << "The problem name is: " << problemName <<"\n";
      problemSpec.reset(new AdamsProblem(vecDir, epsThick));
   }
   else if (problemName == "AlternatingIncident")
   {
      std::cout << "The problem name is: " << problemName << "\n";
      std::cout << "*** Don't forget to select the correct mesh ***\n";
      problemSpec.reset(new AlternatingIncident(vecDir, epsThick));
   }
   else if (problemName == "MMSGleicher")
   {
      std::cout << "The problem name is: " << problemName << "\n";
      const int probSelection = 2;
      double alpha = 0;
      double beta = 0;
      if (probSelection == 1)
      {
         alpha = 1.0;
         beta = 1.0/2;
      }
      else if (probSelection == 2)
      {
         alpha = 2.0;
         beta = 3.0/2;
      }
      else if (probSelection == 3)
      {
         alpha = 3.0;
         beta = 5.0/2;
      }
      else if (probSelection == 4)
      {
         alpha = 4.0;
         beta = 7.0/2;
      }
      problemSpec.reset(new MMSGleicher(vecDir, epsThick, alpha, beta));
   }
   else if (problemName == "MMSLinear")
   {
      std::cout << "The problem name is: " << problemName << "\n";
      const double alpha = 1;
      const double beta = 1;
      const double gamma = 1;
      problemSpec.reset(new MMSLinear(vecDir, epsThick, alpha, beta, gamma));
   }
   else if (problemName == "MMSWarsa")
   {
      std::cout << "The problem name is: " << problemName << "\n";
      const double alpha = 1.1;
      const double beta = 1.25;
      problemSpec.reset(new MMSWarsa(vecDir, epsThick, alpha, beta));
   }
   else if (problemName == "MMSWoodsMS")
   {
      std::cout << "The problem name is: " << problemName << "\n";
      const double alpha = 10.;
      const double beta = 5.;
      const double gamma = 1.;
      const double delta = 5.;
      problemSpec.reset(new MMSWoodsMS(vecDir, epsThick, alpha, beta, gamma, delta));
   }
   else if (problemName == "MultiMaterial")
   {
      std::cout << "The problem name is: " << problemName << "\n";
      std::cout << "*** Don't forget to select the correct mesh ***\n";
      problemSpec.reset(new MultiMaterial(vecDir));
   }
   else if (problemName == "WangRagusa")
   {
      std::cout << "The problem name is: " << problemName << "\n";
      std::cout << "*** Don't forget to select the correct mesh ***\n";
      problemSpec.reset(new WangRagusa(vecDir, epsThick));
   }
   else if (problemName == "uniformInfMed")
   {
      std::cout << "The problem name is: " << problemName << "\n";
      problemSpec.reset(new uniformInfMed(vecDir, epsThick));
   }
   else
   {
      MFEM_ABORT("Unsupported problem specification: " << problemName);
   }
}

///////////////////////////////////////////////////////////////////////////
