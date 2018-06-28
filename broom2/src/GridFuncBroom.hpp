#ifndef GRIDFUNCBROOM_HH__
#define GRIDFUNCBROOM_HH__

#include "mfem.hpp"

class GridFunctionBroom: public mfem::GridFunction
{
   protected:
      mfem::FiniteElementSpace *fes2;
      mfem::FiniteElementCollection *fec2;
      long sequence2;
   public:
   /// Copy Constructor
   GridFunctionBroom(const mfem::GridFunction &orig)
      : GridFunction(orig) { }
   
   virtual double ComputeL2ErrorBroom2(mfem::Coefficient &exsol,
                                 const mfem::IntegrationRule *irs[] = NULL) const
   { return ComputeLpErrorBroom2(2.0, exsol, NULL, irs); }
   
   virtual double ComputeLpErrorBroom2(const double p, mfem::Coefficient &exsol,
                                 mfem::Coefficient *weight = NULL,
                                 const mfem::IntegrationRule *irs[] = NULL) const;

};
#endif // GRIDFUNCTIONBROOM_HH__
