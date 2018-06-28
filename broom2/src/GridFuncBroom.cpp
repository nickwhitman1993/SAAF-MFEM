#include "GridFuncBroom.hpp"
#include "mfem.hpp"

/// This is identical to the mfem::GridFunction version, but this allows for
//    tweaking and toying (to see if we are getting the real L2 norm).
double GridFunctionBroom::ComputeLpErrorBroom2(const double p, mfem::Coefficient &exsol,
                                    mfem::Coefficient *weight,
                                    const mfem::IntegrationRule *irs[]) const
{
   // Test getting the mesh order
   //std::cout << "mesh order:\t" << fes->GetMesh()->GetElementTransformation(0)->Order() << "\n";
   // end test
   
   double error = 0.0;
   const mfem::FiniteElement *fe;
   mfem::ElementTransformation *T;
   mfem::Vector vals;

   // Number of integration points as usual.
   //std::cout << "fes->GetNE():\t" << fes->GetNE() << "\n";
   for (int i = 0; i < fes->GetNE(); i++)
   {
      fe = fes->GetFE(i);
      const mfem::IntegrationRule *ir;
      if (irs)
      {
         ir = irs[fe->GetGeomType()];
      }
      else
      {
         // The FEM integration order we use is:
         //       2 * feOrder + meshOrder * dim - 1
         // Should this be the same?
         int meshOrder = fes->GetMesh()->GetElementTransformation(0)->Order();
         //std::cout << "meshOrder:\t" << meshOrder << "\n";
         int dim = fes->GetMesh()->Dimension();
         //std::cout << "dim:\t" << dim << "\n";
         //int intorder = 2*fe->GetOrder() + meshOrder*dim - 1;
         //intorder = pow(10, dim);
         //std::cout << "intorder:\t" << intorder << "\n";
         int intorder = 2*fe->GetOrder() + 1; // <----------
         ir = &(mfem::IntRules.Get(fe->GetGeomType(), intorder));
      }
      GetValues(i, *ir, vals); // where does it get the vals (phi)?
      T = fes->GetElementTransformation(i);
      
      // Number of integration points per element.
      //std::cout << "ir->GetNPoints():\t" << ir->GetNPoints() << "\n";
      for (int j = 0; j < ir->GetNPoints(); j++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(j);
         double d[2]={0,0};
         ip.Get(d, 2);
         //std::cout << "ip (x,y):\t" << d[0] << "\t" << d[1] << "\n";
         T->SetIntPoint(&ip);
         double err = fabs(vals(j) - exsol.Eval(*T, ip));
         if (p < std::numeric_limits<double>::infinity())
         {
            err = pow(err, p);
            if (weight)
            {
               err *= weight->Eval(*T, ip);
            }
            // ip.weight is the integration point weight.
            // T->Weight() is the dx of the area
            error += ip.weight * T->Weight() * err;
            //std::cout << "ir->GetNPoints():\t" << ir->GetNPoints() << "\tip.weight:\t" << ip.weight << "\tT->Weight():\t" << T->Weight() << "\n";
         }
         else
         {
            if (weight)
            {
               err *= weight->Eval(*T, ip);
            }
            error = std::max(error, err);
         }
      }
   }

   if (p < std::numeric_limits<double>::infinity())
   {
      // negative quadrature weights may cause the error to be negative
      if (error < 0.)
      {
         error = -pow(-error, 1./p);
      }
      else
      {
         error = pow(error, 1./p);
      }
   }

   return error;
}
