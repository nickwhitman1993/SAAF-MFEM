#include "mfem.hpp"
#include "physics.hpp"
#include "diffSolve.hpp"

diffSolve::diffSolve(int dim,
            double (*sigma_a_function)(mfem::Vector &x),
            double (*diff_Q_function)(mfem::Vector &x),
            double (*diff_function)(mfem::Vector &x),
            mfem::FiniteElementSpace &fes,
            mfem::GridFunction &diffPhi,
            bool directSolve,
            int myRank)
      : sigma_a_vis(&fes), diff(&fes), DiffCoeff(&fes)
{

   bool DirichletBC(true);
   bool IncidentCurrentBC(false);

   int feOrder = fes.GetOrder(0);
   int C = 2;
   double kappa = C * feOrder * (feOrder + 1);
   double kappab = C * feOrder * (feOrder + 1);
   
   // kappa = 1/4 * h/D = 1/4 * h*3*sig_t
   //kappa = 0.25*(1./12.)*3*1/1e-1;
   //kappa = 1;
   kappab = kappa;
   
   double sigma = -1.0;
   double sigmab = -1.0;
   
   
   mfem::ConstantCoefficient zeroCoef(0.0);
   mfem::ConstantCoefficient oneCoef(1.0);
   mfem::ConstantCoefficient halfCoef(0.5); // for Method 3
   mfem::ConstantCoefficient quartCoef(0.25); // for Methods 1,2
   mfem::FunctionCoefficient diffCoef(diff_function);
   mfem::FunctionCoefficient qCoef(diff_Q_function);
   mfem::FunctionCoefficient sigmaACoef(sigma_a_function);
   //mfem::FunctionCoefficient inflow(qext_function);
   //std::cout << "Created functionCoefficient\n";
   
   //std::cout << "elementGeometry\n";
   int elementGeometry = fes.GetMesh()->GetElement(0)->GetGeometryType();
   // GetBdrElement(0) screams if there aren't actually any Bdrs (periodic mesh)
   int checkBE = fes.GetMesh()->GetNBE();
   //std::cout << "checkBE: " << checkBE << "\n";
   int borderElementGeometry;
   if(checkBE > 0)
   {
      borderElementGeometry =
         fes.GetMesh()->GetBdrElement(0)->GetGeometryType();
   }
   
   //std::cout << "Get transformed mesh\n";
   int meshOrder = fes.GetMesh()->GetElementTransformation(0)->Order();
   int elementIntegrationOrder = 0;
   if (feOrder == 0) {
      elementIntegrationOrder = 2 * meshOrder;
   }
   else if(meshOrder == 0) {
      elementIntegrationOrder = 2 * feOrder;
   } else {
      elementIntegrationOrder = 2 * (feOrder)+dim * meshOrder - 1;
   }
   int borderElementIntegrationOrder = elementIntegrationOrder;
   
   //std::cout << "integration orders\n";
   const mfem::IntegrationRule &elIntRule =
      mfem::IntRules.Get(elementGeometry, elementIntegrationOrder);
   const mfem::IntegrationRule &bdrIntRule = mfem::IntRules.Get(elementGeometry, elementIntegrationOrder);
   if(checkBE > 0)
   {
      const mfem::IntegrationRule &bdrIntRule =
         mfem::IntRules.Get(borderElementGeometry, borderElementIntegrationOrder);
   }
   
    // Construct system of equations
    //std::cout << "Construct system of equations\n";

    mfem::LinearForm *diffRHS = new mfem::LinearForm(&fes);
    mfem::BilinearForm *diffMatrix = new mfem::BilinearForm(&fes);
    
    // RHS = diff_Q_function
    mfem::LinearFormIntegrator *LFinteriorFlux = new mfem::DomainLFIntegrator(qCoef);
    LFinteriorFlux->SetIntRule(&elIntRule);
    diffRHS->AddDomainIntegrator(LFinteriorFlux);
   
   if(DirichletBC) {
    // For imposing non-zero Dirichlet Bdry Conds.
    // TODO: there is an issue with non-zero Dirichlet BCs if problem is
    //       optically thick and diffusive. Consider modifying
    //       DGDirichletLFIntegrator to something similar to MIP "switch".
    mfem::LinearFormIntegrator *DirichletBdr = new mfem::DGDirichletLFIntegrator( zeroCoef, diffCoef, sigmab, kappab);
    DirichletBdr->SetIntRule(&bdrIntRule);
    diffRHS->AddBdrFaceIntegrator(DirichletBdr);
    }
    
   // LHS = -D*grad.grad.f + sig_a*f; solve for f
   mfem::BilinearFormIntegrator *diffInterior = new mfem::DiffusionIntegrator(diffCoef);
   diffInterior->SetIntRule(&elIntRule);
   diffMatrix->AddDomainIntegrator(diffInterior);
   
   mfem::BilinearFormIntegrator *mass = new mfem::MassIntegrator(sigmaACoef);
   mass->SetIntRule(&elIntRule);
   diffMatrix->AddDomainIntegrator(mass);
   
   mfem::BilinearFormIntegrator *DGInterior = new mfem::DGDiffusionMIPIntegrator(diffCoef, sigma, C);
   DGInterior->SetIntRule(&elIntRule);
   diffMatrix->AddInteriorFaceIntegrator(DGInterior);
   
   if(DirichletBC) {
   // Dirichlet BCs (may need to add linear term)
   mfem::BilinearFormIntegrator *DGBdr = new mfem::DGDiffusionMIPIntegrator(diffCoef, sigmab, C);
   DGBdr->SetIntRule(&bdrIntRule);
   diffMatrix->AddBdrFaceIntegrator(DGBdr);
   }
   
    if(IncidentCurrentBC) {
    // Incident Current BC for Methods 1,2 only
    mfem::BilinearFormIntegrator *DGDiffusionCurrent = new mfem::DGCurrentIntegrator(diffCoef, sigmab, kappab);
    DGDiffusionCurrent->SetIntRule(&bdrIntRule);
    //diffMatrix->AddBdrFaceIntegrator(DGDiffusionCurrent);
   
   // For imposing incident current BCs (add linear term if not zero current)

    mfem::BilinearFormIntegrator *DGCurrent = new mfem::BoundaryMassIntegrator(halfCoef);
    DGCurrent->SetIntRule(&bdrIntRule);
    diffMatrix->AddBdrFaceIntegrator(DGCurrent);
    }
    
    //std::cout << "Assemble\n";
    diffRHS->Assemble();
    diffMatrix->Assemble();
    diffMatrix->Finalize();
    //std::cout << diffMatrix->Size() << "\n";
    
    if(myRank==0){
    //diffMatrix->PrintMatlab(std::cout, 64,64);
    }
    
    const mfem::SparseMatrix &DiffA = diffMatrix->SpMat();
    
    if(myRank == 0) {
      //DiffA.PrintMatlab();
    }
    
   // Solve the system
   //std::cout << "Solve\n";
   mfem::UMFPackSolver solver;
   solver.SetOperator(DiffA);
   solver.Mult(*diffRHS, diffPhi);
   
}
