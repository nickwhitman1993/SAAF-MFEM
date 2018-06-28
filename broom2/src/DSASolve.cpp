#include "mfem.hpp"
#include "physics.hpp"
#include "DSASolve.hpp"
#include "HPDouble.hh"
#include "HPDoubleOps.hh"

DSASolve::DSASolve(int dim,
            double (*sigma_s_function)(mfem::Vector &x),
            double (*sigma_a_function)(mfem::Vector &x),
            double (*diff_function)(mfem::Vector &x),
            mfem::FiniteElementSpace &fes,
            mfem::GridFunction &phi_old,
            mfem::GridFunction &lDsaPhi,
            mfem::GridFunction &gphi,
            bool directSolve,
            int myRank)
      : error(&fes), sigma_s_vis(&fes), sigma_a_vis(&fes), diff(&fes), DiffCoeff(&fes), scatCorr(&fes), diff_vis(&fes)
{
   
   bool DirichletBC(true);
   bool IncidentCurrentBC(false);
   
   const double scale = 1.;

   int feOrder = fes.GetOrder(0);
   int C = 2;
   double kappa = C * feOrder * (feOrder + 1);
   double kappab = C * feOrder * (feOrder + 1);
   double sigma = -1.0;
   double sigmab = -1.0;
   
   mfem::FunctionCoefficient diffCoef(diff_function);
   diff_vis.ProjectCoefficient(diffCoef);
   //kappa = 150; // TP2
   //kappa = 4.5e2; // TP3 --- needs varying kappa
   //kappa = 0.25*(1./16.)/diff_vis[0]; // (1/4 * h/D)
   //kappa = 4.6875e1;
   //kappa = 1;
   kappab = kappa;
   
   
   mfem::ConstantCoefficient zeroCoef(0.0);
   mfem::ConstantCoefficient oneCoef(1.0);
   mfem::ConstantCoefficient halfCoef(0.5);
   mfem::FunctionCoefficient sigmaSCoef(sigma_s_function);
   mfem::FunctionCoefficient sigmaACoef(sigma_a_function);
   
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
   elementIntegrationOrder += 0;
   int borderElementIntegrationOrder = elementIntegrationOrder;
   
   const mfem::IntegrationRule &elIntRule =
      mfem::IntRules.Get(elementGeometry, elementIntegrationOrder);
   const mfem::IntegrationRule &bdrIntRule = mfem::IntRules.Get(elementGeometry, elementIntegrationOrder);
   if(checkBE > 0)
   {
      const mfem::IntegrationRule &bdrIntRule =
         mfem::IntRules.Get(borderElementGeometry, borderElementIntegrationOrder);
   }

    // Diffusion Solve

    mfem::LinearForm *diffRHS = new mfem::LinearForm(&fes);
    mfem::BilinearForm *diffMatrix = new mfem::BilinearForm(&fes);
    
    
#if 1
   // basic scatCorr calculation
   sigma_s_vis.ProjectCoefficient(sigmaSCoef);
   subtract(gphi, phi_old, error);
   for(int i=0; i<gphi.Size(); i++) {
      scatCorr[i] = error[i]*sigma_s_vis[i];
   }
#endif
    
#if 0
   // Use HPDouble for this subtraction.
   double tmpgphi, tmpphi_old;
   sigma_s_vis.ProjectCoefficient(sigmaSCoef);
   for(int i=0; i<gphi.Size(); i++) {
      tmpgphi = gphi[i]*sigma_s_vis[i]*scale;
      tmpphi_old = phi_old[i]*sigma_s_vis[i]*scale;
      HPDouble HPgphi(tmpgphi);
      HPDouble HPphi_old(tmpphi_old);
      HPDouble HPerror(0.0);
      HPerror -= HPphi_old;
      HPerror += HPgphi;
      scatCorr[i] = HPerror.getValue();
   }
   //scatCorr.Print(std::cout, 8);
#endif
   
#if 0
   // RHS = sig_s * error
   sigma_s_vis.ProjectCoefficient(sigmaSCoef);
   //subtract(tmpgphi, tmpphi_old, error);
   //error.Print(std::cout, 8);
   for(int i=0; i<gphi.Size(); i++) {
      double tmpgphi = gphi[i]*sigma_s_vis[i]*scale;
      double tmpphi_old = phi_old[i]*sigma_s_vis[i]*scale;
      if(tmpgphi < 0) {
         std::cout << std::setprecision(19) << tmpgphi << "\n";
         tmpgphi = fmod(-tmpgphi,1);
         tmpgphi *= -1;
         std::cout << std::setprecision(19) << tmpgphi << "\n";
      }
      else {
         std::cout << std::setprecision(19) << tmpgphi << "\n";
         tmpgphi = fmod(tmpgphi,1);
         std::cout << std::setprecision(19) << tmpgphi << "\n";
      }
      if(tmpphi_old < 0) {
         tmpphi_old = fmod(-tmpphi_old,1);
         tmpphi_old *= -1;
      }
      else {
         tmpphi_old = fmod(tmpphi_old,1);
      }
      scatCorr[i] = (tmpgphi) - (tmpphi_old);
      //scatCorr[i] = error[i] * sigma_s_vis[i];
      //scatCorr[i] = gphi[i]*sigma_s_vis[i] - phi_old[i]*sigma_s_vis[i];
   }
   //scatCorr.Print(std::cout, 8);
   //std::cout << "phi_ol:\t" << std::setprecision(19) <<  phi_old[24] << "\n";
#endif
   //mfem::FunctionCoefficient *rhs = new mfem::GridFunction(&scatCorr);
   // scatCorr needs to be a Vector
   // rhs needs to be a Coefficient:
   mfem::GridFunctionCoefficient rhs(&scatCorr);
   mfem::LinearFormIntegrator *LFinteriorFlux = new mfem::DomainLFIntegrator(rhs);
   LFinteriorFlux->SetIntRule(&elIntRule);
   diffRHS->AddDomainIntegrator(LFinteriorFlux);
   
   if(DirichletBC) {
      // For non-zero Dirichlet BCs
      mfem::LinearFormIntegrator *DirichletBdr = new mfem::DGDirichletLFIntegrator( zeroCoef, diffCoef, sigmab, kappab);
      DirichletBdr->SetIntRule(&bdrIntRule);
      diffRHS->AddBdrFaceIntegrator(DirichletBdr);
   }
    
   // TODO Only need to form this matrix once. Can we save it?
   // LHS = -grad.D grad f + sig_a*f; solve for f
   mfem::BilinearFormIntegrator *diffInterior = new mfem::DiffusionIntegrator(diffCoef);
   diffInterior->SetIntRule(&elIntRule);
   diffMatrix->AddDomainIntegrator(diffInterior);
   
   mfem::BilinearFormIntegrator *mass = new mfem::MassIntegrator(sigmaACoef);
   mass->SetIntRule(&elIntRule);
   diffMatrix->AddDomainIntegrator(mass);
   
   // DGDiffusionMIPIntegrator requires using C instead of kappa.
   mfem::BilinearFormIntegrator *DGInterior = new mfem::DGDiffusionMIPIntegrator(diffCoef, sigma, C);
   DGInterior->SetIntRule(&elIntRule);
   diffMatrix->AddInteriorFaceIntegrator(DGInterior);
   
   if(DirichletBC) {
      // Dirichlet BCs (may need to add linear term)
      mfem::BilinearFormIntegrator *DGBdr = new mfem::DGDiffusionMIPIntegrator(diffCoef, sigmab, C);
      DGBdr->SetIntRule(&bdrIntRule);
      diffMatrix->AddBdrFaceIntegrator(DGBdr);
   }
    
   else if(IncidentCurrentBC) {
      // Incident Current BC
      //mfem::BilinearFormIntegrator *DGDiffusionCurrent = new mfem::DGCurrentIntegrator(diffCoef, sigmab, kappab);
      //DGDiffusionCurrent->SetIntRule(&bdrIntRule);
      //diffMatrix->AddBdrFaceIntegrator(DGDiffusionCurrent);
   
      // For incident current BCs (may need to add linear term)
      mfem::BilinearFormIntegrator *DGCurrent = new mfem::BoundaryMassIntegrator(halfCoef);
      DGCurrent->SetIntRule(&bdrIntRule);
      diffMatrix->AddBdrFaceIntegrator(DGCurrent);
   }
    
   diffRHS->Assemble();
   //diffRHS->Print(std::cout, 8);
   diffMatrix->Assemble();
   diffMatrix->Finalize();
   //diffMatrix->PrintMatlab(std::cout);
   //std::cout << diffMatrix->Size() << "\n";
    
   //if(myRank==0){
      //std::cout << "break\n";
      //diffMatrix->PrintMatlab(std::cout, 256, 256);
      //std::cout << "break\n";
   //}
    
   const mfem::SparseMatrix &DiffA = diffMatrix->SpMat();
    
   if(myRank == 0) {
      //std::cout << "break\n";
      //DiffA.PrintMatlab();
      //std::cout << "break\n";
   }
    
   // Solve the system
   lDsaPhi = 0;
   
   mfem::UMFPackSolver solver;
   solver.SetOperator(DiffA);
   solver.Mult(*diffRHS, lDsaPhi);
   //lDsaPhi /= scale;
    
   if(myRank == 0) {
      for(int i=0;i<3;i++){
         //std::cout << diffPhi[i] << "\n";
      }
   }
}
