#include "mfem.hpp"
#include "physics.hpp"
#include "Solve.hpp"

Solve::Solve(int dim, int &qloop, int &q,
             void (*omega_function)(const mfem::Vector &x, mfem::Vector &v),
             double (*inflow_function)(mfem::Vector &x),
             double (*psi0_function)(mfem::Vector &x),
             double (*sigma_t_function)(mfem::Vector &x),
             double (*sigma_s_function)(mfem::Vector &x),
             double (*S0_function)(mfem::Vector &x),
             mfem::FiniteElementSpace &fes, mfem::GridFunction &phi_old,
             mfem::GridFunction &phi, int myRank, int numAngles, int numProcs,
             bool directSolve,
             std::unique_ptr<ProblemSpecification> &problemSpec)
    : psi(&fes), sigma_t_vis(&fes), sigma_s_vis(&fes), scat_src(&fes),
      scat_src_vis(&fes), qext_vis(&fes), psi0_vis(&fes), inflow_vis(&fes) {

  q = myRank * (numAngles / numProcs) + qloop;
  //std::cout << "q:\t" << q << "\tmyRank:\t" << myRank << "\n";
  // 5. Set up and assemble the bilinear and linear forms corresponding to
  // the
  //    DG discretization. The DGTraceIntegrator involves integrals over
  //    mesh
  //    interior faces.
  mfem::VectorFunctionCoefficient omega(dim, omega_function);
  mfem::FunctionCoefficient inflow(inflow_function);
  mfem::FunctionCoefficient psi0(psi0_function);
  mfem::FunctionCoefficient sigma_t(sigma_t_function);

  mfem::BilinearForm L(&fes);
  mfem::LinearForm b(&fes);

  // Make Integrators for the bilinear form.  The L object will
  // take ownership of the pointers and delete them.  The integrators
  // DO NOT take ownership of the integration rule, so we can safely
  // pass one to multiple integrators.

  // Assume that all zones in the mesh are the same type (all quads, triangles,
  // etc.)
  // We find the geometry type by looking at the first element and first border
  // element.
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

  // From the guts of all the different integrators showing their defaults.
  // For Quads and hexes, less for triangles/tets.  Line numbers for mfem
  // version 6472008bca376a6083ccbdc02fb1ada07e14d0d1
  // k = mesh-order
  // d = element dimension
  // l = finite element order
  // OrderW = meshOrder * dim - 1                      fem/eltrans.cpp, 103
  // OrderGrad = meshOrder * (dim-1) + (feOrder-1)     fem/eltrans.cpp, 122
  // OrderJ = meshOrder                                fem/eltrans.cpp, 89
  // Order = meshOrder                                 fem/eltrans.hpp, 78
  // a = 2, b = 0 (user set)                           fem/lininteg.hpp, 54

  // mass: 2 * feOrder + OrderW                        fem/bilininteg.cpp, 501
  //           = 2 * feOrder + meshOrder * dim - 1
  // convection: OrderGrad + Order + feOrder           fem/bilininteg.cpp, 587
  //           = meshOrder * dim + 2 * feOrder - 1
  // trace: OrderW + 2 * feOrder                       fem/bilininteg.cpp, 1564
  //           = meshOrder * dim -1 + 2 * feOrder
  // boundary: OrderW + 2 * feOrder                    fem/lininteg.cpp, 392
  //           = meshOrder * dim - 1 + 2 * feOrder
  // domain(source): a * meshOrder + b                 fem/lininteg.cpp, 42
  //           = 2 * meshOrder
  //
  //
  // DIFFUSION OPERATORS
  //
  // DGDiffusion: 2 * feOrder                          fem/bilininteg.cpp
  // Diffusion: 2 * feOrder + dim - 1                  fem/bilininteg.cpp

  int meshOrder = fes.GetMesh()->GetElementTransformation(0)->Order();
  int feOrder = fes.GetOrder(0);
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

  // MFEM ISSUE: Can we set this on the integrator?  Or even the bilinear or
  // linear forms?
  // Even better, specify elementIntegrationOrder, and let it get the rule for
  // the geometry that
  // exists, instead of assuming they are all the same on the mesh.
  const mfem::IntegrationRule &elIntRule =
      mfem::IntRules.Get(elementGeometry, elementIntegrationOrder);
   const mfem::IntegrationRule &bdrIntRule = mfem::IntRules.Get(elementGeometry, elementIntegrationOrder);
   if(checkBE > 0)
   {
      const mfem::IntegrationRule &bdrIntRule =
         mfem::IntRules.Get(borderElementGeometry, borderElementIntegrationOrder);
   }

  // Add the removal term from the zones.
  mfem::BilinearFormIntegrator *mass = new mfem::MassIntegrator(sigma_t);
  mass->SetIntRule(&elIntRule);
  L.AddDomainIntegrator(mass);

  // Add convection term within the zones
  mfem::BilinearFormIntegrator *convection =
      new mfem::ConvectionIntegrator(omega, +1.0);
  convection->SetIntRule(&elIntRule);
  L.AddDomainIntegrator(convection);

  // Upwind fluxes between zones on each face interior to mesh
  // For surface, we need a boundary geometry integration rule and order
  mfem::BilinearFormIntegrator *traceInterior =
      new mfem::DGTraceIntegrator(omega, -1.0, +0.5);
  traceInterior->SetIntRule(&bdrIntRule);

  // MFEM ISSUE: Maybe have the transpose operator pass along the setting of the
  // IntRule
  // to its owned worker.  Or at least throw an error.
  mfem::BilinearFormIntegrator *transInterior =
      new mfem::TransposeIntegrator(traceInterior);
  L.AddInteriorFaceIntegrator(transInterior);

  // Now add the inflow boundary conditions to the matrix and RHSs
  // need to SetIntRule on the DGTraceIntegrator rather than the
  // TransposeIntegrator
  mfem::BilinearFormIntegrator *traceBdr =
      new mfem::DGTraceIntegrator(omega, -1.0, +0.5);
  traceBdr->SetIntRule(&bdrIntRule);
  mfem::BilinearFormIntegrator *transBdr =
      new mfem::TransposeIntegrator(traceBdr);
  L.AddBdrFaceIntegrator(transBdr);

   mfem::Array<int> in_bdr(4);
   in_bdr[0] = 1; // bottom
   in_bdr[1] = 1; // right
   in_bdr[2] = 1; // top
   in_bdr[3] = 1; // left
  mfem::LinearFormIntegrator *BndryFlow =
      new mfem::BoundaryFlowIntegrator(inflow, omega, -1.0, -0.5);
  BndryFlow->SetIntRule(&bdrIntRule);
  b.AddBdrFaceIntegrator(BndryFlow); // MFEM 3.2
  //b.AddBdrFaceIntegrator(BndryFlow, in_bdr); // MFEM 3.3
  inflow_vis.ProjectCoefficient(inflow);

  // SOURCE TERM //
  // TODO: Move this term outside the angle loop and pass it in. Only need to calc once per source iteration.
  mfem::FunctionCoefficient sigma_s(sigma_s_function);
  sigma_s_vis.ProjectCoefficient(sigma_s);
  for (int i = 0; i < phi_old.Size(); i++) {
    scat_src[i] = phi_old[i] * sigma_s_vis[i] / (4.0 * M_PI);
  }
  // Why are we taking the address of this?  It should take a reference.

  mfem::GridFunctionCoefficient flux(&scat_src);
  scat_src_vis.ProjectCoefficient(flux);
  
  mfem::LinearFormIntegrator *LFinteriorFlux =
      new mfem::DomainLFIntegrator(flux);
  LFinteriorFlux->SetIntRule(&elIntRule);
  b.AddDomainIntegrator(LFinteriorFlux);
     
   mfem::FunctionCoefficient qext(S0_function);
   qext_vis.ProjectCoefficient(qext);
   
  mfem::LinearFormIntegrator *LFinteriorQext =
      new mfem::DomainLFIntegrator(qext);
  LFinteriorQext->SetIntRule(&elIntRule);
  b.AddDomainIntegrator(LFinteriorQext);

  // Assemble and build the matrix and RHS
  L.Assemble();
  L.Finalize();
  b.Assemble();

  // 6. Define the initial conditions, save the corresponding grid
  // function to
  //    a file and (optionally) save data in the VisIt format and
  //    initialize
  //    GLVis visualization.
  psi.ProjectCoefficient(psi0);
  psi0_vis.ProjectCoefficient(psi0);

  sigma_t_vis.ProjectCoefficient(sigma_t);

  // Get the actual matrix
  const mfem::SparseMatrix &A = L.SpMat();

// 7. Solve the system
#if defined(MFEM_USE_SUITESPARSE)
  if (!directSolve) {
    // Use KLU as a preconditioner in the GMRES solve
    //mfem::KLUSolver M;
    //M.SetOperator(A);
//#else
#endif
    mfem::DSmoother M(
        A, 1); // l1-Jacobi, we can try l1-GS and polynomial in parallel
//#endif
    int maxit = 1000;
    double tol = 1e-14;
    double atol = 0.0;
    aGMRES(A, psi, b, M, maxit, 20, 1, 5, 0.4, tol, atol, 1);
#if defined(MFEM_USE_SUITESPARSE)
  } else {
    //mfem::KLUSolver klu_solver;
    //klu_solver.SetOperator(A);
    //klu_solver.Mult(b,psi);
    mfem::UMFPackSolver umf_solver;
    umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
    umf_solver.SetOperator(A);
    umf_solver.Mult(b, psi);
  }
#endif

  // ADD ANGULAR FLUX TO SCALAR FLUX
  for (int i = 0; i < psi.Size(); i++) {
    phi[i] += w[q] * psi[i];
  }
}
