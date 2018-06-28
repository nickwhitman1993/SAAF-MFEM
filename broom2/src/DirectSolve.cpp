//                                MFEM **NAME**
//
// Compile with: make ex13
//
// Sample runs:
//    ex11 -m ~/broom/Meshes/inline-quad.mesh
//
// These keep the number of degrees of freedom fixed while decreasing
// the element count and increasing the feOrder of the mesh and discretization.
//    ex11 -m ~/broom/Meshes/zsine-q1-r5.mesh -o 1 -r 2
//    ex11 -m ~/broom/Meshes/zsine-q3-r4.mesh -o 3 -r 2
//    ex11 -m ~/broom/Meshes/zsine-q7-r3.mesh -o 7 -r 2
//
// Description:  This example code solves the steady state linear advection
//               equation omega.grad(u) +sigma_t psi = 0, where omega is the
//               direction of photon travel and sigma_t is the total removal
//               term.
//
//               The example demonstrates the use of Discontinuous Galerkin (DG)
//               bilinear forms in MFEM (face integrators) with a linear solve
//               for an advection-like equation.

#include "mfem.hpp"
#include "physics.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <sstream>
#include <cstring>
#include <vector>
#include "BroomConfig.hh"
#include "Quadrature.hpp"
#include "physics.hpp"
#include "MeshLoader.hh"
#include "tryOptions.hh"
#include <stdio.h>
#include <cstdio>
#include "diffSolve.hpp"
#include "klu.h"
#if defined(BROOM_USE_MPI)
#include <mpi.h>
#endif

#if defined(__linux__)
#include <fenv.h>
#endif

using namespace std;
using namespace mfem;

// The discrete ordinates direction.  It is really constant.
void omega_function(const Vector &x, Vector &v);

// Initial condition
double psi0_function(Vector &x);
double sigma_t_function(Vector &x);
double sigma_s_function(Vector &x);
double qext_function(Vector &x);
double sigma_a_function(Vector &x);
double diff_Q_function(Vector &x);
double diff_function(Vector &x);

double probEpsilon;

// Inflow boundary condition
double inflow_function(Vector &x);

int main(int argc, char *argv[]) {

  int myRank = 0;
  int numProcs = 1;
   int provided;
   std::cout.precision(8);

#if defined(BROOM_USE_MPI)
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
#endif

#if defined(__linux__)
  feenableexcept(FE_DIVBYZERO);
  feenableexcept(FE_INVALID);
  feenableexcept(FE_OVERFLOW);
#endif

  /// CONSTRUCTOR
   int ref_levels, SnOrder, feOrder, meshOrder, scatIter, quadType;
   double meshDistortionAlpha, tolerance, initialGuess;
   bool visit, directSolve, DSAOption, restart, timer;
   std::string meshFileStr, meshDistortionStr;
   int BuffSize(100);
   char meshFile[BuffSize], meshDistortion[BuffSize], quadFile[BuffSize];
   
   if(myRank == 0) {
      // 1. Parse command-line options.
      ref_levels = tryOption("-r", "refinement levels", argc, argv, 0);
      SnOrder = tryOption("-s", "Sn quadrature order", argc, argv, 4);
      feOrder = tryOption("-e", "Finite element order", argc, argv, 1);
      meshOrder = tryOption("-g", "Geometry (mesh) order", argc, argv, -1);
      meshDistortionAlpha = tryOption("-a", "mesh distortion magnitude", argc, argv, -1.0);
      visit = tryOption("--visit", "--no-visit", "Dump vis files", argc, argv, true);
      directSolve = tryOption("--directSolve", "--no-directSolve", "Do a direct sparse solve using UMFPACK", argc, argv, true);
      DSAOption = tryOption("--DSA", "--no-DSA", "Use the diffusion synthetic acceleration", argc, argv, true);
      meshFileStr = tryOption("-m", "Mesh file", argc, argv, "../../../broom/Meshes/inline-quad.mesh");
      // conver to char to broadcast
      sprintf(meshFile, meshFileStr.c_str());
      meshDistortionStr = tryOption("-d", "Mesh distortion", argc, argv, "none");
      // conver to char to broadcast
      sprintf(meshDistortion, meshDistortionStr.c_str());
      tolerance = tryOption("-t", "Source iteration tolerance", argc, argv, 1.0e-10);
      initialGuess = tryOption("-p", "Initial guess for phi", argc, argv, 0.3);
      probEpsilon = tryOption("-probEpsilon", "factor of epsilon in diffusion limit", argc, argv, 0.1);
      scatIter = tryOption("-i", "Max scattering source iterations", argc, argv, 100);
      restart = tryOption( "--restart","--no-restart", "Use restart file", argc, argv, false);
      quadType = tryOption("--quadType", "Quadrature order type [1=LS from Cheuck, 2=LS from Ardra]", argc, argv, 2);
      timer =     tryOption("--timer", "--no-timer", "Prints timing statistics", argc, argv, true);
      if( quadType == 1) {
         sprintf(quadFile, "../../../broom/src/LS Quadratures/LS_%d.txt", SnOrder / 2);
      } else {
         sprintf(quadFile, "../../../broom/src/LS Quadratures2/S%d_LS_from_ARDRA.txt", SnOrder);
      }
   } // end myRank==0
   
  // TODO The initialGuess actually gets doubled in the first iteration.
  // My first thought is it gets doubled by the weights (since this is a 1D
  // problem; I assume it would be 4x in 2D)
  // initialGuess = 0.5 * initialGuess;
  // TODO check the ^^ previous ^^ TODO
   
#if defined(BROOM_USE_MPI)
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Bcast(&ref_levels, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&SnOrder, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&feOrder, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&meshOrder, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&meshDistortionAlpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(&visit, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&directSolve, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&meshFile, BuffSize, MPI_CHAR, 0, MPI_COMM_WORLD);
   MPI_Bcast(&meshDistortion, BuffSize, MPI_CHAR, 0, MPI_COMM_WORLD);
   MPI_Bcast(&tolerance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(&initialGuess, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(&probEpsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(&scatIter, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&restart, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&quadType, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(quadFile, BuffSize, MPI_CHAR, 0, MPI_COMM_WORLD);
   MPI_Bcast(&timer, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
   
   // Timer initialization
   double t0(0), t1(0), t2(0), t3(0), t4(0), t5(0), t6(0), t7(0);
   double transportTime(0), DSAtime(0), saveTime(0);
   if(timer)
   {
      t0 = MPI_Wtime();
   }

  Quadrature quadrature(quadFile, SnOrder);

  const int numAngles = quadrature.getNumAngles();
  //cout << "numAngles: " << numAngles << "\n";

  // Hack to copy data into global vairables so that mfem function's can access
  // them.
  mu = quadrature.getMu();
  eta = quadrature.getEta();
  w = quadrature.getWeight();

  MeshLoader meshLoader(meshFile, meshOrder, feOrder, ref_levels,
                        meshDistortion, meshDistortionAlpha);

  Mesh &mesh(meshLoader.getMesh());
  if (myRank == 0) {
    mesh.PrintCharacteristics();
  } 
  int dim = mesh.Dimension();

  // 4. Define the discontinuous DG finite element space of the given
  //    polynomial feOrder on the refined mesh.
  DG_FECollection fec(feOrder, dim);
  FiniteElementSpace fes(&mesh, &fec);

  if (myRank == 0) {
    std::cout << "Number of unknowns: " << fes.GetVSize() << '\n';
    std::cout << "Number of zones: " << mesh.GetNE() << '\n';
  }

  GridFunction gphi(&fes);
  GridFunction phi_old(&fes);
  GridFunction diffPhi(&fes);
  vector<double> sigma_s(gphi.Size());
  GridFunction sigma_t_vis(&fes);
  GridFunction sigma_s_vis(&fes);
  GridFunction scat_src(&fes);
  GridFunction scat_src_vis(&fes);
  GridFunction qext_vis(&fes);
  GridFunction psi0_vis(&fes);
  
  
   phi = 0.0;
         
#if defined(BROOM_USE_MPI)
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(gphi, gphi.Size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(phi_old, phi_old.Size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
  
  if(timer)
  {
      t1 = MPI_Wtime();
  }

   // Diffusion Equation Solve
   // TODO only one process needs to do this
   diffPhi = 0;
  
   diffSolve diffSolve(dim, sigma_a_function, diff_Q_function, diff_function, fes, diffPhi, directSolve, myRank);
   
   if(timer)
   {
      t2 = MPI_Wtime();
   }
   
   mfem::GridFunctionCoefficient diffPhiCoef(&diffPhi);
  std::cout << "Diffusion solve complete\n";

  // Define the BlockStructure of the problem, i.e. define the array of
  //    offsets for each variable. The last component of the Array is the sum
  //    of the dimensions of each block.
  Array<int> block_offsets(SnOrder * (SnOrder + 2) +
                           1); // number of variables + 1
  block_offsets[0] = 0;
  for (int i = 1; i < SnOrder * (SnOrder + 2) + 1; i++) {
    block_offsets[i] = fes.GetVSize();
    //cout << block_offsets[i] << "\n";
  }
  block_offsets.PartialSum();
   //block_offsets.Print(cout,1);

  BlockMatrix hugeMatrix(block_offsets);
  BlockVector psi(block_offsets);
  BlockVector rhs(block_offsets);
  rhs=0.0;
  
#if 0
  std::cout << "***********************************************************\n";
  std::cout << "dim(fes) = " << block_offsets[1] - block_offsets[0] << "\n";
  std::cout << "dim(last) = " << block_offsets.Last() << "\n";
  std::cout << "***********************************************************\n";
#endif
   mfem::BilinearForm *L = new BilinearForm(&fes);
   mfem::LinearForm *b(new LinearForm(&fes));
   mfem::SparseMatrix &A = L->SpMat();

   mfem::SparseMatrix* matrixList[SnOrder*(SnOrder+2)];
   mfem::LinearForm* vectorList[SnOrder*(SnOrder+2)];
   mfem::SparseMatrix* weightList[SnOrder*(SnOrder+2)];
   
   for (int i=0;i<SnOrder*(SnOrder+2);i++) {
     //matrixList[i] = 0;
     //matrixList[i]->Print();
   }

  /// SOLVE function ( all the physics functions)
  // LOOP OVER OMEGA
  if(timer){
      t3 = MPI_Wtime();
   }
  
  for (int qloop = 0; qloop < numAngles / numProcs; qloop++) {
  
   mfem::BilinearForm *W(new BilinearForm(&fes));
  
    q = myRank * (numAngles / numProcs) + qloop;
    
    //std::cout << mu[q] << " " << eta[q] << " " << w[q] << "\n";

  mfem::VectorFunctionCoefficient omega(dim, omega_function);
  mfem::FunctionCoefficient inflow(inflow_function);
  mfem::FunctionCoefficient psi0(psi0_function);
  mfem::FunctionCoefficient sigma_t(sigma_t_function);
  mfem::FunctionCoefficient sigma_s(sigma_s_function);
    
 int elementGeometry = fes.GetMesh()->GetElement(0)->GetGeometryType();
  int borderElementGeometry =
      fes.GetMesh()->GetBdrElement(0)->GetGeometryType();

  meshOrder = fes.GetMesh()->GetElementTransformation(0)->Order();
  feOrder = fes.GetOrder(0);
  int elementIntegrationOrder = 0;
  
  if (feOrder == 0) {
    elementIntegrationOrder = 2 * meshOrder;
  } else if (meshOrder == 0) {
    elementIntegrationOrder = 2 * feOrder;
  } else {
    elementIntegrationOrder = 2 * (feOrder)+dim * meshOrder - 1;
  }
  int borderElementIntegrationOrder = elementIntegrationOrder;
  
  const mfem::IntegrationRule &elIntRule =
      mfem::IntRules.Get(elementGeometry, elementIntegrationOrder);
  const mfem::IntegrationRule &bdrIntRule =
      mfem::IntRules.Get(borderElementGeometry, borderElementIntegrationOrder);
      
    mfem::BilinearForm *L(new BilinearForm(&fes));
    mfem::LinearForm *b(new LinearForm(&fes));

    // Add the removal term from the zones.
    mfem::BilinearFormIntegrator *mass = new mfem::MassIntegrator(sigma_t);
    mass->SetIntRule(&elIntRule);
    L->AddDomainIntegrator(mass);

    // Add convection term within the zones
    mfem::BilinearFormIntegrator *convection =
        new mfem::ConvectionIntegrator(omega, +1.0);
    convection->SetIntRule(&elIntRule);
    L->AddDomainIntegrator(convection);

    // Upwind fluxes between zones on each face interior to mesh
    mfem::BilinearFormIntegrator *traceInterior =
        new mfem::DGTraceIntegrator(omega, -1.0, +0.5);
    traceInterior->SetIntRule(&bdrIntRule);

    mfem::BilinearFormIntegrator *transInterior =
        new mfem::TransposeIntegrator(traceInterior);
    L->AddInteriorFaceIntegrator(transInterior);

    // Now add the inflow boundary conditions to the matrix and RHSs
    mfem::BilinearFormIntegrator *traceBdr =
        new mfem::DGTraceIntegrator(omega, -1.0, +0.5);
    traceBdr->SetIntRule(&bdrIntRule);
    mfem::BilinearFormIntegrator *transBdr =
        new mfem::TransposeIntegrator(traceBdr);
    L->AddBdrFaceIntegrator(transBdr);
    
    // Subtract the scattering term of this angle
   mfem::BilinearFormIntegrator *scatter = new mfem::MassIntegrator(sigma_s);
   scatter->SetIntRule(&elIntRule);
   L->AddDomainIntegrator(scatter);
   
   // Add scattering term to weightList for off-diagonals of hugeMatrix
   W->AddDomainIntegrator(scatter);
   W->Assemble();
   W->Finalize();
   
   mfem::SparseMatrix *myW = W->LoseMat();
   mfem::SparseMatrix &spW = *myW;
   weightList[qloop] = &spW;
   //weightList[qloop]->PrintMatlab();

    mfem::LinearFormIntegrator *BndryFlow =
        new mfem::BoundaryFlowIntegrator(inflow, omega, -1.0, -0.5);
    BndryFlow->SetIntRule(&bdrIntRule);
    b->AddBdrFaceIntegrator(BndryFlow);

    mfem::FunctionCoefficient qext(qext_function);
    qext_vis.ProjectCoefficient(qext);

    mfem::LinearFormIntegrator *LFinteriorQext =
        new mfem::DomainLFIntegrator(qext);
    LFinteriorQext->SetIntRule(&elIntRule);
    b->AddDomainIntegrator(LFinteriorQext);

    // Assemble and build the matrix and RHS
    L->Assemble();
    L->Finalize();

    b->Assemble();

    //b->Print(cout,1);
    for (int i = 0; i < fes.GetVSize(); i++) {
      //rhs[qloop * fes.GetVSize() + i] = *b[i];
    }
    vectorList[qloop] = b;
    //b->Print(cout,1);
    //b->Update(&fes,rhs.GetBlock(qloop),0);
    //rhs.Update(b->GetData(),block_offsets);
    //rhs.GetBlock(qloop).Print(cout,1);

    // cout << rhs.Size() << "\n";
    // cout << "rhs blocksize " << qloop << "  " << rhs.BlockSize(qloop) <<
    // "\n";

    //rhs.Print();

    // 6. Define the initial conditions
    // psi.ProjectCoefficient(psi0);
    psi0_vis.ProjectCoefficient(psi0);
    sigma_t_vis.ProjectCoefficient(sigma_t);
    sigma_s_vis.ProjectCoefficient(sigma_s);

    // Get the actual matrix
    mfem::SparseMatrix *myA = L->LoseMat();
    //cout << L << "\n";
    mfem::SparseMatrix &A = *myA;
    //A = L.LoseMat();
    
    //A.LoseData();
    matrixList[qloop] = &A;
    //cout << &matrixList[qloop] << "\n";

    // fill up hugeMatrix with A matrices on diagonal
    //hugeMatrix.SetBlock(qloop, qloop, matrixList[qloop]);
    hugeMatrix.SetBlock(qloop,qloop,&A);
    
    // Print individual block matrices: WORKS individually...
    // hugeMatrix.GetBlock(qloop,qloop).PrintMatlab();

    // cout << "test 5\n";

    // mfem::SparseMatrix S = hugeMatrix;

    //if (qloop == SnOrder * (SnOrder + 2) - 1) {
      //hugeMatrix.PrintMatlab();
      //cout << qloop << "\n";
      //hugeMatrix.GetBlock(qloop,qloop).PrintMatlab();
    //}

  } // END LOOP OVER OMEGA
  
      Vector rhsVec(SnOrder*(SnOrder+2)*fes.GetVSize());
      double *intermediateRHS;
      
   for (int i = 0; i < SnOrder * (SnOrder + 2); i++) {
      //weightList[i]->PrintMatlab();
      for(int j = 0; j < SnOrder * (SnOrder + 2); j++) {
         if(i != j) {
            //cout << i << "  " << j <<"\n";
            hugeMatrix.SetBlock(i,j,weightList[j]); // works.
            //hugeMatrix.GetBlock(i,j).PrintMatlab(); // works.
         }
      }
      //cout << i << "\n";
      //cout << &matrixList[i] << "\n";
      hugeMatrix.SetBlock(i,i,matrixList[i]); // works.
      //hugeMatrix.GetBlock(i,i).PrintMatlab(); // works.
   }
      //cout << "\nblock " << i << "\n";
      //vectorList[i]->Print(cout,1);
      for( int j=0;j<SnOrder*(SnOrder+2);j++) {
            //intermediateRhs.StealData(vectorList[j]);
            intermediateRHS = vectorList[j]->GetData();
         for(int i=0;i<fes.GetVSize();i++) {
            rhsVec[j*fes.GetVSize()+i] = intermediateRHS[i];
            
         }
      }
    ofstream hugeMatrixFile;
    hugeMatrixFile.open ("hugeMatrix.txt");
    hugeMatrix.PrintMatlab(hugeMatrixFile);
    hugeMatrixFile.close();
      
      //hugeMatrix.GetBlock(23,23).PrintMatlab();
      
      //rhsVec.AddElementVector(block_offsets,vectorList[i]->GetData());
      //rhsVec.SetVector(vectorList[i]->GetData(),i*16);
      
      
   
      //rhs.Update(vectorList->GetData(),block_offsets);
      //cout << "\nrhs\n";
      //rhsVec.Print(cout,1);
  //cout << "test out\n";

   //hugeMatrix.GetBlock(0,0).PrintMatlab();
   //rhs.GetBlock(0).Print(cout,1);

  //cout << "This is hugeMatrix after the loop: \n";
  //hugeMatrix.PrintMatlab();

// only solve after the last quadrature angle


// 7. Solve the system
#if defined(MFEM_USE_SUITESPARSE)
  if (!directSolve) {
    cout << "aGMRES solve\n";
#endif
    mfem::DSmoother M(
        A, 1); // l1-Jacobi, we can try l1-GS and polynomial in parallel
    int maxit = 1000;
    double tol = 1e-12;
    double atol = 0.0;
    aGMRES(hugeMatrix, psi, rhsVec, M, maxit, 20, 1, 5, 0.4, tol, atol, 1);
#if defined(MFEM_USE_SUITESPARSE)
  } else {
    cout << "\nMINRESSolve\n";
    int maxIter(1000000000);
    double rtol(1.e-10);
    double atol(1.e-12);
    //MINRESSolver solver;
    GMRESSolver solver;
    solver.SetAbsTol(atol);
    solver.SetRelTol(rtol);
    solver.SetMaxIter(maxIter);
    solver.SetOperator(hugeMatrix);
    solver.SetPrintLevel(1);
    psi = 0.0;

    //cout << rhsVec.Size() << "\n";
    //cout << psi.Size() << "\n";
    //cout << rhsVec.Size() << "\n";
    //rhsVec.Print(cout,1);
    
    //for (int i = 0; i < SnOrder * (SnOrder + 2) + 1; i++) {
    //hugeMatrix.PrintMatlab();
    //}
    //hugeMatrix.GetBlock(0,0).PrintMatlab();

    solver.Mult(rhsVec, psi);

#if 0
  if (solver.GetConverged())
      std::cout << "MINRES converged in " << solver.GetNumIterations()
                << " iterations with a residual norm of " << solver.GetFinalNorm() << ".\n";
   else
      std::cout << "MINRES did not converge in " << solver.GetNumIterations()
                << " iterations. Residual norm is " << solver.GetFinalNorm() << ".\n";
#endif

#if 0
    mfem::UMFPackSolver umf_solver;
    umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
    umf_solver.SetOperator(hugeMatrix);
    umf_solver.Mult(rhsVec, psi);
#endif
  }
  
#if 0
  // KLU Solver
  mfem::KLUSolver klu_solver;
  klu_solver.SetOperator(hugeMatrix);
  
  klu_symbolic *Symbolic;
  klu_numeric *Numeric;
  klu_common *Common;
  int ok;
  
  ok = klu_tsolve (Symbolic, Numeric, psi.Size(), 1, rhsVec, Common);
#endif
#endif

         if(timer)
         {
            t4 = MPI_Wtime();
            transportTime = t4-t3;
         }

   // Get Scalar Flux
   // cout << "\nAngular Flux\n";
   //psi.Print(cout,1);
   for(int j = 0;j< SnOrder*(SnOrder+2);j++) {
      for(int i = 0;i<gphi.Size();i++) {
         gphi[i] += w[j]*psi[j * gphi.Size() + i];
         //cout << "j: " << j << "  i: " << i << "  w: " << w[j] << "  gphi: " << gphi[i] << "\n";
      }
   }
   // cout << gphi.Size() << "\n";
   // cout << gphi << "\n";
   
   

#if 1
   if(timer) {
      t6 = MPI_Wtime();
   }
  /// SAVE function
  VisItDataCollection visit_dc("DirectSolve", &mesh);
  visit_dc.RegisterField("phi", &gphi);
  visit_dc.RegisterField("diffPhi", &diffPhi);
  //visit_dc.RegisterField("psi", &psi);
  visit_dc.RegisterField("sigma_t", &sigma_t_vis);
  visit_dc.RegisterField("sigma_s", &sigma_s_vis);
  visit_dc.RegisterField("scat_src", &scat_src_vis);
  visit_dc.RegisterField("S_0", &qext_vis);
  visit_dc.RegisterField("psi0", &psi0_vis);
  if (visit) {
    visit_dc.SetCycle(0);
    visit_dc.SetTime(0.0);
    visit_dc.Save();
      if (myRank == 0) {
        cout << "L2 error: " << gphi.ComputeL2Error(diffPhiCoef) << "\n";
        //cout << "L2 error: " << gphi.ComputeL2Error(exSol) << "\n";
      }
  }
#endif

      if(timer)
      {
         t7 = MPI_Wtime();
         saveTime = t7-t6;
      }

      if(timer)
            {
               cout << "\n=============================================\n";
               cout << "======           timing stats          ======\n";
               cout << "=============================================\n";
               cout << "Set up time:               " << std::setprecision(3) << t1-t0 << " s. " << "(" << (t1-t0)/(t7-t0)*100 << "%)" << "\n";
               cout << "diffSolve time:            " << t2-t1 << " s. " << "(" << (t2-t1)/(t7-t0)*100 << "%)" << "\n";
               cout << "transport solve time:  " << transportTime << " s. " << "(" << (transportTime)/(t7-t0)*100 << "%)" << "\n";
               cout << "postprocess/save time: " << saveTime << " s. " << "(" << (saveTime)/(t7-t0)*100 << "%)" << "\n";
               cout << "total time:                " << t7-t0 << " s. " << "(" << (t7-t0)/(t7-t0)*100 << "%)" << "\n";
            }

  return 0;
}

// The direction of the discrete ordinate
void omega_function(const Vector &x, Vector &v) {
  int dim = x.Size();
  if (dim == 1) {
    v(0) = mu[q];
  } else if (dim == 2) {
    v(0) = mu[q];
    v(1) = eta[q];
  } else {
    v(0) = mu[q];
    v(1) = eta[q];
    v(2) = xi[q];
  } 
}

// Initial condition
double psi0_function(Vector &x) { return 0.0; }

// The absorption term
double sigma_t_function(Vector &x) { return 1.0 / probEpsilon; }
double sigma_s_function(Vector &x) { 
   double scat = -w[q]*(sigma_t_function(x) - sigma_a_function(x))/(4*M_PI);
   return scat;
}

double qext_function(Vector &x) {
  //     [q = S_0 / (4*pi)]
  return probEpsilon / (4 * M_PI);
}

// Inflow boundary condition (zero for the problems considered in this
// example)
double inflow_function(Vector &x) {
  // psi_inc = phi / sigma_t = s / (4 * pi * sigma_a)
  // return 1/(4*M_PI);
  return 0;
}

double sigma_a_function(Vector &x) { return 1*probEpsilon; }

double diff_function(Vector &x) { 
   double D = 1/(3*sigma_t_function(x));
   return D;
}

double diff_Q_function(Vector &x) {
  //     [q = S_0]
  return qext_function(x)*4*M_PI;
}
