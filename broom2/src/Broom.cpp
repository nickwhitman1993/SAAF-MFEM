//                                Broom driver
//
// Description:  This code solves the steady state linear advection
//               equation omega.grad(u) +sigma_t psi = 0, where omega is the
//               direction of photon travel and sigma_t is the total removal
//               term.
//
//               The example demonstrates the use of Discontinuous Galerkin (DG)
//               bilinear forms in MFEM (face integrators) with a linear solve
//               for an advection-like equation.
//
//               The problem is called from ProblemSpecifications.cpp
//

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include "BroomConfig.hh"
#include "Quadrature.hpp"
#include "physics.hpp"
#include "Solve.hpp"
#include "DSASolve.hpp"
#include "diffSolve.hpp"
#include "MeshLoader.hh"
#include "tryOptions.hh"
#include "ProblemSpecifications.hpp"
#include "GridFuncBroom.hpp"
#include <stdio.h>
#include <cstdio>
#include <cmath>
#if defined(BROOM_USE_MPI)
#include <mpi.h>
#endif

#if defined(__linux__)
#include <fenv.h>
#endif

using namespace std;
using namespace mfem;

// The discrete ordinates directions.
void omega_function(const Vector &x, Vector &v);

// Initial condition
double psi0_function(Vector &x);
double sigma_t_function(Vector &x);
double sigma_s_function(Vector &x);
double sigma_a_function(Vector &x);
double diff_function(Vector &x);
double diff_Q_function(Vector &x);
double qext_function(Vector &x);
double exSol_function(Vector &x);

// Inflow boundary condition
double inflow_function(Vector &x);
double probEpsilon;

///////////////////////////////////////////////////////////////////
//    This hack makes it a global variable. Need it for *_function.
// This instance of the ProblemSpecification has all the details
// needed to run a particular problem, like material properties,
// boundary conditions, etc.  We use an auto_ptr so we don't have
// to remember to delete it.
std::unique_ptr<ProblemSpecification> problemSpec;


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
   int refLevels, SnOrder, feOrder, meshOrder, scatIter, quadType;
   double meshDistortionAlpha, tolerance, initialGuess;
   bool visit, directSolve, DSAOption, restart, timer;
   std::string meshFileStr, meshDistortionStr, getProblemStr;
   int BuffSize(100);
   char meshFile[BuffSize], meshDistortion[BuffSize], quadFile[BuffSize], getProblem[BuffSize];

   if(myRank == 0) {
      // 1. Parse command-line options.
      refLevels = tryOption("-r", "refinement levels", argc, argv, 0);
      feOrder = tryOption("-e", "Finite element order", argc, argv, 1);
      SnOrder = tryOption("-s", "Sn quadrature order", argc, argv, 4);
      // TODO choices should now be LS or QR
      quadType = tryOption("--quadType", "Quadrature order type [1=LS from Cheuck, 2=LS from Ardra]", argc, argv, 2);
      timer =     tryOption("--timer", "--no-timer", "Prints timing statistics", argc, argv, false);
      if( quadType == 1) {
         sprintf(quadFile, "../../../broom2/src/LS Quadratures/LS_%d.txt", SnOrder / 2);
      } else {
         sprintf(quadFile, "../../../broom2/src/LS Quadratures2/S%d_LS_from_ARDRA.txt", SnOrder);
      meshOrder = tryOption("-g", "Geometry (mesh) order", argc, argv, -1);
      meshDistortionAlpha = tryOption("-a", "mesh distortion magnitude", argc, argv, -1.0);
      meshFileStr = tryOption("-m", "Mesh file", argc, argv, "../../../broom2/Meshes/inline-quad.mesh");
      // conver to char to broadcast
      sprintf(meshFile, meshFileStr.c_str());
      meshDistortionStr = tryOption("-d", "Mesh distortion", argc, argv, "none");
      // conver to char to broadcast
      sprintf(meshDistortion, meshDistortionStr.c_str());
      visit = tryOption("--visit", "--no-visit", "Dump vis files", argc, argv, true);
      directSolve = tryOption("--directSolve", "--no-directSolve", "Do a direct sparse solve using UMFPACK", argc, argv, true);
      DSAOption = tryOption("--DSA", "--no-DSA", "Use the diffusion synthetic acceleration", argc, argv, true);
      tolerance = tryOption("-t", "Source iteration tolerance", argc, argv, 1.0e-12);
      // TODO remove initial guess - all should be random now. Maybe allow to use
      //    analytic solution for initial guess for debugging.
      initialGuess = tryOption("-p", "Initial guess for phi", argc, argv, 0.3);
      scatIter = tryOption("-i", "Max scattering source iterations", argc, argv, 500);
      getProblemStr = tryOption("-b", "Problem Specification", argc, argv, "MMSGleicher");
      sprintf(getProblem, getProblemStr.c_str());
      probEpsilon = tryOption("-probEpsilon", "factor of epsilon in diffusion limit", argc, argv, 0.1);
      restart = tryOption( "--restart","--no-restart", "Use restart file", argc, argv, false);
      }
   } // end myRank==0



  // TODO The initialGuess actually gets doubled in the first iteration.
  // My first thought is it gets doubled by the weights (since this is a 1D
  // problem; I assume it would be 4x in 2D)
  // initialGuess = 0.5 * initialGuess;
  // TODO check the ^^ previous ^^ TODO

#if defined(BROOM_USE_MPI)
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Bcast(&refLevels, 1, MPI_INT, 0, MPI_COMM_WORLD);
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
   MPI_Bcast(&getProblem, BuffSize, MPI_CHAR, 0, MPI_COMM_WORLD);
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

  // Hack to copy data into global vairables so that mfem function's can access
  // them.
  mu = quadrature.getMu();
  eta = quadrature.getEta();
  xi = quadrature.getXi();
  w = quadrature.getWeight();

  // Put all angle data into a vetor to be passed to ProblemSpecification
  std::vector<std::vector<double> > angleData(4);
  angleData[0] = quadrature.getMu();
  angleData[1] = quadrature.getEta();
  angleData[2] = quadrature.getXi();
  angleData[3] = quadrature.getWeight();

   MeshLoader meshLoader(meshFile, meshOrder, feOrder, refLevels,
                         meshDistortion, meshDistortionAlpha);
   Mesh &mesh(meshLoader.getMesh());
   if (myRank == 0) {
     mesh.PrintCharacteristics();
   }
   int dim = mesh.Dimension();

   // Get vector direction
   std::vector<double> vecDirection(2);
   vecDirection[0] = mu[q];
   vecDirection[1] = eta[q];

   // Parse the input and figure out what problem we're running.
   //setProblemSpecification(getProblem, probEpsilon, vecDirection, problemSpec);
   setProblemSpecification(getProblem, probEpsilon, angleData, problemSpec);
   mfem::Vector test(1);
   if (myRank == 0) {
      //std::cout << problemSpec->computeS0(test) << "\n";
      //std::cout << "computeSigmaT:\t" << problemSpec->computeSigmaT(test) << "\n";
   }



  // 4. Define the discontinuous DG finite element space of the given
  //    polynomial feOrder on the refined mesh.
  enum BasisType{GaussLegendre=0, GaussLobatto=1, Positive=2};
  //L2_FECollection fec(feOrder, dim);
  DG_FECollection fec(feOrder, dim, GaussLegendre);
  //LinearFECollection fec;
  //GaussLinearDiscont2DFECollection fec;
  FiniteElementSpace fes(&mesh, &fec);

  if (myRank == 0) {
    std::cout << "Number of unknowns: " << fes.GetVSize() << '\n';
    std::cout << "Number of zones: " << mesh.GetNE() << '\n';
  }

  GridFunction phi_old(&fes);
  GridFunction phi_conv(&fes);
  GridFunction gphiAbs(&fes);
  GridFunction phi_oldAbs(&fes);
  GridFunction lphi(&fes);
  GridFunction gphi(&fes);
  GridFunction diffPhi(&fes);
  GridFunction lDsaPhi(&fes);
  GridFunction gDsaPhi(&fes);
  GridFunction phi_new(&fes);
  double phi_conv_num(1), phi_conv_denom(1e-20), spec_radius(1);
  vector<double> sigma_s(phi_old.Size());
  GridFunction sigma_t_vis(&fes);
  GridFunction sigma_s_vis(&fes);
  GridFunction sigma_a_vis(&fes);
  GridFunction scat_src_vis(&fes);
  GridFunction qext_vis(&fes);
  GridFunction psi0_vis(&fes);
  FunctionCoefficient exSol(exSol_function);
  GridFunction exSol_vis(&fes);
  exSol_vis.ProjectCoefficient(exSol);
  GridFunction inflow_vis(&fes);

   GridFunction psi_vis0(&fes);
   GridFunction psi_vis1(&fes);
   GridFunction psi_vis2(&fes);
   GridFunction psi_vis3(&fes);

  // Initialize memory: otherwise it's random
  for (int i = 0; i < phi_old.Size(); i++) {
    phi_conv[i] = 0;
    lphi[i] = 0;
  }

   // Generate random initial guess
   double randGuess[gphi.Size()];
   //mfem::GridFunction randGuess;
   if(myRank == 0) {
      double low(1.0);
      double high(1.0e3);
      for(int i=0; i<gphi.Size(); i++) {
         randGuess[i] = exp(log(low) + (log(high) - log(low))*rand()/RAND_MAX);
      }
   }

   // read in restart file
   // TODO restart file variable name
   if (restart) {
      std::ifstream readResFile;
      readResFile.open("./Adams2D001.res");
      std::string line("");
      int count = 0;
      std::string phiInitStr;
      if (readResFile.is_open()) {
         for (int i = 0; i < phi_old.Size(); i++) {
            std::getline(readResFile, line);
            std::istringstream iss(line);
            iss >> phiInitStr;
            phi_old[count] = atof(phiInitStr.c_str());
            gphi[count] = atof(phiInitStr.c_str());
            count++;
         }
      readResFile.close();
      } else {
         std::cout << "Unable to open file: \n";
      }
   } else {
      for (int i = 0; i < phi_old.Size(); i++) {
         phi_old[i] = randGuess[i];
         gphi[i] = randGuess[i];
      }
   }

#if defined(BROOM_USE_MPI)
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(gphi, gphi.Size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(phi_old, phi_old.Size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  // OR we could do this in three loops, but less writing.
  // phi_old = initialGuess;
  // phi = initialGuess;
  // phi_conv = 0.0;

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

   // If problemName == WangRagusa. Used to track spectral radius
   double sr1, sr2, sr3;
   sr1=10;
   sr2=10;
   sr3=10;

   /// SOLVE function ( all the physics functions)
   // LOOP OVER PHI UNTIL CONVERGED //
   int converge = 0; // TODO change this to bool
   int cntconv = 0;
   while (converge == 0 && cntconv < scatIter) {
      // LOOP OVER OMEGA
      // TODO MFEM_Assert that SnOrder is divisible by numProcs
      for (int qloop = 0; qloop < numAngles / numProcs; qloop++) {
         if(timer)
         {
            t3 = MPI_Wtime();
         }

         Solve solve(dim, qloop, q, omega_function, inflow_function, psi0_function, sigma_t_function, sigma_s_function, qext_function, fes, phi_old, lphi, myRank, numAngles, numProcs, directSolve, problemSpec);


         sigma_t_vis = solve.getSigt();
         sigma_s_vis = solve.getSigs();
         qext_vis = solve.getQext();
         scat_src_vis = solve.getScatSrc();
         psi0_vis = solve.getPsi0();
         inflow_vis = solve.getInflow();
         if (q == 0) {
            psi_vis0 = solve.getPsi();
         } else if (q == 1) {
            psi_vis1 = solve.getPsi();
         } else if (q == 2) {
            psi_vis2 = solve.getPsi();
         } else if (q == 3) {
            psi_vis3 = solve.getPsi();
         }

      } // END LOOP OVER OMEGA

         if(timer)
         {
            t4 = MPI_Wtime();
            transportTime += t4-t3;
         }

#if defined(BROOM_USE_MPI)
      gphi = 0;
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Allreduce(lphi.GetData(), gphi.GetData(), phi_old.Size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
      gphi = lphi;
#endif

      // DSA Solve
      gDsaPhi = 0;
      lDsaPhi = 0;
      if(timer)
      {
         t5 = MPI_Wtime();
      }
      if(DSAOption) {
         if(myRank == 0) {
            DSASolve DSASolve(dim, sigma_s_function, sigma_a_function, diff_function, fes, phi_old, lDsaPhi, gphi, directSolve, myRank);
            gDsaPhi = lDsaPhi;
            gphi += gDsaPhi;
         }
      }
      if(timer)
      {
         t6 = MPI_Wtime();
         DSAtime += t6-t5;
      }

#if defined(BROOM_USE_MPI)
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(gphi, gphi.Size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(lDsaPhi, lDsaPhi.Size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

#if 0 // switch to debug and not overwrite.
      // Restart File
      ofstream writeResFile("./Adams2D001.res");
      if (writeResFile.is_open()) {
         for (int i = 0; i < phi_old.Size(); i++) {
            writeResFile << std::setprecision(20) << phi_old[i] << "\n";
         }
      writeResFile.close();
      } else {
         std::cout << "Unable to open writeResFile\n";
      }
#endif

      // TODO only need one process to perform this
      // SAVE function
      char* meshName = new char[getProblemStr.length()];
      std::strcpy(meshName, getProblemStr.c_str());
      VisItDataCollection visit_dc(meshName, &mesh);
      visit_dc.RegisterField("phi", &gphi);
      visit_dc.RegisterField("diffPhi", &diffPhi);
      visit_dc.RegisterField("lDsaPhi", &lDsaPhi);
      visit_dc.RegisterField("exSol", &exSol_vis);
      // visit_dc.RegisterField("psi", &psi);
      visit_dc.RegisterField("sigma_t", &sigma_t_vis);
      visit_dc.RegisterField("sigma_s", &sigma_s_vis);
      visit_dc.RegisterField("scat_src", &scat_src_vis);
      visit_dc.RegisterField("qext", &qext_vis);
      visit_dc.RegisterField("inflow", &inflow_vis);
      visit_dc.RegisterField("psi_0", &psi_vis0);
      visit_dc.RegisterField("psi_1", &psi_vis1);
      visit_dc.RegisterField("psi_2", &psi_vis2);
      visit_dc.RegisterField("psi_3", &psi_vis3);
      if (visit)
      {
         visit_dc.SetCycle(cntconv);
         visit_dc.SetTime(1.0 * cntconv);
         // visit_dc.SetCycle(0.0);
         // visit_dc.SetTime(0.0);
         visit_dc.Save();
         // Only calculate L2 error if we haveAnalyticSolution()
         if (problemSpec->haveAnalyticSolution())
         {
            if (myRank == 0)
            {
               if(dim == 1)
               {
                  //cout << "L2 error: " << gphi.ComputeL2Error(exSol) << "\n";
               }
               else
               {
                  GridFunctionBroom gphiBroom(gphi);
                  cout << "L2 error: " << gphiBroom.ComputeL2ErrorBroom2(exSol) << "\n";
               }
            }
         }
      }

      /// POST-SWEEP UPDATE FUNCTION
      swap(gphi, phi_old);
      for (int i = 0; i < gphi.Size(); i++) { // calculate delta phi between iterations
         phi_conv[i] = phi_old[i] - gphi[i];
         if (phi_conv[i] < 0) { // only want positive magnitudes
            phi_conv[i] = -phi_conv[i];
         }
      }

      // Check for Convergence
      double phiMax = phi_conv.Max();
      phi_conv_num = phiMax;
      for(int i=0; i<gphi.Size(); i++) {
         gphiAbs[i] = std::abs(gphi[i]);
         phi_oldAbs[i] = std::abs(phi_old[i]);
      }
      double gphiMax = gphiAbs.Max();
      double phi_oldMax = phi_oldAbs.Max();

      for (int i = 0; i < gphi.Size(); i++) {
         gphi[i] = 0;
         lphi[i] = 0;
      }

#if 1
      // normal spectral radius
      spec_radius = phi_conv_num / phi_conv_denom;

      if(getProblemStr == "WangRagusa")
      {
         // Keep third-to-last spec_radius
         sr3 = sr2;
         sr2 = sr1;
         sr1 = spec_radius;
      }

      if(myRank == 0) {
         cout << "Itr = " << cntconv << ", phi_conv.max =  " << phiMax << ", spectral radius: " << spec_radius << "\n";
      }
#else
      // This spectral radius is for zero solution.
      if(cntconv > 0) {
         if(gphiMax == 0) {
            if(myRank == 0) {
               mfem_error("ERROR: gphiMax=0. Cannot divide by zero.");
            }
         }
         spec_radius = phi_oldMax / gphiMax;
      } else {
         spec_radius = 1.0e20;
      }

      if(myRank == 0) {
         cout << "Itr = " << cntconv << ", phi_conv.max =  " << phi_oldMax << ", spectral radius: " << spec_radius << "\n";
      }
#endif

      if(timer)
      {
         t7 = MPI_Wtime();
         saveTime += t7-t6;
      }

      phi_conv_denom = phi_conv_num;
      if (phi_conv_num < tolerance * (1 - spec_radius) * std::abs(phi_old.Max())) {
         if(myRank == 0) {
            cout << "phi_conv_num:\t" << phi_conv_num << "\n";
         }
         converge = 1;
         if (myRank == 0) {
            cout << "BINGO, converged" << "\n";
            cout << "scattering source iterations: " << cntconv << "\n";
            if(timer)
            {
               cout << "\n=============================================\n";
               cout << "======           timing stats          ======\n";
               cout << "=============================================\n";
               cout << "Set up time:               " << std::setprecision(3) << t1-t0 << " s. " << "(" << (t1-t0)/(t7-t0)*100 << "%)" << "\n";
               cout << "diffSolve time:            " << t2-t1 << " s. " << "(" << (t2-t1)/(t7-t0)*100 << "%)" << "\n";
               cout << "avg transport solve time:  " << transportTime/cntconv << " s. " << "(" << (transportTime/cntconv)/(t7-t0)*100 << "%)" << "\n";
               cout << "avg DSA solve time:        " << DSAtime/cntconv << " s. " << "(" << (DSAtime/cntconv)/(t7-t0)*100 << "%)" << "\n";
               cout << "avg postprocess/save time: " << saveTime/cntconv << " s. " << "(" << (saveTime/cntconv)/(t7-t0)*100 << "%)" << "\n";
               cout << "total time:                " << t7-t0 << " s. " << "(" << (t7-t0)/(t7-t0)*100 << "%)" << "\n";
            }
         }

         if (getProblemStr == "WangRagusa") {
            bool spec_rad_study(true);
            if(myRank == 0) {
               if(spec_rad_study) {
                  std::ostringstream oss;
                  oss << "WR_specRad_" << feOrder << ".out";
                  std::string strSRFile = oss.str();
                  // Use third-to-last spec_radius as
                  std::ofstream specRad;
                  specRad.open(strSRFile.c_str(), std::ofstream::app);
                  specRad << sr3 << "\n";
                  specRad.close();
               }
            }
         }
      }

      cntconv++;

   } // END CONVERGE LOOP

#if defined(BROOM_USE_MPI)
  MPI_Finalize();
#endif

  return 0;
}

// The directions of the discrete ordinates
void omega_function(const Vector &x, Vector &v) {
  int dim = x.Size();
  if (dim == 1) {
    v(0) = mu[q];
  } else if (dim == 2) {
      // MFEM uses these backwards from what is expected.
      //    v(0) is in the y-direction
      //    v(1) is in the x-direction
      // Compare to x[0] is x-direction and x[1] is y-direction.
    v(0) = eta[q];
    v(1) = mu[q];
  } else {
      // Not sure about the ordering of v(i) in 3-D...
    v(0) = mu[q];
    v(1) = eta[q];
    v(2) = xi[q];
  }
}

// Exact solution
double exSol_function(Vector &x){
   return problemSpec->computeAnalyticSolution(x);
}

// Initial condition
double psi0_function(Vector &x){
   return problemSpec->computeInitialCond(x);
}

// Total cross section
double sigma_t_function(Vector &x) {
   //std::cout << "computeSigmaT:\t" << problemSpec->computeSigmaT(x) << "\n";
   return problemSpec->computeSigmaT(x);
}

// Scattering cross section
double sigma_s_function(Vector &x) {
   return problemSpec->computeSigmaS(x);
}

// Isotropic source
double qext_function(Vector &x) {
  return problemSpec->computeS0(x);
}

// Absorption cross section
double sigma_a_function(Vector &x) {
   return problemSpec->computeSigmaA(x);
}

// Diffusion coefficient
double diff_function(Vector &x) {
   return problemSpec->computeD(x);
}

// Diffusion equation isotropic source
double diff_Q_function(Vector &x) {
  return problemSpec->computeDiffQ(x);
}

// Inflow boundary condition
double inflow_function(Vector &x) {
   return problemSpec->computeInflow(x);
}
