//                                Swiffer driver
//
// Description:  This code solves the steady state linear diffusion
//               equation grad.D grad(u) + sigma_a psi = Q, where D is the
//               diffusion coefficient, sigma_a is the absorption term, and
//               Q is a volumetric source
//
//               The example demonstrates the use of Discontinuous Galerkin (DG)
//               bilinear forms in MFEM (face integrators) with a linear solve
//               for a diffusion equation.
//
//               The problem is called from SwifferProblems.cpp
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
#include "SwifferProblems.hpp"
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
double sigma_t_function(Vector &x);
double sigma_s_function(Vector &x);
double sigma_a_function(Vector &x);
double diff_function(Vector &x);
double diff_Q_function(Vector &x);
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
std::unique_ptr<SwifferProblemSpec> problemSpec;


int main(int argc, char *argv[]) {

   std::cout.precision(8);

#if defined(__linux__)
   feenableexcept(FE_DIVBYZERO);
   feenableexcept(FE_INVALID);
   feenableexcept(FE_OVERFLOW);
#endif

   /// CONSTRUCTOR
   int refLevels, feOrder, meshOrder;
   double meshDistortionAlpha;
   bool visit, directSolve, restart, timer;
   std::string meshFileStr, meshDistortionStr, getProblemStr;
   int BuffSize(100);
   char meshFile[BuffSize], meshDistortion[BuffSize], getProblem[BuffSize];
   
   // 1. Parse command-line options.
   refLevels = tryOption("-r", "refinement levels", argc, argv, 0);
   feOrder = tryOption("-e", "Finite element order", argc, argv, 1);
   timer =     tryOption("--timer", "--no-timer", "Prints timing statistics", argc, argv, false);
   meshOrder = tryOption("-g", "Geometry (mesh) order", argc, argv, -1);
   meshDistortionAlpha = tryOption("-a", "mesh distortion magnitude", argc, argv, -1.0);
   meshFileStr = tryOption("-m", "Mesh file", argc, argv, "../../../broom/Meshes/inline-quad.mesh");
   // conver to char to broadcast
   sprintf(meshFile, meshFileStr.c_str());
   meshDistortionStr = tryOption("-d", "Mesh distortion", argc, argv, "none");
   // conver to char to broadcast
   sprintf(meshDistortion, meshDistortionStr.c_str());
   visit = tryOption("--visit", "--no-visit", "Dump vis files", argc, argv, true);
   directSolve = tryOption("--directSolve", "--no-directSolve", "Do a direct sparse solve using UMFPACK", argc, argv, true);
   getProblemStr = tryOption("-b", "Problem Specification", argc, argv, "uniformInfMed");
   sprintf(getProblem, getProblemStr.c_str());
   probEpsilon = tryOption("-probEpsilon", "factor of epsilon in diffusion limit", argc, argv, 0.1);
   restart = tryOption( "--restart","--no-restart", "Use restart file", argc, argv, false);
   
   // Timer initialization
   double t0(0), t1(0), t2(0), t3(0), t4(0), t5(0), t6(0), t7(0);
   double saveTime(0);
   if(timer)
   {
      t0 = MPI_Wtime();
   }

   MeshLoader meshLoader(meshFile, meshOrder, feOrder, refLevels,
                         meshDistortion, meshDistortionAlpha);
   Mesh &mesh(meshLoader.getMesh());
   mesh.PrintCharacteristics();
   int dim = mesh.Dimension();
   
   // Parse the input and figure out what problem we're running.
   setProblemSpecification(getProblem, probEpsilon, problemSpec);

   // 4. Define the discontinuous DG finite element space of the given
   //    polynomial feOrder on the refined mesh.
   enum BasisType{GaussLegendre=0, GaussLobatto=1, Positive=2};
   //L2_FECollection fec(feOrder, dim);
   DG_FECollection fec(feOrder, dim, GaussLegendre);
   //LinearFECollection fec;
   //GaussLinearDiscont2DFECollection fec;
   FiniteElementSpace fes(&mesh, &fec);

   std::cout << "Number of unknowns: " << fes.GetVSize() << '\n';
   std::cout << "Number of zones: " << mesh.GetNE() << '\n';

   GridFunction phi_old(&fes);
   GridFunction phi_conv(&fes);
   GridFunction gphiAbs(&fes);
   GridFunction phi_oldAbs(&fes);
   GridFunction lphi(&fes);
   GridFunction gphi(&fes);
   GridFunction diffPhi(&fes);
   GridFunction phi_new(&fes);
   vector<double> sigma_s(phi_old.Size());
   GridFunction sigma_t_vis(&fes);
   GridFunction sigma_s_vis(&fes);
   GridFunction sigma_a_vis(&fes);
   FunctionCoefficient exSol(exSol_function);
   GridFunction exSol_vis(&fes);
   exSol_vis.ProjectCoefficient(exSol);
   GridFunction inflow_vis(&fes);

  // Initialize memory: otherwise it's random
  for (int i = 0; i < phi_old.Size(); i++) {
    phi_conv[i] = 0;
    lphi[i] = 0;
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
   }

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
   
   int tmpMyRank = 0;
   
   diffSolve diffSolve(dim, sigma_a_function, diff_Q_function, diff_function, fes, gphi, directSolve, tmpMyRank);
   
   if(timer)
   {
      t2 = MPI_Wtime();
   }
   
   mfem::GridFunctionCoefficient diffPhiCoef(&diffPhi);
         
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
   visit_dc.RegisterField("exSol", &exSol_vis);
   visit_dc.RegisterField("sigma_t", &sigma_t_vis);
   visit_dc.RegisterField("sigma_s", &sigma_s_vis);
   visit_dc.RegisterField("inflow", &inflow_vis);
   if (visit)
   {
      visit_dc.SetCycle(0);
      visit_dc.SetTime(0);
      // visit_dc.SetCycle(0.0);
      // visit_dc.SetTime(0.0);
      visit_dc.Save();
      // Only calculate L2 error if we haveAnalyticSolution()
      if (problemSpec->haveAnalyticSolution())
      {
         if(dim == 1)
         {
            cout << "L2 error: " << gphi.ComputeL2Error(exSol) << "\n";
         }
         else
         {
            GridFunctionBroom gphiBroom(gphi);
            cout << "L2 error: " << gphiBroom.ComputeL2ErrorBroom2(exSol) << "\n";
         }
         
      }
   }

   if(timer)
   {
      t7 = MPI_Wtime();
      saveTime += t7-t6;
   }
   
   cout << "BINGO, solved\n";
   
   if(timer)
   {
      cout << "\n=============================================\n";
      cout << "======           timing stats          ======\n";
      cout << "=============================================\n";
      cout << "Set up time:               " << std::setprecision(3) << t1-t0 << " s. " << "(" << (t1-t0)/(t7-t0)*100 << "%)" << "\n";
      cout << "diffSolve time:            " << t2-t1 << " s. " << "(" << (t2-t1)/(t7-t0)*100 << "%)" << "\n";
      cout << "avg postprocess/save time: " << saveTime << " s. " << "(" << (saveTime)/(t7-t0)*100 << "%)" << "\n";
      cout << "total time:                " << t7-t0 << " s. " << "(" << (t7-t0)/(t7-t0)*100 << "%)" << "\n";
   }

  return 0;
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
