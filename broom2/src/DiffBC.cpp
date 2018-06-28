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
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <sstream>
#include <cstring>
#include <vector>
#include "BroomConfig.hh"
#include "Quadrature.hpp"
#include "diffSolve.hpp"
#include "MeshLoader.hh"
#include "tryOptions.hh"
#include <cstdio>
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
double sigma_a_function(Vector &x);
double qext_function(Vector &x);
double diff_function(Vector &x);
double diff_Q_function(Vector &x);
double exSol_function(Vector &x);

// Inflow boundary condition
double inflow_function(Vector &x);
double probEpsilon;

int main(int argc, char *argv[]) {

  int myRank = 0;
  int numProcs = 1;

#if defined(BROOM_USE_MPI)
  MPI_Init(&argc, &argv);
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
   bool visit, directSolve, DSAOption, restart;
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
      meshFileStr = tryOption("-m", "Mesh file", argc, argv, "../../../broom/Meshes/oneD.mesh");
      // conver to char to broadcast
      sprintf(meshFile, meshFileStr.c_str());
      meshDistortionStr = tryOption("-d", "Mesh distortion", argc, argv, "none");
      // conver to char to broadcast
      sprintf(meshDistortion, meshDistortionStr.c_str());
      tolerance = tryOption("-t", "Source iteration tolerance", argc, argv, 1.0e-10);
      initialGuess = tryOption("-p", "Initial guess for phi", argc, argv, 0.3);
      probEpsilon = tryOption("-probEpsilon", "factor of epsilon in diffusion limit", argc, argv, 0.1);
      scatIter = tryOption("-i", "Max scattering source iterations", argc, argv, 1000);
      restart = tryOption( "--restart","--no-restart", "Use restart file", argc, argv, false);
      quadType = tryOption("--quadType", "Quadrature order type [1=LS from Cheuck, 2=LS from Ardra]", argc, argv, 2);

      if( quadType == 1) {
         sprintf(quadFile, "../../../broom/src/LS Quadratures/LS_%d.txt", SnOrder / 2);
      } else {
         sprintf(quadFile, "../../../broom/src/LS Quadratures2/S%d_LS_from_ARDRA.txt", SnOrder);
      }
   } // end myRank==0
   
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
#endif

  MeshLoader meshLoader(meshFile, meshOrder, feOrder, ref_levels, meshDistortion, meshDistortionAlpha);

   Mesh &mesh(meshLoader.getMesh());
   if(myRank == 0) {
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

  GridFunction phi_conv(&fes);
  GridFunction lphi(&fes);
  GridFunction gphi(&fes);
  GridFunction diffPhi(&fes);
  vector<double> sigma_s(gphi.Size());
  GridFunction sigma_t_vis(&fes);
  GridFunction sigma_s_vis(&fes);
  GridFunction scat_src_vis(&fes);
  GridFunction qext_vis(&fes);
  GridFunction psi0_vis(&fes);
  FunctionCoefficient exSol(exSol_function);
  GridFunction exSol_vis(&fes);
  exSol_vis.ProjectCoefficient(exSol);

  // Move this into a function at the bottom
  for (int i = 0; i < gphi.Size(); i++) {
    phi_conv[i] = 0;
    lphi[i] = 0;
  }
  
   // Diffusion Equation Solve
   diffPhi = 0;

   diffSolve diffSolve(dim, sigma_a_function, diff_Q_function, diff_function, fes, diffPhi, directSolve, myRank);

    /// SAVE function
    VisItDataCollection visit_dc("DiffBC", &mesh);
    visit_dc.RegisterField("diffPhi", &diffPhi);
    visit_dc.RegisterField("exSol", &exSol_vis);
    // visit_dc.RegisterField("psi", &psi);
    //visit_dc.RegisterField("sigma_t", &sigma_t_vis);
    //visit_dc.RegisterField("sigma_s", &sigma_s_vis);
    //visit_dc.RegisterField("scat_src", &scat_src_vis);
    //visit_dc.RegisterField("qext", &qext_vis);
    //visit_dc.RegisterField("psi0", &psi0_vis);
    if (visit) {
      visit_dc.SetCycle(0);
      visit_dc.SetTime(0);
      visit_dc.Save();
      if (myRank == 0) {
        cout << "L2 error: " << diffPhi.ComputeL2Error(exSol) << "\n";
      }
    }

#if defined(BROOM_USE_MPI)
  MPI_Finalize();
#endif

  return 0;
}

// Exact Solution for 1-D DiffBC
double exSol_function(Vector &x){
   double exSol(0);
   double sigma_a(sigma_a_function(x)/1);
   double D(1/(3*sigma_t_function(x)*1));
   double L(sqrt(D / sigma_a));
   double Q(diff_Q_function(x)/1);
  
   double c1(1/4.*Q/sigma_a*((1/4. + 1/2.*D/L)*exp(1/L) - (1/4. - 1/2.*D/L))/ ((1/4. - 1/2.*D/L)*(1/4. - 1/2.*D/L) - (1/4. + 1/2.*D/L)*(1/4. + 1/2.*D/L)*exp(2/L)));
  
   double c2((-c1*(1/4. + 1/2.*D/L)*exp(1/L) - 1/4.*Q/sigma_a)/(1/4. - 1/2.*D/L)*exp(1/L));
  
   exSol = c1*exp(x[0]/L) + c2*exp(-x[0]/L) + Q/sigma_a;
  return exSol;
}

// The absorption term
double sigma_t_function(Vector &x) {
   return 1. / probEpsilon;
}

double qext_function(Vector &x) {
  //     [q = S_0 / (4*pi)]
  return diff_Q_function(x) / (4 * M_PI);
}

double sigma_a_function(Vector &x) {
   return 1*probEpsilon;
}

double diff_function(Vector &x) { 
   return 1/(3*sigma_t_function(x));
}

double diff_Q_function(Vector &x) {
  //     [q = S_0]
  return 1*probEpsilon;
}
