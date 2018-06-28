// This file is part of BLAST - a high-order finite element hydrocode

#include <fstream>

#include "mfem.hpp"
using namespace mfem;

using namespace std;

int main(int argc, char *argv[])
{
   if (argc != 3)
   {
      cout << "usage: " << argv[0]
           << " in_mesh_file in_gf_file" << endl;
      return 1;
   }

   Mesh *mesh;

   {
      ifstream imesh (argv[1]);
      if (!imesh)
      {
         cout << "can not open mesh file: " << argv[1] << endl;
         return 2;
      }
      mesh = new Mesh (imesh, 1, 1);
   }

   int dim = mesh -> Dimension();
   cout << "mesh dimension: " << dim << endl;

   GridFunction *gf;
   {
      ifstream igf(argv[2]);
      gf = new GridFunction(mesh, igf);
   }

   //////////////////////////////////
   ofstream         scatter_file("scatter.dat");
   int              sd = 4;
   RefinedGeometry *RefG;
   Vector           values;
   DenseMatrix      pointmat, vec_values;
   const IntegrationRule *ir;

   cout << "Enter subdivision factor : " << flush;
   cin >> sd;
   bool vector_field = (gf->VectorDim() > 1);
   scatter_file.precision(16);
   for (int i = 0; i < mesh->GetNE(); i++)
   {
      if (sd > 0)
      {
         RefG = GlobGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
                                           sd, 1);
         ir = &(RefG->RefPts);
      }
      else
      {
         ir = &(IntRules.Get(mesh->GetElementBaseGeometry(i), -sd));
      }
      if (vector_field)
      {
         gf->GetVectorValues (i, *ir, vec_values, pointmat);
      }
      else
      {
         gf->GetValues(i, *ir, values, pointmat);
      }
      for (int j = 0; j < pointmat.Width(); j++)
      {
         for (int k = 0; k < pointmat.Height(); k++)
         {
            scatter_file << pointmat(k, j) << ' ';
         }
         if (vector_field)
         {
            scatter_file << vec_values(0, j);
            for (int k = 1; k < vec_values.Height(); k++)
            {
               scatter_file << ' ' << vec_values(k, j);
            }
            scatter_file << '\n';
         }
         else
         {
            scatter_file << values(j) << '\n';
         }
      }
   }
   scatter_file.close();
   //////////////////////////////////

   delete gf;
   delete mesh;

   return 0;
}

