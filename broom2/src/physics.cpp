
#include "physics.hpp"
#include "mfem.hpp"

int qloop, q;
std::vector<double> mu, eta, xi, w;
mfem::GridFunction phi;
