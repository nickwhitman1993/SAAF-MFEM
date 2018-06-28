#ifndef QUADRATURE_HPP
#define QUADRATURE_HPP

#include <vector>
#include <cassert>
#include <string>

class Quadrature {

public:
  explicit Quadrature(const char *quadType, const int order);

  int getOrder() const { return order; }
  int getNumAngles() const { return numAngles; }

  const std::vector<double> &getMu() const { return mu; }
  const std::vector<double> &getEta() const { return eta; }
  const std::vector<double> &getXi() const { return xi; }
  const std::vector<double> &getWeight() const { return w; }

  double getMu(int i) const {
    assert(i < numAngles and i >= 0);
    return mu[i];
  }
  double getEta(int i) const {
    assert(i < numAngles and i >= 0);
    return eta[i];
  }
  double getXi(int i) const {
    assert(i < numAngles and i >= 0);
    return xi[i];
  }
  double getWeight(int i) const {
    assert(i < numAngles and i >= 0);
    return w[i];
  }

private:
  const int order;
  const int numAngles;
  std::vector<double> mu;
  std::vector<double> eta;
  std::vector<double> xi;
  std::vector<double> w;
};

#endif
