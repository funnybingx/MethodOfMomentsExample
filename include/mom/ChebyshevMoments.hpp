
#ifndef MOMENTSCALCULATOR_HPP
#define MOMENTSCALCULATOR_HPP

#include "mom/IMomentsCalculator.hpp"

#include <vector>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooChebychev.h>
#include <RooArgList.h>
#include <TTree.h>

class ChebyshevMoments : protected IMomentsCalculator
{
 public:
  ChebyshevMoments(unsigned int order, RooRealVar* x, RooRealVar* w=NULL);
  ~ChebyshevMoments();
  void run(const RooDataSet& data);
  RooChebychev* getRooPdf();
  void setDebug(bool b=true) { _debug = b; }
  std::vector<double> getMoments(){ return _chebyshevBasisCoefficients; }
  std::vector<std::vector<double>> getVariances(){ return _chebyshevBasisVariances; }
 private:
  std::vector<double> _chebyshevBasisCoefficients;
  std::vector<std::vector<double>> _chebyshevBasisVariances;
  std::vector<RooRealVar*> _coefficientsRRVs;
  void calculateMoments(const RooDataSet&);   // run the Method of Moments
  void calculateVariances(const RooDataSet&); // run the Method of Moments for variances
  RooArgList convertToRooArgList();
};

#endif // MOMENTSCALCULATOR_HPP
