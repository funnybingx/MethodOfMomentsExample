
#ifndef IMOMENTSCALCULATOR_HPP
#define IMOMENTSCALCULATOR_HPP

#include <vector>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooAbsPdf.h>

class IMomentsCalculator
{
 public:
  IMomentsCalculator(unsigned int order, RooRealVar* x, RooRealVar* w=NULL);
  ~IMomentsCalculator() {}
  virtual void run(const RooDataSet& data) = 0;
  virtual RooAbsPdf* getRooPdf() = 0;
  virtual std::vector<double> getMoments() = 0;
  virtual std::vector<std::vector<double>> getVariances() = 0;
  void setDebug(bool b=true) { _debug = b; }
 protected:
  unsigned int _order;
  RooRealVar* _xvar;
  RooRealVar* _wvar; // optional weight, NULL <--> weight = 1
  bool _debug;
  virtual void calculateMoments(const RooDataSet&)   = 0;
  virtual void calculateVariances(const RooDataSet&) = 0;
};

#endif // IMOMENTSCALCULATOR_HPP
