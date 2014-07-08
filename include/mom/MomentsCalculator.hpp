
#ifndef MOMENTSCALCULATOR_HPP
#define MOMENTSCALCULATOR_HPP

#include <vector>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooChebychev.h>
#include <RooArgList.h>
#include <TTree.h>

typedef std::vector<std::vector<double> > vvd;
class MomentsCalculator
{
  public:
   MomentsCalculator(unsigned int order, RooRealVar* x, RooRealVar* w=NULL);
   ~MomentsCalculator();
   void run(const RooDataSet& data);
   RooChebychev* getRooPdf();
   void setDebug(bool b=true) { _debug = b; }
   vvd getVariances(){ return _chebyshevBasisVariances; }
  std::vector<double> getMoments(){ return _chebyshevBasisCoefficients; }
  private:
   unsigned int _order;
   RooRealVar* _xvar;
   RooRealVar* _wvar; // optional weight, NULL <--> weight = 1
   bool _debug;
   std::vector<double> _chebyshevBasisCoefficients;
   vvd _chebyshevBasisVariances;
   std::vector<double> _simplePolyBasisCoefficients;
   std::vector<RooRealVar*> _coefficientsRRVs;
   void calculateMoments(const RooDataSet&);   // run the Method of Moments
   void calculateVariances(const RooDataSet&); // run the Method of Moments for variances
   void convertToSimplePoly(); // switch the basis to simple polynomials
   RooArgList convertToRooArgList(TString which="cheby");
};

#endif // MOMENTSCALCULATOR_HPP
