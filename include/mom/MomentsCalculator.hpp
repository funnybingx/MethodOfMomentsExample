
#ifndef MOMENTSCALCULATOR_HPP
#define MOMENTSCALCULATOR_HPP

#include <vector>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooChebychev.h>
#include <RooArgList.h>

class MomentsCalculator
{
  public:
   MomentsCalculator(unsigned int order, RooRealVar* x, RooRealVar* w=NULL);
   ~MomentsCalculator();
   std::vector<double> run(const RooDataSet& data);
   RooChebychev* getRooPdf();
   void setDebug(bool b=true) { _debug = b; }
  private:
   unsigned int _order;
   RooRealVar* _xvar;
   RooRealVar* _wvar; // optional weight, NULL <--> weight = 1
   bool _debug;
   std::vector<double> _chebyshevBasisCoefficients;
   std::vector<double> _simplePolyBasisCoefficients;
   //std::vector<RooRealVar*> _chebyshevBasisRRVs;
   void calculateMoments(const RooDataSet&);    // run the Method of Moments
   void convertToSimplePoly(); // switch the basis to simple polynomials
   RooArgList convertToRooArgList(std::vector<RooRealVar*>& v, TString which);
};

#endif // MOMENTSCALCULATOR_HPP
