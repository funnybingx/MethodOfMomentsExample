
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
   MomentsCalculator(unsigned int order, RooRealVar* x);
   ~MomentsCalculator();
   std::vector<double> run(const RooDataSet& data);
   RooChebychev* getRooPdf();
  private:
   unsigned int _order;
   RooRealVar* _xvar;
   std::vector<double> _chebyshevBasisCoefficients;
   std::vector<double> _simplePolyBasisCoefficients;
   void calculateMoments(const RooDataSet&);    // run the Method of Moments
   void convertToSimplePoly(); // switch the basis to simple polynomials
   RooArgList convertToRooArgList(); // make a RooArgList iiof the llist
};

#endif // MOMENTSCALCULATOR_HPP
