
#ifndef LEGENDREMOMENTSCALCULATOR_HPP
#define LEGENDREMOMENTSCALCULATOR_HPP

#include <vector>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooLegendre.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooArgList.h>
#include <TTree.h>

typedef std::vector<std::vector<double> > vvd;
class LegendreMomentsCalculator
{
  public:
   LegendreMomentsCalculator(unsigned int order, RooRealVar* x, RooRealVar* w=NULL);
   ~LegendreMomentsCalculator();
   void run(const RooDataSet& data);
   RooAddPdf* getRooPdf();
   void setDebug(bool b=true) { _debug = b; }
   vvd getVariances(){ return _legendreBasisVariances; }
   std::vector<double> getMoments(){ return _legendreBasisCoefficients; }
  private:
   unsigned int _order;
   RooRealVar* _xvar;
   RooRealVar* _wvar; // optional weight, NULL <--> weight = 1
   bool _debug;
   std::vector<double> _legendreBasisCoefficients;
   vvd _legendreBasisVariances;
   std::vector<RooRealVar*> _coefficientsRRVs;
   std::vector<RooLegendre*> _legendreTerms;
   void calculateMoments(const RooDataSet&);   // run the Method of Moments
   void calculateVariances(const RooDataSet&); // run the Method of Moments for variances
   RooArgList convertToRooArgList();
};

#endif // LEGENDREMOMENTSCALCULATOR_HPP
