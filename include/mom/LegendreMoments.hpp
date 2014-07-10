
#ifndef LEGENDREMOMENTSCALCULATOR_HPP
#define LEGENDREMOMENTSCALCULATOR_HPP

#include "mom/IMomentsCalculator.hpp"

#include <vector>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooLegendre.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooArgList.h>
#include <TTree.h>

class LegendreMoments : protected IMomentsCalculator
{
  public:
   LegendreMoments(unsigned int order, RooRealVar* x, RooRealVar* w=NULL);
   ~LegendreMoments();
   void run(const RooDataSet& data);
   RooAddPdf* getRooPdf();
   void setDebug(bool b=true) { _debug = b; }
   std::vector<double> getMoments(){ return _legendreBasisCoefficients; }
   std::vector<std::vector<double>> getVariances(){ return _legendreBasisVariances; }
  private:
   std::vector<double> _legendreBasisCoefficients;
   std::vector<std::vector<double>> _legendreBasisVariances;
   std::vector<RooRealVar*> _coefficientsRRVs;
   std::vector<RooLegendre*> _legendreTerms;
   void calculateMoments(const RooDataSet&);   // run the Method of Moments
   void calculateVariances(const RooDataSet&); // run the Method of Moments for variances
   RooArgList convertToRooArgList();
};

#endif // LEGENDREMOMENTSCALCULATOR_HPP
