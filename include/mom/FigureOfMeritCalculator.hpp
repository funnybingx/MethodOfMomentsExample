/******************************************************************************
 * Calculates Sam and Mike's figures of merit
 * Hopefully can be extended to include more FoMs and goodness of fit tests.
 * Only implemented in 1D at the moment.
 *****************************************************************************/
#ifndef FIGUREOFMERITCALCULATOR_HPP
#define FIGUREOFMERITCALCULATOR_HPP

#include <vector>
#include <RooAbsPdf.h>
#include <RooDataSet.h>
#include <RooRealVar.h>

class FigureOfMeritCalculator
{
 public:
  FigureOfMeritCalculator(int nbins, RooDataSet* ds, RooRealVar* x, RooRealVar* w=NULL);
  ~FigureOfMeritCalculator();
  void print(RooAbsPdf* pdf);
  double pvalue(RooAbsPdf* pdf, int extraNdof=0);
  double chi2(RooAbsPdf* pdf);
  double chi2Ndof(RooAbsPdf* pdf, int extraNdof=0);
 private:
  RooDataSet* _ds;
  RooRealVar* _xvar;
  RooRealVar* _wvar;
  int _nbins;
  std::vector<double> _binnedData;
  std::vector<double> _binnedSumW2;
  std::vector<double> _binEdgesLo;
  std::vector<double> _binEdgesHi;
  double _chi2;
  double _pval;
  int _ndof;
  void populateBinnedDataArrays();
  double integratePdfOverBin(int, RooAbsPdf*);
};

#endif // FIGUREOFMERITCALCULATOR_HPP
