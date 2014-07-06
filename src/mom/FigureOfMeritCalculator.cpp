
#include <TMath.h>         // for TMath::Prob
#include <RooMsgService.h> // make RooFit quiet down
#include <iostream>
#include "mom/FigureOfMeritCalculator.hpp"

using std::cout;
using std::endl;
using std::cerr;
using std::vector;

FigureOfMeritCalculator::FigureOfMeritCalculator(int nbins, RooDataSet* d, 
                                                 RooRealVar* x, RooRealVar* w) :
    _ds(d),
    _xvar(x),
    _wvar(w),
    _nbins(nbins),
    _chi2(-1.0),
    _pval(-1.0),
    _ndof(-1)
{
  // define binning scheme
  _binEdgesLo  = vector<double> (_nbins, 0.0);
  _binEdgesHi  = vector<double> (_nbins, 0.0);
  for (int ibin=0; ibin<_nbins; ++ibin) {
    _binEdgesLo[ibin] = -1.0 + ibin*(2.0)/_nbins;
    _binEdgesHi[ibin] = -1.0 + (ibin+1)*(2.0)/_nbins;
  }
}

FigureOfMeritCalculator::~FigureOfMeritCalculator()
{
  _binEdgesHi.clear();
  _binEdgesLo.clear();
  _binnedData.clear();
  _binnedSumW2.clear();
}

void FigureOfMeritCalculator::print(RooAbsPdf* pdf)
{
  this->pvalue(pdf); // re-calculates everything
  cout << "chi2: " << _chi2
      << " ndof: " << _ndof 
      << " chi2/ndof: " << _chi2/_ndof
      << " p-value: " << _pval << endl;
  return;
}

double FigureOfMeritCalculator::pvalue(RooAbsPdf* pdf, int extraNdof)
{
  _ndof = _nbins + extraNdof;
  this->chi2(pdf);
  _pval = TMath::Prob(_chi2, _ndof);
  return _pval;
}

double FigureOfMeritCalculator::chi2(RooAbsPdf* pdf)
{
  // make roofit shut up about defining bins and calculating integrals
  RooMsgService::instance().setStreamStatus(1,false);

  // populate the data array
  if ( _binnedData.size() == 0 )
    this->populateBinnedDataArrays();

  // loop over bins, calculate chi2
  _chi2 = 0.0;
  double sumEntries = _ds->sumEntries();
  for (int ibin=0; ibin<_nbins; ++ibin) {
    double fractionalInt = this->integratePdfOverBin(ibin, pdf);
    double predictedOccupancy = fractionalInt * sumEntries;
    double diff = predictedOccupancy - _binnedData[ibin];
    _chi2 += diff*diff / _binnedSumW2[ibin];
  }
  RooMsgService::instance().setStreamStatus(1,true);
  return _chi2; // is also now saved in member variable
}

double FigureOfMeritCalculator::chi2Ndof(RooAbsPdf* pdf, int extraNdof)
{
  _ndof = _nbins + extraNdof;
  return this->chi2(pdf) / _ndof;
}


void FigureOfMeritCalculator::populateBinnedDataArrays()
{
  cout << "FigureOfMeritCalculator::populateBinnedDataArrays INFO: filling" << endl;

  // reserve empty arrays for the data
  _binnedData  = vector<double> (_nbins, 0.0);
  _binnedSumW2 = vector<double> (_nbins, 0.0);

  for (int idata=0; idata<_ds->numEntries(); ++idata) {

    // get datum out of RooDataSet
    const RooArgSet* ras = _ds->get(idata);
    double x = ras->getRealValue(_xvar->GetName());
    double w = 1.0; // optional weight
    if (_wvar) w = ras->getRealValue(_wvar->GetName());

    // fill the correct bin
    for (int ibin=0; ibin<_nbins; ++ibin) {
      if ( x > _binEdgesLo[ibin] && x < _binEdgesHi[ibin] ) {
        _binnedData[ibin] += w;
        _binnedSumW2[ibin] += w*w;
      }
    }
  }
  return;
}

double FigureOfMeritCalculator::integratePdfOverBin(int ibin, RooAbsPdf* pdf)
{
  //---------------------------------------------------------------------------
  // Have to calculate fractional integral in this bin range. The full pdf is 
  // not nessicarily normalised so do fraction by hand.
  //---------------------------------------------------------------------------
  double unNormalisedTotal = pdf->createIntegral(RooArgSet(*_xvar))->getVal();
  TString r="bin_"; r+=ibin;
  _xvar->setRange(r, _binEdgesLo[ibin], _binEdgesHi[ibin]);
  double rangeInt = pdf->createIntegral(RooArgSet(*_xvar), r)->getVal();
  return rangeInt / unNormalisedTotal;
}
