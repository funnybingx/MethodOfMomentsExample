
#include "mom/MomentsCalculator.hpp"
#include "poly/chebyshev.hpp"
#include "utils/help.hpp"
#include <iostream>
#include <assert.h>
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

MomentsCalculator::MomentsCalculator(unsigned int order, RooRealVar* x, RooRealVar* w) :
    _order(order),
    _xvar(x),
    _wvar(w)
{
  // reserve memory for coefficients in chebyshev basis
  _chebyshevBasisCoefficients = std::vector<double> ( _order, 0.0 );
  _coefficientsRRVs = std::vector<RooRealVar*>( _order, NULL );
}

MomentsCalculator::~MomentsCalculator()
{
  // dealocate the RooRealVal coefficients
  for (unsigned int irrv=0; irrv<_coefficientsRRVs.size(); ++irrv)
    delete _coefficientsRRVs[irrv];
}

vector<double> MomentsCalculator::run(const RooDataSet& dataset)
{
  // runs all steps in order
  this->calculateMoments(dataset);
  this->convertToSimplePoly();
  return _simplePolyBasisCoefficients;
}

void MomentsCalculator::calculateMoments(const RooDataSet& dataset)
{
  //---------------------------------------------------------------------------
  // actually run the method of moments (in the Chebyshev basis) over the data
  // populates the vector of Chebyshev coefficients
  //---------------------------------------------------------------------------
  double entries_norm = dataset.sumEntries(); // takes care of weights if exist
  // loop over all orders (==coefficients, since 1D)
  for (unsigned int i=0; i<_order; ++i)
  {
    // calculate this coefficient
    double this_coefficient = 0.0;
    for (int idata=0; idata<dataset.numEntries(); ++idata)
    { 
      const RooArgSet *ras = dataset.get(idata);
      double x = ras->getRealValue(_xvar->GetName());
      double w = 1.0; // optional weight
      if (_wvar) w = ras->getRealValue(_wvar->GetName());
      double sqrt_den = sqrt(1 - x*x); // funny term for cheby orthog relation
      double orth = poly::modified_heaviside(i); // othogonality value
      this_coefficient += w*( poly::chebyshev(x, i) / sqrt_den * orth );
    }
    // persist normalised coefficients
    _chebyshevBasisCoefficients[i] = this_coefficient / entries_norm;
    utils::help::reportOnLoop(i, _order);
  }
  cout << endl; // clear up from loop reporting
  return;
}

void MomentsCalculator::convertToSimplePoly()
{
  return;
}

RooChebychev* MomentsCalculator::getRooPdf()
{
  // get a RooChebychev pdf with these coefficients
  RooArgList chebyCoefficientsRAL = this->convertToRooArgList("cheby");
  if (_debug) chebyCoefficientsRAL.Print();
  return new RooChebychev("pdf", "pdf", *_xvar, chebyCoefficientsRAL);
}

RooArgList MomentsCalculator::convertToRooArgList(TString which)
{
  // go through vector of coefficients and build a RooArgList of RooRealVars
  RooArgList out;
  assert(_chebyshevBasisCoefficients.size() > 1);

  for (unsigned int icoeff=1; icoeff<_order; ++icoeff)
  {
    // need to rescale everything by the zeroth coefficient to make a RooFit 
    // format coefficient list ... RF counts from 1st order, we count from 0th
    TString n="a"; n += icoeff;
    double coeff_RooFit_counting; 
    if (which == "cheby")
      coeff_RooFit_counting = _chebyshevBasisCoefficients[icoeff]/_chebyshevBasisCoefficients[0];
    else if (which == "simple")
      coeff_RooFit_counting = _simplePolyBasisCoefficients[icoeff]/_simplePolyBasisCoefficients[0];
    else
    {
      cerr << "MomentsCalculator::convertToRooArgList ERROR:"
          " unrecognised option: " << which << endl;
      assert(false);
    }
    _coefficientsRRVs[icoeff] = new RooRealVar(n, n, coeff_RooFit_counting);
    out.add(*_coefficientsRRVs[icoeff]);
  }
  return out;
}
