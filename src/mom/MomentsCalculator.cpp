
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
  std::vector<double> temp ( _order, 0.0 );
  _chebyshevBasisCoefficients = temp;
}

MomentsCalculator::~MomentsCalculator()
{
}

std::vector<double> MomentsCalculator::run(const RooDataSet& dataset)
{
  // runs all steps in order
  this->calculateMoments(dataset);
  this->convertToSimplePoly();
  return _simplePolyBasisCoefficients;
}

void MomentsCalculator::calculateMoments(const RooDataSet& dataset)
{
  //
  // actually run the method of moments (in the Chebyshev basis) over the data
  // populates the vector of Chebyshev coefficients
  //
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
  vector<RooRealVar*> vecRRVs ( _order, NULL ); // FIXME MEM LEAK
  RooArgList chebyCoefficients = this->convertToRooArgList(vecRRVs, "cheby");
  if (_debug) chebyCoefficients.Print();
  RooChebychev* out = new RooChebychev("pdf", "pdf", *_xvar, chebyCoefficients);
  //delete here?
  return out;
}

RooArgList MomentsCalculator::convertToRooArgList(vector<RooRealVar*>& v, TString which)
{
  // go through vector of coefficients and build a RooArgList of RooRealVars
  RooArgList out;
  for (unsigned int icoeff=0; icoeff<_order; ++icoeff)
  {
    TString name="a"; name += icoeff;
    if (which == "cheby")
      v[icoeff] = new RooRealVar(name, name, _chebyshevBasisCoefficients[icoeff]);
    else if (which == "simple")
      v[icoeff] = new RooRealVar(name, name, _simplePolyBasisCoefficients[icoeff]);
    else
    {
      cerr << "MomentsCalculator::convertToRooArgList ERROR:"
          " unrecognised option: " << which << endl;
      assert(false);
    }
    out.add(*v[icoeff]);
  }
  return out;
}
