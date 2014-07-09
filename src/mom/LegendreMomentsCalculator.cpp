
#include "mom/LegendreMomentsCalculator.hpp"
#include <boost/math/special_functions/legendre.hpp>
#include "utils/help.hpp"
#include <RooArgSet.h>
#include <iostream>
#include <assert.h>
using std::vector;
using std::cout;
using std::cerr;
using std::endl;

//namespace boost::math 

LegendreMomentsCalculator::LegendreMomentsCalculator(unsigned int order, RooRealVar* x, RooRealVar* w) :
    _order(order),
    _xvar(x),
    _wvar(w)
{
  // reserve memory for coefficients in chebyshev basis
  _legendreBasisCoefficients = std::vector<double> ( _order, 0.0 );
  // reserve memory for variances in chebyshev basis making 
  // use of the initialisation of the basis coeffs for simplicity
  _legendreBasisVariances = vvd( _order, _legendreBasisCoefficients);
  _coefficientsRRVs = std::vector<RooRealVar*>( _order, NULL );
  _legendreTerms = std::vector<RooLegendre*>( _order, NULL );
}

LegendreMomentsCalculator::~LegendreMomentsCalculator()
{
  // dealocate the RooRealVar coefficients and the RooLegendre terms
  for (unsigned int irrv=0; irrv<_coefficientsRRVs.size(); ++irrv)
    delete _coefficientsRRVs[irrv];
  for (unsigned int irlg=0; irlg<_legendreTerms.size(); ++irlg)
    delete _legendreTerms[irlg];
}

void LegendreMomentsCalculator::run(const RooDataSet& dataset)
{
  // runs all steps in order
  this->calculateMoments(dataset);
  this->calculateVariances(dataset);
}

void LegendreMomentsCalculator::calculateMoments(const RooDataSet& dataset)
{
  //---------------------------------------------------------------------------
  // actually run the method of moments (in the Legendre basis) over the data
  // populates the vector of coefficients for Legendre polynomials in a series
  //---------------------------------------------------------------------------
  double entries_norm = dataset.sumEntries(); // takes care of weights if exist
  // loop over all orders (==coefficients, since 1D)
  for (unsigned int i=0; i<_order; ++i)
  {
    // calculate this coefficient and its variance
    double this_coefficient = 0.0;
    for (int idata=0; idata<dataset.numEntries(); ++idata)
    { 
      const RooArgSet *ras = dataset.get(idata);
      double x = ras->getRealValue(_xvar->GetName());
      double w = 1.0; // optional weight
      if (_wvar) w = ras->getRealValue(_wvar->GetName());
      double orthog = 2.0 / (2.0*i + 1.0);
      this_coefficient += w*( boost::math::legendre_p(i, x) ) / orthog;
    }
    // persist normalised coefficients
    _legendreBasisCoefficients[i] = this_coefficient / entries_norm;
    if (_debug) utils::help::reportOnLoop(i, _order);
  }
  if (_debug) cout << endl; // clear up from loop reporting
  return;
}

void LegendreMomentsCalculator::calculateVariances(const RooDataSet& dataset)
{
  //---------------------------------------------------------------------------
  // now calculate the variances i.e. second moments (in the Legendre basis)
  // populates the vector of variances of Legendre coefficients
  //---------------------------------------------------------------------------
  double entries_norm = dataset.sumEntries(); // takes care of weights if exist
  // loop over all orders (==coefficients, since 1D)
  for (unsigned int i=0; i<_order; ++i)
  {
    for (unsigned int j=0; j<_order; ++j)
    {
      // calculate this coefficient and its variance
      double this_variance = 0.0;
      for (int idata=0; idata<dataset.numEntries(); ++idata)
      { 
        const RooArgSet *ras = dataset.get(idata);
        double x = ras->getRealValue(_xvar->GetName());
        double w = 1.0; // optional weight
        if (_wvar) w = ras->getRealValue(_wvar->GetName());
        double term1 = (boost::math::legendre_p(i, x))-_legendreBasisCoefficients[i];
        double term2 = (boost::math::legendre_p(j, x))-_legendreBasisCoefficients[j];
        this_variance += w*( term1*term2 );
      }
      // persist normalised coefficients
      _legendreBasisVariances[i][j] = this_variance / entries_norm;
      if (_debug) utils::help::reportOnLoop(i, _order);
    }
  }
  if (_debug) cout << endl; // clear up from loop reporting
  return;
}

RooAddPdf* LegendreMomentsCalculator::getRooPdf()
{
  // get a RooAddPdf of a series of Legendres with these coefficients
  RooArgList legendreCoefficientsRAL = this->convertToRooArgList();
  RooArgList legendreTermsRAL;
  for (unsigned int icoeff=0; icoeff<_order; ++icoeff) {
    TString n="legendreTerm"; n+=icoeff;
    _legendreTerms[icoeff] = new RooLegendre(n, n, *_xvar, icoeff); 
    _legendreTerms[icoeff]->Print(); // half completes
    legendreTermsRAL.add( *_legendreTerms[icoeff] );
    // a new P_l^m (x) where  l=icoeff (this term's order) m=0
    // RooFit's RooLegendre implementation swaps the order to that of the boost
    // libraries
  }
  if (_debug) legendreCoefficientsRAL.Print();
  if (_debug) legendreTermsRAL.Print();
  return new RooAddPdf("pdf", "pdf", legendreTermsRAL, legendreCoefficientsRAL);
}

RooArgList LegendreMomentsCalculator::convertToRooArgList()
{
  // go through vector of coefficients and build a RooArgList of RooRealVars
  RooArgList out;
  assert(_legendreBasisCoefficients.size() > 1);

  for (unsigned int icoeff=0; icoeff<_order; ++icoeff)
  {
    // need to rescale everything by the zeroth coefficient to make a RooFit 
    // format coefficient list ... RF counts from 1st order, we count from 0th
    TString n="a"; n += icoeff;
    double coeff_RooFit_counting; 
    coeff_RooFit_counting = _legendreBasisCoefficients[icoeff]; ///_legendreBasisCoefficients[0];
    _coefficientsRRVs[icoeff] = new RooRealVar(n, n, coeff_RooFit_counting);
    out.add(*_coefficientsRRVs[icoeff]);
  }
  return out;
}
