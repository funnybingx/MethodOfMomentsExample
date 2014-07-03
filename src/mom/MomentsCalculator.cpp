
#include "mom/MomentsCalculator.hpp"

MomentsCalculator::MomentsCalculator(unsigned int order, RooRealVar* x) :
    _order(order),
    _xvar(x)
{
  std::vector<double> temp ( _order, 0.0 );
  _chebyshevBasisCoefficients = temp;

  std::vector<double> ttemp ( _order, 0.0 );
  _simplePolyBasisCoefficients = ttemp;
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
  for (unsigned int iterm=0; iterm<_order; ++iterm)
  {
    double coefficient = 0.0;
    for (int idata=0; idata<dataset.numEntries(); ++idata)
    { 
      const RooArgSet *ras = dataset.get(idata);
      double x = ras->getRealValue(_xvar->GetName());
      x+=1;// do MoM here
    }
    _chebyshevBasisCoefficients[iterm] = coefficient;
  }
  return;
}

void MomentsCalculator::convertToSimplePoly()
{
  return;
}

RooChebychev* MomentsCalculator::getRooPdf()
{
  RooArgList chebyCoefficients = this->convertToRooArgList();
  return new RooChebychev("pdf", "pdf", *_xvar, chebyCoefficients);
}

RooArgList MomentsCalculator::convertToRooArgList()
{
  RooArgList out;
  for (unsigned int icoeff=0; icoeff<_order; ++icoeff)
  {
    TString name="a"; name += icoeff;
    out.add(RooRealVar(name, name, _chebyshevBasisCoefficients[icoeff]));
  }
  return out;
}
