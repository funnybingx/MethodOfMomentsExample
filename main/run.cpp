
// RooFit
#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooDataSet.h>
#include <RooArgSet.h>

// MoM framework
#include "mom/MomentsCalculator.hpp"

int main(int argc, char *argv[])
{
  // the 1D variable
  RooRealVar* x = new RooRealVar("x", "x", -1, 1);

  // build a toy pdf to get toy data
  RooRealVar mean("mean", "mean", 0);
  RooRealVar sigma("sigma", "sigma", 0.1);
  RooGaussian toypdf("toypdf", "Toy pdf Gaussian", *x, mean, sigma);
  RooDataSet* toydata = toypdf.generate(RooArgSet(*x), 5000);

  // run the method of moments
  MomentsCalculator mc(2, x);
  mc.run(*toydata);

  delete toydata;
  return 0;
}
