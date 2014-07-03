
// ROOT
#include <TCanvas.h>

// RooFit
#include <RooPlot.h>
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooDataSet.h>
#include <RooChebychev.h>
#include <RooFit.h>

// MoM framework
#include "mom/MomentsCalculator.hpp"

namespace rf = RooFit; 

void save(RooPlot*& p)
{
  TCanvas canvas;
  p->Draw();
  canvas.SaveAs("plot.pdf");
  return;
}

int main(int argc, char *argv[])
{
  // the 1D variable
  RooRealVar* x = new RooRealVar("x", "x", -1, 1);

  // build a toy pdf to get toy data
  RooRealVar mean("mean", "mean", 0);
  RooRealVar sigma("sigma", "sigma", 0.1);
  RooGaussian toyPdf("toypdf", "Toy pdf Gaussian", *x, mean, sigma);
  RooDataSet* toyData = toyPdf.generate(RooArgSet(*x), 5000);

  // run the method of moments
  MomentsCalculator mc(2, x);
  mc.setDebug();
  mc.run(*toyData);
  RooChebychev* momPdf = mc.getRooPdf();
  momPdf->Print();

  // now plot
  RooPlot* frame = x->frame(50);
  toyData->plotOn(frame);
  toyPdf.plotOn(frame, rf::LineColor(kBlue));
  momPdf->plotOn(frame, rf::LineColor(kRed));
  save(frame);

  delete toyData;
  return 0;
}
