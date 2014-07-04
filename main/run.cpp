
// ROOT
#include <TCanvas.h>
#include <TLegend.h>

// RooFit
#include <RooPlot.h>
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <RooPolynomial.h>
#include <RooDataSet.h>
#include <RooChebychev.h>
#include <RooFit.h>
#include <RooRandom.h>

// MoM framework
#include "mom/MomentsCalculator.hpp"

namespace rf = RooFit; 

struct config
{
  bool debug; // run in debug mode
  int ndata;  // toy data stats
};

void save(RooPlot*& p, TLegend* l=NULL)
{
  TCanvas canvas;
  p->Draw();
  p->SetTitle("Sam\'s 1D MoM implementation");
  //TLegend leg(0.1,0.7,0.48,0.9);
  ////leg.AddEntry(, "toy data");
  //leg.AddEntry("toypdf", "original pdf", "l");
  //leg.AddEntry("moMPdf", "MoM", "l");
  //leg.Draw();
  canvas.SaveAs("plot.pdf");
  return;
}

int main(int argc, char *argv[])
{
  config c;
  c.debug=true;

  // the 1D variable
  RooRealVar* x = new RooRealVar("x", "x", -1, 1);
  x->setRange("reduced_xrange", -0.5, 0.5);

  // build a toy pdf to get toy data
  RooRealVar c1("c1", "linear order coeff", 0.0);
  RooRealVar c2("c2", "quadratic order coeff", 0.5);
  RooPolynomial toyPdf("toypdf", "Toy pdf polynomial", *x, RooArgList(c1, c2));
  RooRandom::randomGenerator()->SetSeed(0);
  RooDataSet* toyData = toyPdf.generate(RooArgSet(*x), 100000);
  //RooDataSet* reducedToy = (RooDataSet*)toyData->reduce("x > -0.5 && x < 0.5");

  // run the method of moments
  MomentsCalculator mc(6, x);
  if (c.debug) mc.setDebug();
  mc.run(*toyData);
  RooChebychev* momPdf = mc.getRooPdf();
  if (c.debug) momPdf->Print();

  // now plot
  RooPlot* frame = x->frame(50);
  toyData->plotOn(frame);
  toyPdf.plotOn(frame, rf::LineColor(kBlue));
  momPdf->plotOn(frame, rf::LineColor(kRed));
  save(frame);

  delete toyData;
  return 0;
}
