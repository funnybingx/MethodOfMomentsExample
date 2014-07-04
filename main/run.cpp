
// 
#include <iostream>
#include "boost/program_options.hpp"

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
using std::cout;
using std::cerr;
using std::endl;

struct config
{
  bool debug;
  int ndata;  
  int ntoys;
  bool fit;
  int ordermom;
};

int parseOptions(config &c, int argc, char *argv[])
{
  // alias for namespace
  namespace po = boost::program_options;

  // declare options
  po::options_description desc("Allowed options") ;
  desc.add_options()
    ("help,h", "show this help")
    ("debug", po::bool_switch(&c.debug)->default_value(true), "run in debug")
    ("ndata", po::value<int>(&c.ndata)->default_value(5000), "toy generated stats")
    ("ntoys", po::value<int>(&c.ntoys)->default_value(1), "number of toys to run")
    ("fit", po::bool_switch(&c.fit), "also run a fit over each toy")
    ("order-mom,o", po::value<int>(&c.ordermom)->default_value(3), "max order of the MoM")
    ;

  // actually do the parsing
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  // show help and exit
  if ((vm.count("help"))) {
    cout << desc << "\n";
    return 1;
  }
  return 0;
}

double chi2FromRooPlot(RooPlot* p, TString rfname1, TString rfname2)
{
  // wrap roofit shite to get the chi2 from a RooPlot. you need to set the names
  // of the components of the RooPlot using the RooFit::Name cmd arg.
  //  i.e. data->plotOn(frame, RooFit::Name("thisname")
  double chi2 = p->chiSquare(rfname1, rfname2);
  cout << "chi2: " << chi2 << endl;
  return chi2;
}

void plot(RooRealVar* x, RooDataSet* d, RooAbsPdf* orig, RooAbsPdf* mom,
          RooAbsPdf* fitted=NULL)
{
  TCanvas canvas;
  TLegend leg(0.1,0.7,0.48,0.9);
  RooPlot* frame = x->frame(50);
  d->plotOn(frame, rf::Name("data"));
  orig->plotOn(frame, rf::LineColor(kBlue), rf::Name("orig"));
  mom->plotOn(frame, rf::LineColor(kRed), rf::Name("mom"));
  chi2FromRooPlot(frame, "mom", "data");
  if (fitted) fitted->plotOn(frame, rf::LineColor(kGreen), rf::Name("fitted"));
  frame->Draw();
  frame->SetTitle("Sam\'s 1D MoM implementation");
  leg.AddEntry("data", "toy data", "p");
  leg.AddEntry("orig", "original pdf", "l");
  leg.AddEntry("mom", "MoM", "l");
  if (fitted) leg.AddEntry("fitted", "RooFit", "l");
  leg.Draw();
  canvas.SaveAs("plot.pdf");
  return;
}

int main(int argc, char *argv[])
{
  // deal with options parsing
  config c;
  int returnCode = parseOptions(c, argc, argv);
  if (returnCode!= 0) { return returnCode; }

  // the 1D variable
  RooRealVar* x = new RooRealVar("x", "x", -1, 1);

  // build a toy pdf to get toy data
  RooRealVar c1("c1", "linear order coeff", 0.0);
  RooRealVar c2("c2", "quadratic order coeff", 0.5);
  RooPolynomial toyPdf("toypdf", "Toy pdf polynomial", *x, RooArgList(c1, c2));
  RooRandom::randomGenerator()->SetSeed(0);
  RooDataSet* toyData = toyPdf.generate(RooArgSet(*x), 100000);

  // run the method of moments
  for (int itoy=0; itoy<c.ntoys; ++itoy) {
    MomentsCalculator moments(c.ordermom, x);
    if (c.debug) moments.setDebug();
    moments.run(*toyData);
    RooChebychev* momPdf = moments.getRooPdf();
    if (c.debug) momPdf->Print();

    // also run a RooFit fit
    RooChebychev* fitPdf=NULL;
    if (c.fit) {
      // implement me!
    }

    plot(x, toyData, &toyPdf, momPdf, fitPdf);
    delete toyData;
  }
  return 0;
}
