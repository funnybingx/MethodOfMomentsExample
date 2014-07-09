
// 
#include <iostream>
#include <vector>
#include "boost/program_options.hpp"

// ROOT
#include <TFile.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TAxis.h>

// RooFit
#include <RooPlot.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooRealVar.h>
#include <RooAddPdf.h>
#include <RooDataSet.h>
#include <RooChebychev.h>
#include <RooFit.h>
#include <RooRandom.h>
#include <RooMsgService.h> // make RooFit quiet down

// MoM framework
#include "utils/help.hpp"
#include "mom/MomentsCalculator.hpp"
#include "mom/LegendreMomentsCalculator.hpp"
#include "mom/FigureOfMeritCalculator.hpp"

namespace rf = RooFit; 
using std::cout;
using std::cerr;
using std::endl;
using std::vector;

struct config // configuration for this executable
{
  bool debug;
  int ndata;  
  int ntoys;
  bool fit;
  bool save;
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
      ("debug", po::bool_switch(&c.debug)->default_value(false), "run in debug")
      ("ndata", po::value<int>(&c.ndata)->default_value(5000), "toy generated stats")
      ("ntoys", po::value<int>(&c.ntoys)->default_value(1), "number of toys to run")
      ("fit", po::bool_switch(&c.fit), "also run a fit over each toy")
      ("save", po::bool_switch(&c.save)->default_value(false), "save dataset with bad pvalue")  
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

void binLogX(TH1D*& h)
{
  //---------------------------------------------------------------------------
  // Bins a histogram logarithmically on the x-axis. Code stolen from:
  // http://root.cern.ch/root/roottalk/roottalk06/1213.html
  //---------------------------------------------------------------------------
  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();
  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = (to - from) / bins;
  Axis_t *new_bins = new Axis_t[bins + 1];
  for (int i = 0; i <= bins; i++) 
    new_bins[i] = TMath::Power(10, from + i * width);
  axis->Set(bins, new_bins);
  delete new_bins;
  return;
} 

void plot(TH1D* orig, TH1D* mom, TH1D* mom2, TH1D* fitted, TString title)
{
  gStyle->SetOptStat(0);
  TCanvas canvas;
  TLegend leg(0.1,0.7,0.48,0.9);
  orig->SetLineColor(kBlue);
  mom->SetLineColor(kRed);
  mom2->SetLineColor(kOrange);
  if (fitted) fitted->SetLineColor(kGreen);
  orig->Draw();
  mom->Draw("same");
  mom2->Draw("same");
  if (fitted) fitted->Draw("same");
  if (title=="pvals") {
    canvas.SetLogx();
    canvas.SetLogy();
  }
  leg.AddEntry(orig, "original pdf", "l");
  leg.AddEntry(mom, "MoM - Chebyshevs", "l");
  leg.AddEntry(mom2, "MoM - Legendres", "l");
  if (fitted) leg.AddEntry(fitted, "RooFit", "l");
  leg.Draw();
  TString name = "plots/"+title;
  canvas.SaveAs(name+".pdf");
  return;
}

void plot(RooRealVar* x, RooDataSet* d, RooAbsPdf* orig, RooAbsPdf* mom,
          RooAbsPdf* mom2, RooAbsPdf* fitted=NULL, int itoy=0)
{
  TCanvas canvas;
  TLegend leg(0.1,0.7,0.48,0.9);
  RooPlot* frame = x->frame(50);
  d->plotOn(frame, rf::Name("data"));
  orig->plotOn(frame, rf::LineColor(kBlue), rf::Name("orig"));
  mom->plotOn(frame, rf::LineColor(kRed), rf::Name("mom"));
  mom2->plotOn(frame, rf::LineColor(kOrange), rf::Name("mom2"));
  if (fitted) fitted->plotOn(frame, rf::LineColor(kGreen), rf::Name("fitted"));
  frame->Draw();
  frame->SetTitle("Sam\'s 1D MoM implementation");
  leg.AddEntry("data", "toy data", "p");
  leg.AddEntry("orig", "original pdf", "l");
  leg.AddEntry("mom", "MoM - Chebyshevs", "l");
  leg.AddEntry("mom2", "MoM - Legendres", "l");
  if (fitted) leg.AddEntry("fitted", "RooFit", "l");
  leg.Draw();
  TString name = "plots/toy";
  name+=itoy;
  canvas.SaveAs(name+".pdf");
  delete frame;
  return;
}

void runFit(RooChebychev*& fitPdf, RooRealVar* x, 
            vector<RooRealVar*>& fitCoeffs, RooDataSet* toyData, bool stfu=false)
{
  RooArgList ral;
  for (unsigned int icoef=0; icoef<fitCoeffs.size(); ++icoef) {
    TString name = "fc"; name+=icoef;
    fitCoeffs[icoef] = new RooRealVar(name, name, 0.5, -1.0, 1.0);
    //fitCoeffs[icoef]->Print();
    ral.add(*fitCoeffs[icoef]);
  }
  fitPdf = new RooChebychev("fitpdf", "Fitted Cheby pdf", *x, ral);
  if (stfu) fitPdf->fitTo( *toyData, rf::PrintLevel(-1000) );
  else fitPdf->fitTo( *toyData );
  return;
}

int main(int argc, char *argv[])
{
  // deal with options parsing
  config c;
  int returnCode = parseOptions(c, argc, argv);
  if (returnCode!= 0) { return returnCode; }
  if (c.ntoys > 200) {
    // then this is sufficiently many toys that I don't want to see any
    // RooFit stdout 
    RooMsgService::instance().setStreamStatus(1,false);
    RooMsgService::instance().setStreamStatus(0,false);
  }

  // the 1D variable
  RooRealVar* x = new RooRealVar("x", "x", -1, 1);

  // build a toy pdf for toy data
  RooRealVar c1("c1", "linear order coeff",0.2);
  RooRealVar c2("c2", "quadratic order coeff", 0.5);
  RooChebychev toyPdf("toypdf","Toy pdf polynomial", *x, RooArgList(c1,c2));
  RooRandom::randomGenerator()->SetSeed(123456789);

  // histograms for FoM's of many toys
  TH1D* hpval_orig = new TH1D("hpval_orig","P-values",100,-6,0); 
  TH1D* hpval_lgnd = new TH1D("hpval_lgnd","P-values",100,-6,0); 
  TH1D* hpval_cheb = new TH1D("hpval_cheb","P-values",100,-6,0); 
  TH1D* hpval_fit  = new TH1D("hpval_fit","P-values",100,-6,0); 
  binLogX(hpval_orig); // make the x-axis binning logarithmic
  binLogX(hpval_lgnd);
  binLogX(hpval_cheb);
  binLogX(hpval_fit);
  TH1D* hchi2_orig = new TH1D("hchi2_orig", "Distributions of #chi^{2}/ndof", 100, 0.0, 2.5);
  TH1D* hchi2_lgnd = new TH1D("hchi2_lgnd", "Distributions of #chi^{2}/ndof", 100, 0.0, 2.5);
  TH1D* hchi2_cheb = new TH1D("hchi2_cheb", "Distributions of #chi^{2}/ndof", 100, 0.0, 2.5);
  TH1D* hchi2_fit  = new TH1D("hchi2_fit", "Distributions of #chi^{2}/ndof", 100, 0.0, 2.5);

  // declare tree storing FoM information // TODO make this cleverer in some way
  double chi2_fit=0.,chi2_cheb=0.,chi2_lgnd=0.,chi2_orig=0.;
  double pval_fit=0.,pval_cheb=0.,pval_lgnd=0.,pval_orig=0.;
  const int const_order=c.ordermom;
  double cheb_vars[const_order][const_order];
  double cheb_vals[const_order];
  double lgnd_vars[const_order][const_order];
  double lgnd_vals[const_order];
  TTree fomTree("fomTree","fomTree");
  fomTree.Branch("pval_orig", &pval_orig);
  fomTree.Branch("pval_cheb",  &pval_cheb);
  fomTree.Branch("pval_lgnd",  &pval_lgnd);
  fomTree.Branch("pval_fit",  &pval_fit);
  fomTree.Branch("chi2_orig", &chi2_orig);
  fomTree.Branch("chi2_cheb",  &chi2_cheb);
  fomTree.Branch("chi2_lgnd",  &chi2_lgnd);
  fomTree.Branch("chi2_fit",  &chi2_fit);

  // declaring branches for the moments and the variances
  for(int iv=0;iv<c.ordermom;++iv)
  {
    char brName[512];
    sprintf(brName,"cheb_val_%d",iv);
    fomTree.Branch(brName,&cheb_vals[iv]);
    sprintf(brName,"lgnd_val_%d",iv);
    fomTree.Branch(brName,&lgnd_vals[iv]);
    for(int jv=0;jv<c.ordermom;++jv)
    {
      sprintf(brName,"cheb_var_%d_%d",iv,jv);
      fomTree.Branch(brName,&cheb_vars[iv][jv]);
      sprintf(brName,"lgnd_var_%d_%d",iv,jv);
      fomTree.Branch(brName,&lgnd_vars[iv][jv]);
    }
  }

  // loop and generate toys
  for (int itoy=0; itoy<c.ntoys; ++itoy) {
    RooDataSet* toyData = toyPdf.generate(RooArgSet(*x), c.ndata);

    // run the method of moment
    MomentsCalculator moments(c.ordermom, x);
    if (c.debug) moments.setDebug();
    moments.run(*toyData);
    RooChebychev* momPdf = moments.getRooPdf();
    if (c.debug) momPdf->Print();

    // get the moments and variances
    vector<double> moment_vals = moments.getMoments();
    vvd variances              = moments.getVariances();
    for(int iv=0;iv<c.ordermom;++iv)
    {
      cheb_vals[iv]=moment_vals[iv];
      for(int jv=0;jv<c.ordermom;++jv)
        cheb_vars[iv][jv]=variances[iv][jv];
    }

    // run the method of moments using the Legendre polynomials
    LegendreMomentsCalculator legendre(c.ordermom, x);
    if (c.debug) moments.setDebug();
    legendre.run(*toyData);
    RooAddPdf* legendrePdf = legendre.getRooPdf();
    if (c.debug) legendrePdf->Print();
  
    // get the moments and variances
    vector<double> legendre_moment_vals = legendre.getMoments();
    vvd legendre_variances              = legendre.getVariances();
    for(int iv=0;iv<c.ordermom;++iv)
    {
      lgnd_vals[iv]=legendre_moment_vals[iv];
      for(int jv=0;jv<c.ordermom;++jv)
        lgnd_vars[iv][jv]=legendre_variances[iv][jv];
    }

    // also run a RooFit fit 
    RooChebychev* fitPdf=NULL;
    vector<RooRealVar*> fitCoeffs ( c.ordermom, NULL );
    if (c.fit) {
      if(c.debug)
        cout << "run::main INFO: will now declare pdf and run a fit" << endl;
      runFit(fitPdf, x, fitCoeffs, toyData, (c.ntoys>200));
    }

    FigureOfMeritCalculator fom(50, toyData, x);    
    if(c.debug){
      cout << "run::main INFO: the FoM w.r.t the original gen pdf" << endl;
      fom.print(&toyPdf);
    }
    chi2_orig = fom.chi2Ndof(&toyPdf);
    pval_orig = fom.pvalue(&toyPdf);
    hpval_orig->Fill( pval_orig );
    hchi2_orig->Fill( chi2_orig );

    if(c.debug){
      cout << "run::main info: the fom w.r.t the mom pdf" << endl;
      fom.print(momPdf);
    }
     chi2_cheb = fom.chi2Ndof(momPdf);
     pval_cheb = fom.pvalue(momPdf);
    hpval_cheb->Fill( pval_cheb );
    hchi2_cheb->Fill( chi2_cheb );
    if(c.debug){
      cout << "run::main info: the fom w.r.t the legendre pdf" << endl;
      fom.print(legendrePdf);
    }
     chi2_lgnd = fom.chi2Ndof(legendrePdf);
     pval_lgnd = fom.pvalue(legendrePdf);
    hpval_lgnd->Fill( pval_lgnd );
    hchi2_lgnd->Fill( chi2_lgnd );
    if (fitPdf) {
      if(c.debug){
        cout << "run::main INFO: the FoM w.r.t the fitted pdf" << endl;
        fom.print(fitPdf);
      }
      chi2_fit = fom.chi2Ndof(fitPdf);
      pval_fit = fom.pvalue(fitPdf);
      hpval_fit->Fill( pval_fit );
      hchi2_fit->Fill( chi2_fit );
    }
    if(pval_cheb<0.0001  || pval_lgnd<0.0001){
      plot(x, toyData, &toyPdf, momPdf, legendrePdf, fitPdf, itoy);
      if(c.save){
        char dataFileName[512];
        sprintf(dataFileName,"plots/bad_pvalue_dataset_%d.root",itoy);
        toyData->SaveAs(dataFileName);
      }
    }
    fomTree.Fill();

    // free up mem
    for (unsigned int icoeff=0; icoeff<fitCoeffs.size(); ++icoeff)
      delete fitCoeffs[icoeff];
    delete momPdf; delete fitPdf; delete toyData; delete legendrePdf;

    // if running for many toys, then just show a percentage completed
    // otherwise there'll be plenty of other stdout each iteration
    if ((c.ntoys>200) && (itoy%100 == 0)) utils::help::reportOnLoop(itoy, c.ntoys);
  }
  if (!c.debug) cout << endl;

  // plot the histograms
  plot(hpval_orig, hpval_cheb, hpval_lgnd, hpval_fit, "pvals");
  plot(hchi2_orig, hchi2_cheb, hchi2_lgnd, hchi2_fit, "chi2s");
  TFile output("plots/histograms.root", "RECREATE");
  output.cd();
  hpval_orig->Write(); hpval_cheb->Write(); hpval_fit->Write();
  hchi2_orig->Write(); hchi2_cheb->Write(); hchi2_fit->Write();
  hpval_lgnd->Write(); hchi2_lgnd->Write();
  fomTree.Write();
  output.Close();
  delete hpval_orig; delete hpval_cheb; delete hpval_lgnd; delete hpval_fit;
  delete hchi2_orig; delete hchi2_cheb; delete hchi2_lgnd; delete hchi2_fit;
  return 0;
}
