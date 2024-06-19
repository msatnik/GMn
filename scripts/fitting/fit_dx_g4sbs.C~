#include <vector>
#include <iostream>
#include <algorithm>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include <fstream>
#include <string>
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TF1.h"

// Global Params
Double_t wide_x_range = 0.4;
Double_t nsigma = 0.6;


// Function Declarations
Double_t gausFit(Double_t *x, Double_t *par);

// MAIN
void fit_dx_g4sbs(){

  // Load rootfiles for histograms
  TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_sbs30p_sf0p359.root"); 

  // Load Histograms
  TH1D *hist1 = (TH1D*)f1->Get("h_dx_HCAL");


  //// Get max bin to find where the proton peak probaby is 
  double maxBin = hist1 -> GetMaximumBin();
  double maxBinValue = hist1 -> GetBinContent(maxBin);
  double maxBinXValue = hist1->GetXaxis()->GetBinCenter(maxBin);
  cout<<"maxBin: "<<maxBin<<endl;
  cout<<"maxBinValue: "<<maxBinValue<<endl;
  cout<<"maxBinXValue: "<<maxBinXValue<<endl;
  /// set up a gaussian fit
  TF1* gaussian_init =  new TF1("gaussian_init",gausFit , maxBinXValue -wide_x_range , maxBinXValue + wide_x_range,3);
  gaussian_init ->SetParameters(maxBinValue, maxBinXValue, 0.2); // amp, mean, sigma 
  /// fit histogram
  gaussian_init ->SetParNames("amp", "mean","sigma");
  hist1->Fit(gaussian_init, "Q", "", maxBinXValue - wide_x_range, maxBinXValue + wide_x_range);
  //get fit results
  double amp_init =  gaussian_init -> GetParameter(0);
  double mean_init = gaussian_init -> GetParameter(1);
  double sigma_init = gaussian_init ->GetParameter(2);
  cout<<"amp init: "<<amp_init<<endl;
  cout<<"mean init: "<<mean_init<<endl;
  cout<<"sigma init: "<<sigma_init<<endl;
  /// fit again with a more contrained gaussian using the paramters from the looser gaussian. 
  TF1* gaussian =  new TF1("gaussian",gausFit , mean_init - nsigma*sigma_init , mean_init + nsigma*sigma_init,3);
  gaussian ->SetParNames("amp", "mean","sigma");
  gaussian->SetParameters(amp_init, mean_init, sigma_init); // amp, mean, sigma 
  hist1->Fit(gaussian, "Q", "", mean_init - nsigma*sigma_init , mean_init + nsigma*sigma_init);
  //get fit results
  double amp =  gaussian -> GetParameter(0);
  double mean = gaussian -> GetParameter(1);
  double sigma = gaussian ->GetParameter(2);
  cout<<"amp: "<<amp<<endl;
  cout<<"mean: "<<mean<<endl;
  cout<<"sigma: "<<sigma<<endl;



  // Set up canvas 
  TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
  hist1->SetLineColor(kBlue);
  hist1->Draw();
  gaussian ->Draw("same");



    // Create a legend.
    //TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9); // Adjust the legend coordinates as needed
   // Add entries to the legend for the fit parameters
    //legend->AddEntry(hist1, "Neutron", "l");
    //legend->Draw();

    c1->Update();
}



Double_t gausFit(Double_t *x, Double_t *par)
{
  Double_t height = par[0];
  Double_t mean = par[1];
  Double_t sigma = par[2];

  Double_t fit = height * exp(-0.5 * pow((x[0] - mean)/sigma,2));
  return fit;
}
