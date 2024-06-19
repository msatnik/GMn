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
Double_t wide_x_range = 0.2;
Double_t nsigma = 0.6;
Double_t xmin = -2;//-1.8; 
Double_t xmax = 1; // 1; 



// Function Declarations
Double_t gausFit(Double_t *x, Double_t *par);
Double_t poly2(Double_t *x, Double_t *par);
Double_t poly3(Double_t *x, Double_t *par);
Double_t poly4(Double_t *x, Double_t *par);
Double_t poly2_2gaus(Double_t *x, Double_t *par);
Double_t poly4_2gaus(Double_t *x, Double_t *par);

// MAIN
void fit_dx_simc(){

// Load rootfiles for histograms
  TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_simc_30p_deen.root"); 
  TFile *f2 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_simc_30p_deep.root"); // Load rootfile

  // Load Histograms
  TH1D *hist3 = (TH1D*)f1->Get("h_dx_HCAL");
  TH1D *hist4 = (TH1D*)f2->Get("h_dx_HCAL");

  // histogram that is a sum of the p and n histograms.
  TH1D *hist1 = (TH1D*)hist3->Clone("hist1");
  hist1->Add(hist4, 1);


  //// Get max bin to find where the proton peak probaby is 
  double maxBin = hist1 -> GetMaximumBin();
  double maxBinValue = hist1 -> GetBinContent(maxBin);
  double maxBinXValue = hist1->GetXaxis()->GetBinCenter(maxBin);
  cout<<"maxBin: "<<maxBin<<endl;
  cout<<"maxBinValue: "<<maxBinValue<<endl;
  cout<<"maxBinXValue: "<<maxBinXValue<<endl;
  cout<<endl;

   /// set up a gaussian fit for the proton peak 
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
  cout<<endl;


  /// set up a gaussian fit for the neutron peak. Mostly the same as proton, but around zero. 
  TF1* gaussian_init2 =  new TF1("gaussian_init2",gausFit , 0-wide_x_range , 0+wide_x_range,3);
  gaussian_init2 ->SetParameters(maxBinValue/2, 0, 0.2); // amp, mean, sigma 
  gaussian_init2 ->SetParLimits(2,0 - wide_x_range, 0 + wide_x_range); // setting limits on neutron peak mean
  /// fit histogram
  gaussian_init2 ->SetParNames("amp", "mean","sigma");
  hist1->Fit(gaussian_init2, "Q", "", 0 - wide_x_range, 0 + wide_x_range);
  //get fit results
  double amp_init2 =  gaussian_init2 -> GetParameter(0);
  double mean_init2 = gaussian_init2 -> GetParameter(1);
  double sigma_init2 = gaussian_init2 ->GetParameter(2);
  cout<<"amp init2: "<<amp_init2<<endl;
  cout<<"mean init2: "<<mean_init2<<endl;
  cout<<"sigma init2: "<<sigma_init2<<endl;
  cout<<endl;



  /// fit again with a more contrained gaussian using the paramters from the looser gaussian for proton peak. 
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
  cout<<endl;


  /// fit again with a more contrained gaussian using the paramters from the looser gaussian for neutron peak. 
  TF1* gaussian2 =  new TF1("gaussian2",gausFit , mean_init2 - nsigma*sigma_init2 , mean_init2 + nsigma*sigma_init2,3);
  gaussian2 ->SetParNames("amp", "mean","sigma");
  gaussian2->SetParameters(amp_init2, mean_init2, sigma_init2); // amp, mean, sigma 
  hist1->Fit(gaussian2, "Q", "", mean_init2 - nsigma*sigma_init2 , mean_init2 + nsigma*sigma_init2);
  //get fit results
  double amp2 =  gaussian2 -> GetParameter(0);
  double mean2 = gaussian2 -> GetParameter(1);
  double sigma2= gaussian2 ->GetParameter(2);
  cout<<"amp2: "<<amp2<<endl;
  cout<<"mean2: "<<mean2<<endl;
  cout<<"sigma2: "<<sigma2<<endl;
  cout<<endl;

  // // Fit the whole distribution to two gaussians plus a polynomail background, using the gauss params we found as initial values. 
  // TF1* poly_2gaus = new TF1("poly_2gaus",poly2_2gaus,xmin,xmax,9);
  // //poly_2gaus ->SetParNames();
  // poly_2gaus ->SetParameters(100,-10,-10,amp,mean,sigma,amp2,mean2,sigma2);
  // poly_2gaus ->SetParLimits(7,mean2 - wide_x_range, mean2 + wide_x_range);
  // hist1->Fit(poly_2gaus, "Q", "", xmin, xmax);
  // double par0 = poly_2gaus -> GetParameter(0);
  // double par0_err = poly_2gaus ->GetParError(0);
  // double par1 = poly_2gaus ->GetParameter(1);
  // double par1_err = poly_2gaus ->GetParError(1);
  // double par2 = poly_2gaus ->GetParameter(2);
  // double par2_err = poly_2gaus ->GetParError(2);
  // double amp_p = poly_2gaus ->GetParameter(3);
  // double amp_p_err = poly_2gaus ->GetParError(3);
  // double mean_p = poly_2gaus ->GetParameter(4);
  // double mean_p_err = poly_2gaus ->GetParError(4);
  // double sigma_p = poly_2gaus ->GetParameter(5);
  // double sigma_p_err = poly_2gaus ->GetParError(5);
  // double amp_n = poly_2gaus ->GetParameter(6);
  // double amp_n_err = poly_2gaus ->GetParError(6);
  // double mean_n = poly_2gaus ->GetParameter(7);
  // double mean_n_err = poly_2gaus ->GetParError(7);
  // double sigma_n = poly_2gaus ->GetParameter(8);
  // double sigma_n_err = poly_2gaus ->GetParError(8);
  // TF1* poly_result = new TF1("poly_result",poly2,xmin,xmax,3);
  // poly_result ->SetParameters(par0,par1,par2);
  // TF1* gaus_result_p = new TF1("gaus_result_p",gausFit, xmin,xmax,3);
  // gaus_result_p ->SetParameters(amp_p, mean_p, sigma_p);
  // TF1* gaus_result_n = new TF1("gaus_result_n",gausFit,xmin,xmax,3);
  // gaus_result_n ->SetParameters(amp_n,mean_n, sigma_n);
  // cout << "Overall Fit"<<endl;
  // cout<<"par0: "<<par0<<" +- "<<par0_err<<endl;
  // cout<<"par1: "<<par1<<" +- "<<par1_err<<endl;
  // cout<<"par2: "<<par2<<" +- "<<par2_err<<endl;
  // cout<<"amp_p: "<<amp_p<<" +- "<<amp_p_err<<endl;
  // cout<<"mean_p: "<<mean_p <<" +- "<<mean_p_err<<endl;
  // cout<<"sigma_p: "<<sigma_p <<" +- "<<sigma_p_err<<endl;
  // cout<<"amp_n: "<<amp_n<<" +- "<<amp_n_err<<endl;
  // cout<<"mean_n: "<<mean_n <<" +- "<<mean_n_err<<endl;
  // cout<<"sigma_n: "<<sigma_n <<" +- "<<sigma_n_err<<endl;
  // cout<<endl;




 // Fit the whole distribution to two gaussians plus a 4th orderpolynomail background, using the gauss params we found as initial values. 
  TF1* poly_2gaus = new TF1("poly_2gaus",poly4_2gaus,xmin,xmax,11);
  //poly_2gaus ->SetParNames();
  poly_2gaus ->SetParameters(100,-10,-10,-10,-10,amp,mean,sigma,amp2,mean2,sigma2);
  poly_2gaus ->SetParLimits(8,amp2 - amp2/4, amp2 + amp2/4); // setting limits on neutron peak amp
  poly_2gaus ->SetParLimits(9,mean2 - wide_x_range, mean2 + wide_x_range); // setting limits on neutron peak mean
  poly_2gaus ->SetParLimits(10,sigma2 - 0.05, sigma2 + 0.05); // setting limits on neutron peak sigma
  poly_2gaus ->SetParLimits(5,amp - amp/4, amp + amp/4); // setting limits on proton peak amp
  poly_2gaus ->SetParLimits(6,mean - wide_x_range, mean + wide_x_range); // setting limits on proton peak mean
  poly_2gaus ->SetParLimits(7,sigma - 0.05, sigma + 0.05); // setting limits on proton peak sigma
  hist1->Fit(poly_2gaus, "Q", "", xmin, xmax);
  double par0 = poly_2gaus -> GetParameter(0);
  double par0_err = poly_2gaus ->GetParError(0);
  double par1 = poly_2gaus ->GetParameter(1);
  double par1_err = poly_2gaus ->GetParError(1);
  double par2 = poly_2gaus ->GetParameter(2);
  double par2_err = poly_2gaus ->GetParError(2);
  double par3 = poly_2gaus ->GetParameter(3);
  double par3_err = poly_2gaus ->GetParError(3);
  double par4 = poly_2gaus ->GetParameter(4);
  double par4_err = poly_2gaus ->GetParError(4);
  double amp_p = poly_2gaus ->GetParameter(5);
  double amp_p_err = poly_2gaus ->GetParError(5);
  double mean_p = poly_2gaus ->GetParameter(6);
  double mean_p_err = poly_2gaus ->GetParError(6);
  double sigma_p = poly_2gaus ->GetParameter(7);
  double sigma_p_err = poly_2gaus ->GetParError(7);
  double amp_n = poly_2gaus ->GetParameter(8);
  double amp_n_err = poly_2gaus ->GetParError(8);
  double mean_n = poly_2gaus ->GetParameter(9);
  double mean_n_err = poly_2gaus ->GetParError(9);
  double sigma_n = poly_2gaus ->GetParameter(10);
  double sigma_n_err = poly_2gaus ->GetParError(10);
  TF1* poly_result = new TF1("poly_result",poly4,xmin,xmax,5);
  poly_result ->SetParameters(par0,par1,par2,par3,par4);
  TF1* gaus_result_p = new TF1("gaus_result_p",gausFit, xmin,xmax,3);
  gaus_result_p ->SetParameters(amp_p, mean_p, sigma_p);
  TF1* gaus_result_n = new TF1("gaus_result_n",gausFit,xmin,xmax,3);
  gaus_result_n ->SetParameters(amp_n,mean_n, sigma_n);
  cout << "Overall Fit"<<endl;
  cout<<"par0: "<<par0<<" +- "<<par0_err<<endl;
  cout<<"par1: "<<par1<<" +- "<<par1_err<<endl;
  cout<<"par2: "<<par2<<" +- "<<par2_err<<endl;
  cout<<"par3: "<<par3<<" +- "<<par3_err<<endl;
  cout<<"par4: "<<par4<<" +- "<<par4_err<<endl;
  cout<<"amp_p: "<<amp_p<<" +- "<<amp_p_err<<endl;
  cout<<"mean_p: "<<mean_p <<" +- "<<mean_p_err<<endl;
  cout<<"sigma_p: "<<sigma_p <<" +- "<<sigma_p_err<<endl;
  cout<<"amp_n: "<<amp_n<<" +- "<<amp_n_err<<endl;
  cout<<"mean_n: "<<mean_n <<" +- "<<mean_n_err<<endl;
  cout<<"sigma_n: "<<sigma_n <<" +- "<<sigma_n_err<<endl;
  cout<<endl;


  // Set up canvas 
  TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
  hist1->SetLineColor(kBlue);
  hist1->Draw();
  //gaussian ->Draw("same");
  //gaussian2 ->Draw("same");
  poly_2gaus ->SetLineColor(kGreen);
  poly_2gaus ->Draw("same");
  poly_result ->SetLineColor(kRed);
  poly_result ->SetLineStyle(2);
  poly_result ->Draw("same");
  gaus_result_p ->SetLineColor(kMagenta); 
  gaus_result_p ->SetLineStyle(2);
  gaus_result_p ->Draw("same");
  gaus_result_n ->SetLineColor(kViolet);
  gaus_result_n ->SetLineStyle(2);
  gaus_result_n ->Draw("same");


  // Create a legend.
  TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9); // Adjust the legend coordinates as needed
  /// Add entries to the legend for the fit parameters
  legend->AddEntry(hist1, "Data", "l");
  legend->AddEntry(poly_2gaus,"Overall Fit","l");
  legend->AddEntry(gaus_result_p,"Proton Peak","l");
  legend->AddEntry("", Form("Mean= %.4f +- %.4f", mean_p,mean_p_err), "");
  legend->AddEntry("", Form("Sigma= %.4f +-%.4f",sigma_p, sigma_p_err), "");
  legend->AddEntry(gaus_result_n,"Neutron Peak","l");
  legend->AddEntry("", Form("Mean= %.4f +- %.4f", mean_n, mean_n_err), "");
  legend->AddEntry("", Form("Sigma= %.4f +- %.4f",sigma_n,sigma_n_err), "");
  legend->AddEntry(poly_result,"Poly Background","l");
  legend->Draw();


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


Double_t poly2(Double_t *x, Double_t *par)
{
  Double_t fit = par[0] + par[1] * x[0] + par[2] * pow(x[0],2);
  return fit;
}// end poly2

Double_t poly3(Double_t *x, Double_t *par)
{
  Double_t fit = par[0] + par[1] * x[0] + par[2] * pow(x[0],2) + par[3] * pow(x[0],3);
  return fit;
}// end poly3

Double_t poly4(Double_t *x, Double_t *par)
{
  Double_t fit = par[0] + par[1] * x[0] + par[2] * pow(x[0],2) + par[3] * pow(x[0],3) + par[4] * pow(x[0],4) ;
  return fit;
}


Double_t poly2_2gaus(Double_t *x, Double_t *par)
{
   Double_t height = par[3];
   Double_t mean = par[4];
   Double_t sigma = par[5];

   Double_t height2 = par[6];
   Double_t mean2 = par[7];
   Double_t sigma2 = par[8];

   Double_t fit = par[0] + par[1] * x[0] + par[2] * pow(x[0],2) +  height * exp(-0.5 * pow((x[0] - mean)/sigma,2)) + height2 * exp(-0.5 * pow((x[0] - mean2)/sigma2,2)) ;
   return fit; 
}

Double_t poly4_2gaus(Double_t *x, Double_t *par)
{
   Double_t height = par[5];
   Double_t mean = par[6];
   Double_t sigma = par[7];

   Double_t height2 = par[8];
   Double_t mean2 = par[9];
   Double_t sigma2 = par[10];

   Double_t fit = par[0] + par[1] * x[0] + par[2] * pow(x[0],2) + par[3] * pow(x[0],3) +  par[4] * pow(x[0],4) +  height * exp(-0.5 * pow((x[0] - mean)/sigma,2)) + height2 * exp(-0.5 * pow((x[0] - mean2)/sigma2,2)) ;
   return fit; 
}


