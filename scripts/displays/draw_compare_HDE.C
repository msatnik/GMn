#include <TSystem.h>
#include <TChain.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
#include <stack>
#include "TLorentzVector.h"
#include "TCut.h"
#include "TLatex.h"
#include "TLine.h"
#include <string>
#include <algorithm>
using namespace std;

#include "/work/halla//sbs/msatnik/GMn/classes/FitHistogramGeneral.h"
#include "/work/halla//sbs/msatnik/GMn/classes/FitHistogramGeneral.cpp"

#include "/work/halla//sbs/msatnik/GMn/classes/Utility.h"
#include "/work/halla//sbs/msatnik/GMn/classes/Utility.cpp"

//// default values to be overwritten 
TString DataFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_2Dhistos.root";
std::pair<double,double> W2_range(0.66, 1.1);
std::pair<double,double> hcal_x_exp_range(-1.15,0.55);
std::pair<double,double> hcal_y_exp_range(-0.3,0.3);
std::pair<double,double> hcal_dy_range(-0.3,0.3);
double hcal_e_min = 0.04;
std::pair<double,double> coin_range(-10,10);
double ps_e_min = 0.2;
double ps_sh_e_min = 3;
double grinch_clus_size_min = 3;
double grinch_clus_adc_min = 25;
std::pair<double,double> e_over_p_range(0.82,1.14);//abs(e_over_p-0.98)<0.16
std::pair<double,double> vz_range(-0.06,0.06);//


// SBS4 30p 
TString DataFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_SpotStudy_sept26.root";
std::pair<double,double> W2_range_sbs4_30p(0.66, 1.1);
std::pair<double,double> hcal_x_exp_range_sbs4_30p(-1.3,0.54);
std::pair<double,double> hcal_y_exp_range_sbs4_30p(-0.35,0.4);
std::pair<double,double> hcal_dy_range_sbs4_30p(-0.3,0.3);
double hcal_e_min_sbs4_30p = 0.02;
std::pair<double,double> coin_range_sbs4_30p(-10,10);
double ps_e_min_sbs4_30p = 0.2;
double ps_sh_e_min_sbs4_30p = 1.7;
double grinch_clus_size_min_sbs4_30p = 3;
double grinch_clus_adc_min_sbs4_30p = 25;
std::pair<double,double> e_over_p_range_sbs4_30p(0.78,1.18);//abs(e_over_p-0.98)<0.2
std::pair<double,double> vz_range_sbs4_30p(-0.06,0.06);//


// SBS4 50p 
TString DataFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_SpotStudy_sept26.root";
std::pair<double,double> W2_range_sbs4_50p(0.66, 1.1);
std::pair<double,double> hcal_x_exp_range_sbs4_50p(-0.83,0.54);
std::pair<double,double> hcal_y_exp_range_sbs4_50p(-0.35,0.4);
std::pair<double,double> hcal_dy_range_sbs4_50p(-0.3,0.3);
double hcal_e_min_sbs4_50p = 0.02;
std::pair<double,double> coin_range_sbs4_50p(-10,10);
double ps_e_min_sbs4_50p = 0.2;
double ps_sh_e_min_sbs4_50p = 1.7;
double grinch_clus_size_min_sbs4_50p = 3;
double grinch_clus_adc_min_sbs4_50p = 25;
std::pair<double,double> e_over_p_range_sbs4_50p(0.78,1.18);//abs(e_over_p-0.98)<0.16
std::pair<double,double> vz_range_sbs4_50p(-0.06,0.06);//

// SBS8 70p
TString DataFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_SpotStudy_sept26.root";
std::pair<double,double> W2_range_sbs8_70p(0.66, 1.1);
std::pair<double,double> hcal_x_exp_range_sbs8_70p(-1.27,0.63);
std::pair<double,double> hcal_y_exp_range_sbs8_70p(-0.4,0.35);
std::pair<double,double> hcal_dy_range_sbs8_70p(-0.3,0.3);
double hcal_e_min_sbs8_70p = 0.04;
std::pair<double,double> coin_range_sbs8_70p(-10,10);
double ps_e_min_sbs8_70p = 0.2;
double ps_sh_e_min_sbs8_70p = 3;
double grinch_clus_size_min_sbs8_70p = 3;
double grinch_clus_adc_min_sbs8_70p = 25;
std::pair<double,double> e_over_p_range_sbs8_70p(0.82,1.14);//abs(e_over_p-0.98)<0.16
std::pair<double,double> vz_range_sbs8_70p(-0.06,0.06);//

// SBS9 70p
TString DataFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_SpotStudy_sept26_1.root";
std::pair<double,double> W2_range_sbs9_70p(0.66, 1.1);
std::pair<double,double> hcal_x_exp_range_sbs9_70p(-0.84,0.27);
std::pair<double,double> hcal_y_exp_range_sbs9_70p(-0.35,0.35);
std::pair<double,double> hcal_dy_range_sbs9_70p(-0.2,0.2);
double hcal_e_min_sbs9_70p = 0.04;
std::pair<double,double> coin_range_sbs9_70p(-10,10);
double ps_e_min_sbs9_70p = 0.2;
double ps_sh_e_min_sbs9_70p = 3;
double grinch_clus_size_min_sbs9_70p = 3;
double grinch_clus_adc_min_sbs9_70p = 25;
std::pair<double,double> e_over_p_range_sbs9_70p(0.82,1.14);//abs(e_over_p-0.98)<0.16
std::pair<double,double> vz_range_sbs9_70p(-0.06,0.06);//


void draw_compare_HDE(TString KineString="sbs4_0p"){// main

  gStyle->SetNumberContours(255); 
  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0);
  gStyle->SetStatTextColor(kBlack);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.06);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
 


  Utility utilityHandler;

  // if (KineString == "sbs4_30p"){
  //   DataFileString = DataFileString_sbs4_30p;
  //   W2_range=W2_range_sbs4_30p;
  //   hcal_x_exp_range=hcal_x_exp_range_sbs4_30p;
  //   hcal_y_exp_range = hcal_y_exp_range_sbs4_30p;
  //   hcal_dy_range = hcal_dy_range_sbs4_30p;
  //   hcal_e_min = hcal_e_min_sbs4_30p;
  //   coin_range = coin_range_sbs4_30p;
  //   ps_e_min = ps_e_min_sbs4_30p;
  //   ps_sh_e_min = ps_sh_e_min_sbs4_30p;
  //   grinch_clus_size_min = grinch_clus_size_min_sbs4_30p;
  //   grinch_clus_adc_min = grinch_clus_adc_min_sbs4_30p;
  //   e_over_p_range= e_over_p_range_sbs4_30p;
  //   vz_range = vz_range_sbs4_30p;
  // }else if (KineString == "sbs4_50p"){
  //   DataFileString = DataFileString_sbs4_50p;
  //   W2_range=W2_range_sbs4_50p;
  //   hcal_x_exp_range=hcal_x_exp_range_sbs4_50p;
  //   hcal_y_exp_range = hcal_y_exp_range_sbs4_50p;
  //   hcal_dy_range = hcal_dy_range_sbs4_50p;
  //   hcal_e_min = hcal_e_min_sbs4_50p;
  //   coin_range = coin_range_sbs4_50p;
  //   ps_e_min = ps_e_min_sbs4_50p;
  //   ps_sh_e_min = ps_sh_e_min_sbs4_50p;
  //   grinch_clus_size_min = grinch_clus_size_min_sbs4_50p;
  //   grinch_clus_adc_min = grinch_clus_adc_min_sbs4_50p;
  //   e_over_p_range= e_over_p_range_sbs4_50p;
  //   vz_range = vz_range_sbs4_50p;
  // }else if (KineString == "sbs8_70p"){
  //   DataFileString = DataFileString_sbs8_70p;
  //   W2_range=W2_range_sbs8_70p;
  //   hcal_x_exp_range=hcal_x_exp_range_sbs8_70p;
  //   hcal_y_exp_range = hcal_y_exp_range_sbs8_70p;
  //   hcal_dy_range = hcal_dy_range_sbs8_70p;
  //   hcal_e_min = hcal_e_min_sbs8_70p;
  //   coin_range = coin_range_sbs8_70p;
  //   ps_e_min = ps_e_min_sbs8_70p;
  //   ps_sh_e_min = ps_sh_e_min_sbs8_70p;
  //   grinch_clus_size_min = grinch_clus_size_min_sbs8_70p;
  //   grinch_clus_adc_min = grinch_clus_adc_min_sbs8_70p;
  //   e_over_p_range= e_over_p_range_sbs8_70p;
  //   vz_range = vz_range_sbs8_70p;
  // }else if (KineString == "sbs9_70p"){
  //   DataFileString = DataFileString_sbs9_70p;
  //   W2_range=W2_range_sbs9_70p;
  //   hcal_x_exp_range=hcal_x_exp_range_sbs9_70p;
  //   hcal_y_exp_range = hcal_y_exp_range_sbs9_70p;
  //   hcal_dy_range = hcal_dy_range_sbs9_70p;
  //   hcal_e_min = hcal_e_min_sbs9_70p;
  //   coin_range = coin_range_sbs9_70p;
  //   ps_e_min = ps_e_min_sbs9_70p;
  //   ps_sh_e_min = ps_sh_e_min_sbs9_70p;
  //   grinch_clus_size_min = grinch_clus_size_min_sbs9_70p;
  //   grinch_clus_adc_min = grinch_clus_adc_min_sbs9_70p;
  //   e_over_p_range= e_over_p_range_sbs9_70p;
  //   vz_range = vz_range_sbs9_70p;
  // }else {
  //   std::cout<<"Error with kinematic setting"<<std::endl;
  //   return;
  // }
  
 
  // cout<<endl;
  // cout<<"What we're running: "<<KineString<<endl;
  // cout<<DataFileString<<endl;
  // cout<<endl;

 
  // load files 
  TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_0p_cuts_HDE_oct11.root"); // Deuterium
  TFile *f2 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_0p_cuts_LH2_HDE_oct14.root"); // Hydrogen 14 has the average mass

  if (!f1||!f2) {
    std::cout<<"Error with loading one of the rootfiles"<<std::endl;
    return;
  }

  // set location and name for output file 
  TString outputfilelocation="../../output/compare_HDE/"+KineString;
  TString outputfilename = outputfilelocation +"/"+ KineString+".root";
  // Declare outfile
  TFile *fout = new TFile(outputfilename,"RECREATE");
  cout<<"writing to file: "<< outputfilename <<endl;

  // Load Histograms

  TH1D *d2_HDE_hcal_x_exp_hist = (TH1D*)f1->Get("HDE_hcal_x_exp_hist");
  TH1D *h2_HDE_hcal_x_exp_hist = (TH1D*)f2->Get("HDE_hcal_x_exp_hist");
  
  TH1D *d2_HDE_hcal_y_exp_hist = (TH1D*)f1->Get("HDE_hcal_y_exp_hist");
  TH1D *h2_HDE_hcal_y_exp_hist = (TH1D*)f2->Get("HDE_hcal_y_exp_hist");

  /// Removing any leftover fits from when these were made
  utilityHandler.RemoveExistingFit(d2_HDE_hcal_x_exp_hist);
  utilityHandler.RemoveExistingFit(h2_HDE_hcal_x_exp_hist);
   
  utilityHandler.RemoveExistingFit(d2_HDE_hcal_y_exp_hist);
  utilityHandler.RemoveExistingFit(h2_HDE_hcal_y_exp_hist);
   

  
  TH1D *d2_HDE_hcal_x_exp_hist_dont_draw = (TH1D*)d2_HDE_hcal_x_exp_hist->Clone("d2_HDE_hcal_x_exp_hist_dont_draw");
  TH1D *h2_HDE_hcal_x_exp_hist_dont_draw = (TH1D*)h2_HDE_hcal_x_exp_hist->Clone("h2_HDE_hcal_x_exp_hist_dont_draw");
  
  TH1D *d2_HDE_hcal_y_exp_hist_dont_draw = (TH1D*)d2_HDE_hcal_y_exp_hist->Clone("d2_HDE_hcal_y_exp_hist_dont_draw");
  TH1D *h2_HDE_hcal_y_exp_hist_dont_draw = (TH1D*)h2_HDE_hcal_y_exp_hist->Clone("h2_HDE_hcal_y_exp_hist_hist_dont_draw");

    

  d2_HDE_hcal_x_exp_hist->SetTitle("HCal Detection Efficiency");
  d2_HDE_hcal_x_exp_hist->SetYTitle("HCal Detection Efficiency");
  h2_HDE_hcal_x_exp_hist->SetTitle("HCal Detection Efficiency");
  h2_HDE_hcal_x_exp_hist->SetYTitle("HCal Detection Efficiency");

  // double x_exp_max = std::max(h2_HDE_hcal_x_exp_hist->GetXaxis()->GetXmax(), d2_HDE_hcal_x_exp_hist->GetXaxis()->GetXmax());
  //  double x_exp_min = std::min(h2_HDE_hcal_x_exp_hist->GetXaxis()->GetXmin(),d2_HDE_hcal_x_exp_hist->GetXaxis()->GetXmin());

  //  double y_exp_max = std::max(h2_HDE_hcal_y_exp_hist->GetXaxis()->GetXmax(), d2_HDE_hcal_y_exp_hist->GetXaxis()->GetXmax());
  // double y_exp_min = std::min(h2_HDE_hcal_y_exp_hist->GetXaxis()->GetXmin(),d2_HDE_hcal_y_exp_hist->GetXaxis()->GetXmin());

  double x_exp_min = -1.5;
  double x_exp_max = 1.5;
  double y_exp_min = -0.65;
  double y_exp_max = 0.8;

  cout<<"x_exp: "<<x_exp_min<<" , "<<x_exp_max<<endl;
  cout<<"y_exp: "<<y_exp_min<<" , "<<y_exp_max<<endl;
    
  FitHistogramGeneral fitHandlerGen_x_exp(d2_HDE_hcal_x_exp_hist, nullptr,-1, 1);

  fitHandlerGen_x_exp.fitData(h2_HDE_hcal_x_exp_hist_dont_draw);// bug in here that corrupts it and crashes if you try and draw

  TF1 *fit_result_x_exp = (TF1*)fitHandlerGen_x_exp.fitFunc->Clone("fit_result_x_exp");// bug in here that corrupts it and crashes if you try and draw

  

  double scale_x_exp = fitHandlerGen_x_exp.scale1;
  double scale_x_exp_err = fitHandlerGen_x_exp.scale1_err;
  double scale2_x_exp = fitHandlerGen_x_exp.scale2;// should be zero
  double scale2_x_exp_err = fitHandlerGen_x_exp.scale2_err; // shoudl be zero
  double ChiSq_x_exp = fitHandlerGen_x_exp.ChiSq;
  double ndf_x_exp = fitHandlerGen_x_exp.NDF;

  if (scale2_x_exp != 0){
    cout<<"Error: the scale on the ``second'' fit histogram should be zero because we are not using it" <<endl;
  }

  cout<<endl<<"scale for d2 over x_exp = "<< scale_x_exp <<" +/- "<<scale_x_exp_err<<endl;
  cout<<endl;

  /// Scale the d2 histogram using the fit result 
  TH1D* hist_result_x_exp = utilityHandler.ScaleAndShiftHistogram(d2_HDE_hcal_x_exp_hist,scale_x_exp ,0);
  hist_result_x_exp->SetLineColor(kGreen+2);

  /// Subtract the scaled d2 from the h2 histos.
  TH1D* hist_resid_x_exp = (TH1D*)h2_HDE_hcal_x_exp_hist->Clone("hist_resid_x_exp");
  hist_resid_x_exp ->Add(hist_result_x_exp,-1);

  hist_resid_x_exp ->GetXaxis() ->SetRangeUser(x_exp_min, x_exp_max);
  // hist_resid_x_exp ->GetYaxis() ->SetRangeUser(-0.5,0.5);
  

  TCanvas *testc = new TCanvas("testc","testc",800,600);
  hist_result_x_exp->Draw("E");
  h2_HDE_hcal_x_exp_hist->Draw("same E");
  hist_resid_x_exp->Draw("E same");


  double markersize =1;
  
  // Convert histogram to TGraphErrors
  TGraphErrors *h2_graph_x_exp = utilityHandler.histogramToGraphErrors(h2_HDE_hcal_x_exp_hist); 
  utilityHandler.customizeGraph(h2_graph_x_exp, 20, kMagenta+2, markersize); // assigns shape, color, and marker size to the points
  h2_graph_x_exp->GetXaxis() ->SetRangeUser(x_exp_min, x_exp_max);
  utilityHandler.RemoveZeroPoints(h2_graph_x_exp);
 
  // Convert histogram to TGraphErrors
  TGraphErrors *result_graph_x_exp = utilityHandler.histogramToGraphErrors(hist_result_x_exp); 
  utilityHandler.customizeGraph(result_graph_x_exp, 20, kGreen+2, markersize); // assigns shape, color, and marker size to the points
  hist_result_x_exp->GetXaxis() ->SetRangeUser(x_exp_min, x_exp_max);

  utilityHandler.RemoveZeroPoints(result_graph_x_exp);
  
  TCanvas *graphcanvas = new TCanvas("graphcanvas","graphcanvas",800,600);   
  //utilityHandler.adjustCanvas(graphcanvas);
  graphcanvas->Divide(1,2);
  // Upper pad
  graphcanvas->cd(1);
  TPad *upperPad1 = (TPad*)gPad;
  upperPad1->SetGrid();
  upperPad1->SetPad(0.01, 0.35, 0.99, 0.99); //Double_t xlow, Double_t ylow, Double_t xup, Double_t yup
  //utilityHandler.adjustPad(upperPad1);
  //upperPad1->SetFillColor(20);
  //upperPad1->SetFrameFillColor(20);
  // Lower pad
  graphcanvas->cd(2);
  TPad *lowerPad1 = (TPad*)gPad;
  lowerPad1->SetGrid();
  lowerPad1->SetPad(0.01, 0.01, 0.99, 0.34);
  //utilityHandler.adjustPad(lowerPad1);
  //  lowerPad->SetFillColor(18);
  //lowe1rPad->SetFrameFillColor(18);

  // utilityHandler.AdjustHistLabelOffset(h2_HDE_hcal_x_exp_hist);
  // utilityHandler.AdjustHistLabelOffset(d2_HDE_hcal_x_exp_hist);
  // utilityHandler.AdjustGraphLabelOffset( h2_graph_x_exp);
  // utilityHandler.AdjustGraphLabelOffset(result_graph_x_exp);


  upperPad1->cd();
  h2_graph_x_exp ->SetTitle("x exp");
  result_graph_x_exp->Draw("AP");
  h2_graph_x_exp->Draw("P SAME"); 
  h2_graph_x_exp->GetXaxis()->SetTitle("hcal_dx");

  int nEntries_data =h2_HDE_hcal_x_exp_hist ->GetEntries();
  
  // graphcanvas->cd(1); 
  // TLatex latex;
  // latex.SetNDC(true);
  // latex.DrawLatex(0.2, 0.85, Form("Scale Elastic = %.5f#pm %.5f", scale_elas, scale_elas_err));
  // latex.DrawLatex(0.2, 0.8, Form("#chi^{2}/ndf = %.2f/%.0f", ChiSq,ndf));
  // latex.DrawLatex(0.2, 0.75, Form("Data Events = %i", nEntries_data));
  // graphcanvas->Update();
  // upperPad1->Update();
  // lowerPad1->Update();

  TLine *zero_line = new TLine(x_exp_min, 0,x_exp_max , 0);
  zero_line ->SetLineColor(kRed);
  zero_line->SetLineWidth(2);

  TLegend *graphlegend = new TLegend(0.6, 0.6, 0.9, 0.9); // Adjust the legend coordinates as needed
  graphlegend->SetTextSize(0.05);  // Adjust this value as needed
  /// Add entries to the legend for the fit parameters
  graphlegend->AddEntry(h2_graph_x_exp, "H2", "P");
  graphlegend->AddEntry(result_graph_x_exp,"Scaled D2","P");
  graphlegend->Draw();

  // draw residual histogram on lower pad
  lowerPad1->cd();
  utilityHandler.AdjustHistLabelOffset(hist_resid_x_exp );
  hist_resid_x_exp  ->SetTitle("Residual");
  hist_resid_x_exp  ->Draw("E sames");
  zero_line ->Draw("same");

  graphcanvas->Update();
  graphcanvas->Draw();


  // canvas
  TCanvas *histcanvas = new TCanvas("histcanvas","histcanvas",800,600);
  //divide canvas 
  histcanvas->Divide(1,2);
  // Upper pad
  histcanvas->cd(1);
  TPad *upperPad = (TPad*)gPad;
  upperPad->SetGrid();
  upperPad->SetPad(0.01, 0.35, 0.99, 0.99); //Double_t xlow, Double_t ylow, Double_t xup, Double_t yup
  //utilityHandler.adjustPad(upperPad);
  //upperPad->SetFillColor(20);
  //upperPad->SetFrameFillColor(20);
    
  // Lower pad
  histcanvas->cd(2);
  TPad *lowerPad = (TPad*)gPad;
  lowerPad->SetGrid();
  lowerPad->SetPad(0.01, 0.01, 0.99, 0.34);
  //utilityHandler.adjustPad(lowerPad);
  //  lowerPad->SetFillColor(18);
  //lowerPad->SetFrameFillColor(18);

  // utilityHandler.AdjustHistLabelOffset(hist_result_x_exp);
  // utilityHandler.AdjustHistLabelOffset(h2_HDE_hcal_x_exp_hist);


  // Draw histogram and fit function on canvas
  upperPad->cd();
  h2_HDE_hcal_x_exp_hist->GetXaxis() ->SetRangeUser(x_exp_min, x_exp_max);
  h2_HDE_hcal_x_exp_hist->Draw("E");
  h2_HDE_hcal_x_exp_hist->SetLineWidth(2);
  //fit_result_x_exp->Draw("same");
  //fitFunc->Draw("same");
  hist_result_x_exp->SetLineColorAlpha(kGreen+2,1);
  hist_result_x_exp ->SetLineWidth(1);
  // hist_result_x_exp ->SetFillColorAlpha(kGreen,0.1);
  hist_result_x_exp ->Draw("E same");

  // Create a legend.
  TLegend *legend = new TLegend(0.55, 0.7, 0.9, 0.9); // Adjust the legend coordinates as needed
  legend->SetTextSize(0.05);  // Adjust this value as needed
  /// Add entries to the legend for the fit parameters
  legend->AddEntry(h2_HDE_hcal_x_exp_hist, "H2", "l");
  //legend->AddEntry(fit_result_x_exp,"Overall Fit","l");
  legend->AddEntry(hist_result_x_exp,"Scaled LD2","f l");
  //legend->AddEntry("", Form("R= %.6f +/- %.6f ", Rsf,Rsf_err), "");
  legend->AddEntry("", Form("#chi^{2}/ndf = %.2f", ChiSq_x_exp/ndf_x_exp), "");
  legend->Draw();

  // draw residual histogram on lower pad
  lowerPad->cd();
  //utilityHandler.AdjustHistLabelOffset(residual_hist);
  hist_resid_x_exp ->Draw("E sames");
  zero_line->Draw("same");

  upperPad->Update(); // this is supposed to update the axis labels but it looks like crap
  lowerPad->Update();
  histcanvas->Update();
  histcanvas ->SaveAs(Form("%s/Rsf_hist.pdf",outputfilelocation.Data() ));

  

  TH1D *d2_h2_ratio_x_exp = (TH1D*)d2_HDE_hcal_x_exp_hist->Clone("d2_h2_ratio_x_exp");
  d2_h2_ratio_x_exp->Divide(h2_HDE_hcal_x_exp_hist);/// hist1/hist2  == hist1->Divide(hist2); 
  d2_h2_ratio_x_exp->SetYTitle("HDE (d2) / HDE (h2)");
  d2_h2_ratio_x_exp->SetTitle("HDE (d2) / HDE (h2)");
  utilityHandler.RemoveExistingFit(d2_h2_ratio_x_exp);
  d2_h2_ratio_x_exp ->GetXaxis() ->SetRangeUser(x_exp_min, x_exp_max);
  // Adjusting histograms to make them look nice
  //utilityHandler.AdjustHistLabelOffset(d2_HDE_hcal_x_exp_hist);
  //utilityHandler.AdjustHistLabelOffset(h2_HDE_hcal_x_exp_hist);
  //utilityHandler.AdjustHistLabelOffset(d2_h2_ratio_x_exp);

  TCanvas* hcal_x_exp_canvas = new TCanvas("hcal_x_exp_canvas", "hcal_x_exp_canvas", 1000, 600);
  hcal_x_exp_canvas ->Divide(1,2);
  hcal_x_exp_canvas->cd(1);
  gPad->SetGrid();
  h2_HDE_hcal_x_exp_hist  ->SetTitle("hcal x expected");
  h2_HDE_hcal_x_exp_hist->SetLineColor(kMagenta);
  h2_HDE_hcal_x_exp_hist->Draw("E");
  d2_HDE_hcal_x_exp_hist  ->SetLineColor(kGreen+2);
  d2_HDE_hcal_x_exp_hist ->Draw("same E");
  TLegend *hcal_x_exp_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  hcal_x_exp_leg ->AddEntry(h2_HDE_hcal_x_exp_hist,"HDE on Hydrogen","l");
  hcal_x_exp_leg ->AddEntry(d2_HDE_hcal_x_exp_hist,"HDE on Deuterium","l");
  hcal_x_exp_leg->Draw();
  hcal_x_exp_canvas->cd(2);
  gPad->SetGrid();
  //utilityHandler.Fit_or_Replace_Pol0(d2_h2_ratio_x_exp, hcal_x_exp_range);
  d2_h2_ratio_x_exp->Draw("E");
  TLegend *hcal_x_exp_leg_r = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  hcal_x_exp_leg_r ->AddEntry(d2_h2_ratio_x_exp,"HDE Ratio: LD2/LH2","l");
  hcal_x_exp_leg_r->Draw();
  hcal_x_exp_canvas->Update();
  hcal_x_exp_canvas ->SaveAs(Form("%s/hcal_x_exp_HDE_ratio.pdf",outputfilelocation.Data() ) );
  hcal_x_exp_canvas->Update();

  
  // if (!d2_HDE_hcal_y_exp_hist || !h2_HDE_hcal_y_exp_hist_dont_draw)
  //   {
  //     cout<<"error, one of the histograms is null!"<<endl;
  //   }
  
  // fitHandlerGen_x_exp.Reset();
  // fitHandlerGen_x_exp.Set(d2_HDE_hcal_y_exp_hist, nullptr,y_exp_min+0.1, y_exp_max-0.1);

  // //fitHandlerGen_x_exp.fitData(h2_HDE_hcal_y_exp_hist_dont_draw);// bug in here that corrupts it and crashes if you try and draw

  // //TF1 *fit_result_y_exp = (TF1*)fitHandlerGen_x_exp.fitFunc->Clone("fit_result_y_exp");// bug in here that corrupts it and crashes if you try and draw

   
  /// Make /// Make a canvas called c1 to overwrite the default canvas and avoid crashes
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  d2_HDE_hcal_x_exp_hist->Draw("E");

 


  fout->Write();
  //f1->Close();
}// end main
