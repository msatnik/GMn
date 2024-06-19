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

// Global Params

// Function Declarations

// MAIN
void plot_dxdy_data(){

  // Load rootfiles for histograms
  TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p.root"); 
  //TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_sim_sf03_fid.root"); 

  // Load Histograms
  TH1D *hist1 = (TH1D*)f1->Get("h_dx_HCAL_noAA");
  TH1D *hist2 = (TH1D*)f1->Get("h_dx_HCAL");
  TH1D *hist3 = (TH1D*)f1->Get("h_dx_HCAL_failedFid");

  TH1D *hist4 = (TH1D*)f1->Get("h_dy_HCAL_noAA");
  TH1D *hist5 = (TH1D*)f1->Get("h_dy_HCAL");
  TH1D *hist6 = (TH1D*)f1->Get("h_dy_HCAL_failedFid");

  //// Normalize the integrals of the histograms to 1
  //hist1->Scale(1.0 / hist1->Integral());
  // hist2->Scale(1.0 / hist2->Integral());

  

  // Create a THStack which we can use to automatically set the scale when we draw the histos
  THStack *stack = new THStack("stack", "");
  // Add histograms to the stack
  stack->Add(hist1);
  stack->Add(hist2);
  stack->Add(hist3);

  THStack *stack2 = new THStack("stack2", "");
  // Add histograms to the stack
  stack2->Add(hist4);
  stack2->Add(hist5);
  stack2->Add(hist6);
 
  // Set up canvas 
  TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
  c1->Divide(2,1);
  c1->cd(1);
  hist1->SetLineColor(kBlue-3);
  //hist1->SetLineWidth(2);
  hist2->SetLineColor(kGreen-3);
  //hist2 ->SetLineWidth(2);
  hist3->SetLineColor(kRed-3);
  //hist3->SetLineWidth(2);
 
  // Draw the stacked histograms
  stack->Draw("nostack"); // "nostack" option prevents histograms from being stacked on top of each other. y-axis is scaled automatically this way. 
  stack->SetTitle("HCal dx");

  // Create a legend.
  TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9); // Adjust the legend coordinates as needed
  // Add entries to the legend for the fit parameters
  legend->AddEntry(hist1, "No Fid Cuts", "l");
  legend->AddEntry(hist2, "Passed Fid Cuts", "l");
  legend ->AddEntry(hist3,"Failed Fid Cuts","l");
  legend->Draw();
 

 c1->cd(2);
  hist4->SetLineColor(kBlue-3);
  //hist4->SetLineWidth(2);
  hist5->SetLineColor(kGreen-3);
  //hist5 ->SetLineWidth(2);
  hist6->SetLineColor(kRed-3);
  // hist6->SetLineWidth(2);

   // Draw the stacked histograms
  stack2->Draw("nostack"); // "nostack" option prevents histograms from being stacked on top of each other. y-axis is scaled automatically this way. 
  stack2->SetTitle("HCal dy");

   // Create a legend.
  TLegend *legend2 = new TLegend(0.6, 0.7, 0.9, 0.9); // Adjust the legend coordinates as needed
  // Add entries to the legend for the fit parameters
  legend2->AddEntry(hist4, "No Fid Cuts", "l");
  legend2->AddEntry(hist5, "Passed Fid Cuts", "l");
  legend2 ->AddEntry(hist6,"Failed Fid Cuts","l");
  legend2->Draw();
  
  c1->Update();
}
