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
using namespace std;

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


void draw_SpotCutStudy(TString KineString="sbs4_30p"){// main

  gStyle->SetNumberContours(255); 
  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(1111);
  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0110);
  gStyle->SetStatTextColor(kBlack);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.06);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);


  Utility utilityHandler;

  if (KineString == "sbs4_30p"){
    DataFileString = DataFileString_sbs4_30p;
    W2_range=W2_range_sbs4_30p;
    hcal_x_exp_range=hcal_x_exp_range_sbs4_30p;
    hcal_y_exp_range = hcal_y_exp_range_sbs4_30p;
    hcal_dy_range = hcal_dy_range_sbs4_30p;
    hcal_e_min = hcal_e_min_sbs4_30p;
    coin_range = coin_range_sbs4_30p;
    ps_e_min = ps_e_min_sbs4_30p;
    ps_sh_e_min = ps_sh_e_min_sbs4_30p;
    grinch_clus_size_min = grinch_clus_size_min_sbs4_30p;
    grinch_clus_adc_min = grinch_clus_adc_min_sbs4_30p;
    e_over_p_range= e_over_p_range_sbs4_30p;
    vz_range = vz_range_sbs4_30p;
  }else if (KineString == "sbs4_50p"){
    DataFileString = DataFileString_sbs4_50p;
    W2_range=W2_range_sbs4_50p;
    hcal_x_exp_range=hcal_x_exp_range_sbs4_50p;
    hcal_y_exp_range = hcal_y_exp_range_sbs4_50p;
    hcal_dy_range = hcal_dy_range_sbs4_50p;
    hcal_e_min = hcal_e_min_sbs4_50p;
    coin_range = coin_range_sbs4_50p;
    ps_e_min = ps_e_min_sbs4_50p;
    ps_sh_e_min = ps_sh_e_min_sbs4_50p;
    grinch_clus_size_min = grinch_clus_size_min_sbs4_50p;
    grinch_clus_adc_min = grinch_clus_adc_min_sbs4_50p;
    e_over_p_range= e_over_p_range_sbs4_50p;
    vz_range = vz_range_sbs4_50p;
  }else if (KineString == "sbs8_70p"){
    DataFileString = DataFileString_sbs8_70p;
    W2_range=W2_range_sbs8_70p;
    hcal_x_exp_range=hcal_x_exp_range_sbs8_70p;
    hcal_y_exp_range = hcal_y_exp_range_sbs8_70p;
    hcal_dy_range = hcal_dy_range_sbs8_70p;
    hcal_e_min = hcal_e_min_sbs8_70p;
    coin_range = coin_range_sbs8_70p;
    ps_e_min = ps_e_min_sbs8_70p;
    ps_sh_e_min = ps_sh_e_min_sbs8_70p;
    grinch_clus_size_min = grinch_clus_size_min_sbs8_70p;
    grinch_clus_adc_min = grinch_clus_adc_min_sbs8_70p;
    e_over_p_range= e_over_p_range_sbs8_70p;
    vz_range = vz_range_sbs8_70p;
  }else if (KineString == "sbs9_70p"){
    DataFileString = DataFileString_sbs9_70p;
    W2_range=W2_range_sbs9_70p;
    hcal_x_exp_range=hcal_x_exp_range_sbs9_70p;
    hcal_y_exp_range = hcal_y_exp_range_sbs9_70p;
    hcal_dy_range = hcal_dy_range_sbs9_70p;
    hcal_e_min = hcal_e_min_sbs9_70p;
    coin_range = coin_range_sbs9_70p;
    ps_e_min = ps_e_min_sbs9_70p;
    ps_sh_e_min = ps_sh_e_min_sbs9_70p;
    grinch_clus_size_min = grinch_clus_size_min_sbs9_70p;
    grinch_clus_adc_min = grinch_clus_adc_min_sbs9_70p;
    e_over_p_range= e_over_p_range_sbs9_70p;
    vz_range = vz_range_sbs9_70p;
  }else {
    std::cout<<"Error with kinematic setting"<<std::endl;
    return;
  }
  
 
  cout<<endl;
  cout<<"What we're running: "<<KineString<<endl;
  cout<<DataFileString<<endl;
  cout<<endl;

 
  // load files 
  TFile *f1 = TFile::Open(DataFileString); // data

  if (!f1) {
    std::cout<<"Error with loading one of the rootfiles"<<std::endl;
    return;
  }

  // set location and name for output file 
  TString outputfilelocation="../../output/draw_SpotCutStudy/"+KineString;
  TString outputfilename = outputfilelocation +"/"+ KineString+".root";
  // Declare outfile
  TFile *fout = new TFile(outputfilename,"RECREATE");
  cout<<"writing to file: "<< outputfilename <<endl;

  // Load Histograms

  // W2
  TH1D *W2_hist = (TH1D*)f1->Get("W2_hist");
  TH1D *W2_hist_proton = (TH1D*)f1->Get("W2_hist_proton");
  TH1D *W2_hist_neutron = (TH1D*)f1->Get("W2_hist_neutron");
  TH1D *np_W2_hist = (TH1D*)f1->Get("np_W2_hist");

  // HCal X Expected
  TH1D *hcal_x_exp_hist = (TH1D*)f1->Get("hcal_x_exp_hist");
  TH1D *hcal_x_exp_hist_proton = (TH1D*)f1->Get("hcal_x_exp_hist_proton");
  TH1D *hcal_x_exp_hist_neutron = (TH1D*)f1->Get("hcal_x_exp_hist_neutron");
  TH1D *np_hcal_x_exp_hist = (TH1D*)f1->Get("np_hcal_x_exp_hist");
   
  // HCal Y Expected
  TH1D *hcal_y_exp_hist = (TH1D*)f1->Get("hcal_y_exp_hist");
  TH1D *hcal_y_exp_hist_proton = (TH1D*)f1->Get("hcal_y_exp_hist_proton");
  TH1D *hcal_y_exp_hist_neutron = (TH1D*)f1->Get("hcal_y_exp_hist_neutron");
  TH1D *np_hcal_y_exp_hist = (TH1D*)f1->Get("np_hcal_y_exp_hist");

  // HCal dy
  TH1D *hcal_dy_hist = (TH1D*)f1->Get("hcal_dy_hist");
  TH1D *hcal_dy_hist_proton = (TH1D*)f1->Get("hcal_dy_hist_proton");
  TH1D *hcal_dy_hist_neutron = (TH1D*)f1->Get("hcal_dy_hist_neutron");
  TH1D *np_hcal_dy_hist = (TH1D*)f1->Get("np_hcal_dy_hist");

  // HCal Energy 
  TH1D *hcal_e_hist = (TH1D*)f1->Get("hcal_e_hist");
  TH1D *hcal_e_hist_proton = (TH1D*)f1->Get("hcal_e_hist_proton");
  TH1D *hcal_e_hist_neutron = (TH1D*)f1->Get("hcal_e_hist_neutron");
  TH1D *np_hcal_e_hist = (TH1D*)f1->Get("np_hcal_e_hist");
  
  // HCal coincidence 
  TH1D *coin_hist = (TH1D*)f1->Get("coin_hist");
  TH1D *coin_hist_proton = (TH1D*)f1->Get("coin_hist_proton");
  TH1D *coin_hist_neutron = (TH1D*)f1->Get("coin_hist_neutron");
  TH1D *np_coin_hist = (TH1D*)f1->Get("np_coin_hist");

  // Preshower Energy 
  TH1D *ps_e_hist = (TH1D*)f1->Get("ps_e_hist");
  TH1D *ps_e_hist_proton = (TH1D*)f1->Get("ps_e_hist_proton");
  TH1D *ps_e_hist_neutron = (TH1D*)f1->Get("ps_e_hist_neutron");
  TH1D *np_ps_e_hist = (TH1D*)f1->Get("np_ps_e_hist");

  // Preshower + Shower Energy 
  TH1D *ps_sh_e_hist = (TH1D*)f1->Get("ps_sh_e_hist");
  TH1D *ps_sh_e_hist_proton = (TH1D*)f1->Get("ps_sh_e_hist_proton");
  TH1D *ps_sh_e_hist_neutron = (TH1D*)f1->Get("ps_sh_e_hist_neutron");
  TH1D *np_ps_sh_e_hist = (TH1D*)f1->Get("np_ps_sh_e_hist");

  // GRINCH cluster size  
  TH1D *grinch_clus_size_hist = (TH1D*)f1->Get("grinch_clus_size_hist");
  TH1D *grinch_clus_size_hist_proton = (TH1D*)f1->Get("grinch_clus_size_hist_proton");
  TH1D *grinch_clus_size_hist_neutron = (TH1D*)f1->Get("grinch_clus_size_hist_neutron");
  TH1D *np_grinch_clus_size_hist = (TH1D*)f1->Get("np_grinch_clus_size_hist");

  // GRINCH cluster adc (ToT sum)  
  TH1D *grinch_clus_adc_hist = (TH1D*)f1->Get("grinch_clus_adc_hist");
  TH1D *grinch_clus_adc_hist_proton = (TH1D*)f1->Get("grinch_clus_adc_hist_proton");
  TH1D *grinch_clus_adc_hist_neutron = (TH1D*)f1->Get("grinch_clus_adc_hist_neutron");
  TH1D *np_grinch_clus_adc_hist = (TH1D*)f1->Get("np_grinch_clus_adc_hist");

  // Total energy over track momentum 
  TH1D *e_over_p_hist = (TH1D*)f1->Get("e_over_p_hist");
  TH1D *e_over_p_hist_proton = (TH1D*)f1->Get("e_over_p_hist_proton");
  TH1D *e_over_p_hist_neutron = (TH1D*)f1->Get("e_over_p_hist_neutron");
  TH1D *np_e_over_p_hist = (TH1D*)f1->Get("np_e_over_p_hist");

  // track vertex 
  TH1D *vz_hist = (TH1D*)f1->Get("vz_hist");
  TH1D *vz_hist_proton = (TH1D*)f1->Get("vz_hist_proton");
  TH1D *vz_hist_neutron = (TH1D*)f1->Get("vz_hist_neutron");
  TH1D *np_vz_hist = (TH1D*)f1->Get("np_vz_hist");

  // HCal Hadron Detection Eff across hcal y expected
  TH1D *expected_hcal_y_exp_hist = (TH1D*)f1->Get("expected_hcal_y_exp_hist");
  TH1D *detected_hcal_y_exp_hist = (TH1D*)f1->Get("detected_hcal_y_exp_hist");
  TH1D *HDE_hcal_y_exp_hist = (TH1D*)f1->Get("HDE_hcal_y_exp_hist");

  // HCal Hadron Detection Eff across hcal x expected 
  TH1D *expected_hcal_x_exp_hist = (TH1D*)f1->Get("expected_hcal_x_exp_hist");
  TH1D *detected_hcal_x_exp_hist = (TH1D*)f1->Get("detected_hcal_x_exp_hist");
  TH1D *HDE_hcal_x_exp_hist = (TH1D*)f1->Get("HDE_hcal_x_exp_hist");

  
  // Drawing W2
  // Adjusting histograms to make them look nice
  utilityHandler.AdjustHistLabelOffset(W2_hist);
  utilityHandler.AdjustHistLabelOffset(W2_hist_proton);
  utilityHandler.AdjustHistLabelOffset(W2_hist_neutron);
  utilityHandler.AdjustHistLabelOffset(np_W2_hist);

  TCanvas* W2_canvas = new TCanvas("W2_canvas", "W2_canvas", 1000, 600);
  W2_canvas ->Divide(1,2);
  W2_canvas->cd(1);
  W2_hist->SetTitle("W^{2}");
  W2_hist->Draw("E");
  W2_hist_proton->SetLineColor(kGreen);
  W2_hist_proton->Draw("same E");
  W2_hist_neutron->SetLineColor(kMagenta);
  W2_hist_neutron->Draw("same E");
  W2_canvas->Update();
  TLegend *W2_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  W2_leg ->AddEntry(W2_hist,"W^{2} w/ no spot cuts","l");
  W2_leg ->AddEntry(W2_hist_proton,"W^{2} w/ proton spot cut","l");
  W2_leg ->AddEntry(W2_hist_neutron,"W^{2} w/ neutron spot cut","l");
  W2_leg->Draw();
  W2_canvas->cd(2);
  utilityHandler.Fit_or_Replace_Pol0(np_W2_hist, W2_range);
  np_W2_hist->Draw("E");
  W2_canvas->Update();
  W2_canvas ->SaveAs(Form("%s/W2_np.pdf",outputfilelocation.Data() ) );
  W2_canvas->Update();

  // Drawing HCal X expected
  // Adjusting histograms to make them look nice
  utilityHandler.AdjustHistLabelOffset(hcal_x_exp_hist);
  utilityHandler.AdjustHistLabelOffset(hcal_x_exp_hist_proton);
  utilityHandler.AdjustHistLabelOffset(hcal_x_exp_hist_neutron);
  utilityHandler.AdjustHistLabelOffset(np_hcal_x_exp_hist);
  TCanvas* hcal_x_exp_canvas = new TCanvas("hcal_x_exp_canvas", "hcal_x_exp_canvas", 1000, 600);
  hcal_x_exp_canvas ->Divide(1,2);
  hcal_x_exp_canvas->cd(1);
  gPad->SetLogy();
  hcal_x_exp_hist->SetTitle("hcal x expected");
  hcal_x_exp_hist->Draw("E");
  hcal_x_exp_hist_proton->SetLineColor(kGreen);
  hcal_x_exp_hist_proton->Draw("same E");
  hcal_x_exp_hist_neutron->SetLineColor(kMagenta);
  hcal_x_exp_hist_neutron->Draw("same E");
  TLegend *hcal_x_exp_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  hcal_x_exp_leg ->AddEntry(hcal_x_exp_hist,"No Spot Cuts","l");
  hcal_x_exp_leg ->AddEntry(hcal_x_exp_hist_proton,"Proton Spot Cut","l");
  hcal_x_exp_leg ->AddEntry(hcal_x_exp_hist_neutron,"Neutron Spot Cut","l");
  hcal_x_exp_leg->Draw();
  hcal_x_exp_canvas->cd(2);
  utilityHandler.Fit_or_Replace_Pol0(np_hcal_x_exp_hist, hcal_x_exp_range);
  np_hcal_x_exp_hist->Draw("E");
  hcal_x_exp_canvas->Update();
  hcal_x_exp_canvas ->SaveAs(Form("%s/hcal_x_exp_np.pdf",outputfilelocation.Data() ) );
  hcal_x_exp_canvas->Update();

  
  // Drawing HCal Y expected
  // Adjusting histograms to make them look nice
  utilityHandler.AdjustHistLabelOffset(hcal_y_exp_hist);
  utilityHandler.AdjustHistLabelOffset(hcal_y_exp_hist_proton);
  utilityHandler.AdjustHistLabelOffset(hcal_y_exp_hist_neutron);
  utilityHandler.AdjustHistLabelOffset(np_hcal_y_exp_hist);
  TCanvas* hcal_y_exp_canvas = new TCanvas("hcal_y_exp_canvas", "hcal_y_exp_canvas", 1000, 600);
  hcal_y_exp_canvas ->Divide(1,2);
  hcal_y_exp_canvas->cd(1);
  gPad->SetLogy();
  hcal_y_exp_hist->SetTitle("hcal y expected");
  hcal_y_exp_hist->Draw("E");
  hcal_y_exp_hist_proton->SetLineColor(kGreen);
  hcal_y_exp_hist_proton->Draw("same E");
  hcal_y_exp_hist_neutron->SetLineColor(kMagenta);
  hcal_y_exp_hist_neutron->Draw("same E");
  TLegend *hcal_y_exp_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  hcal_y_exp_leg ->AddEntry(hcal_y_exp_hist,"No Spot Cuts","l");
  hcal_y_exp_leg ->AddEntry(hcal_y_exp_hist_proton,"Proton Spot Cut","l");
  hcal_y_exp_leg ->AddEntry(hcal_y_exp_hist_neutron,"Neutron Spot Cut","l");
  hcal_y_exp_leg->Draw();
  hcal_y_exp_canvas->cd(2);
  // Replacing fit 
  utilityHandler.Fit_or_Replace_Pol0(np_hcal_y_exp_hist, hcal_y_exp_range);
  np_hcal_y_exp_hist ->GetYaxis() ->SetRangeUser(0, 1);
  np_hcal_y_exp_hist->Draw("E");
  hcal_y_exp_canvas->Update();
  hcal_y_exp_canvas ->SaveAs(Form("%s/hcal_y_exp_np.pdf",outputfilelocation.Data() ) );
  

  // Drawing HCal dy
  // Adjusting histograms to make them look nice
  utilityHandler.AdjustHistLabelOffset(hcal_dy_hist);
  utilityHandler.AdjustHistLabelOffset(hcal_dy_hist_proton);
  utilityHandler.AdjustHistLabelOffset(hcal_dy_hist_neutron);
  utilityHandler.AdjustHistLabelOffset(np_hcal_dy_hist);
  TCanvas* hcal_dy_canvas = new TCanvas("hcal_dy_canvas", "hcal_dy_canvas", 1000, 600);
  hcal_dy_canvas ->Divide(1,2);
  hcal_dy_canvas->cd(1);
  gPad->SetLogy();
  hcal_dy_hist->SetTitle("hcal dy");
  hcal_dy_hist->Draw("E");
  hcal_dy_hist_proton->SetLineColor(kGreen);
  hcal_dy_hist_proton->Draw("same E");
  hcal_dy_hist_neutron->SetLineColor(kMagenta);
  hcal_dy_hist_neutron->Draw("same E");
  TLegend *hcal_dy_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  hcal_dy_leg ->AddEntry(hcal_dy_hist,"No Spot Cuts","l");
  hcal_dy_leg ->AddEntry(hcal_dy_hist_proton,"Proton Spot Cut","l");
  hcal_dy_leg ->AddEntry(hcal_dy_hist_neutron,"Neutron Spot Cut","l");
  hcal_dy_leg->Draw();
  hcal_dy_canvas->cd(2);
  // Replacing fit 
  utilityHandler.Fit_or_Replace_Pol0(np_hcal_dy_hist, hcal_dy_range);
  np_hcal_dy_hist->Draw("E");
  hcal_dy_canvas->Update();
  hcal_dy_canvas ->SaveAs(Form("%s/hcal_dy_np.pdf",outputfilelocation.Data() ) );

  // Drawing HCal Energy
  // Adjusting histograms to make them look nice
  utilityHandler.AdjustHistLabelOffset(hcal_e_hist);
  utilityHandler.AdjustHistLabelOffset(hcal_e_hist_proton);
  utilityHandler.AdjustHistLabelOffset(hcal_e_hist_neutron);
  utilityHandler.AdjustHistLabelOffset(np_hcal_e_hist);
  TCanvas* hcal_e_canvas = new TCanvas("hcal_e_canvas", "hcal_e_canvas", 1000, 600);
  hcal_e_canvas ->Divide(1,2);
  hcal_e_canvas->cd(1);
  hcal_e_hist->SetTitle("HCal Energy");
  hcal_e_hist->Draw("E");
  hcal_e_hist_proton->SetLineColor(kGreen);
  hcal_e_hist_proton->Draw("same E");
  hcal_e_hist_neutron->SetLineColor(kMagenta);
  hcal_e_hist_neutron->Draw("same E");
  TLegend *hcal_e_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  hcal_e_leg ->AddEntry(hcal_e_hist,"No Spot Cuts","l");
  hcal_e_leg ->AddEntry(hcal_e_hist_proton,"Proton Spot Cut","l");
  hcal_e_leg ->AddEntry(hcal_e_hist_neutron,"Neutron Spot Cut","l");
  hcal_e_leg->Draw();
  hcal_e_canvas->cd(2);
  // Replacing fit 
  utilityHandler.Fit_or_Replace_Pol0_xmin(np_hcal_e_hist, hcal_e_min);
  np_hcal_e_hist->Draw("E");
  hcal_e_canvas->Update();
  hcal_e_canvas ->SaveAs(Form("%s/hcal_e_np.pdf",outputfilelocation.Data() ) );

 
  // Draw HCal coincidence
  // Adjusting histograms to make them look nice
  utilityHandler.AdjustHistLabelOffset(coin_hist);
  utilityHandler.AdjustHistLabelOffset(coin_hist_proton);
  utilityHandler.AdjustHistLabelOffset(coin_hist_neutron);
  utilityHandler.AdjustHistLabelOffset(np_coin_hist);
  TCanvas* coin_canvas = new TCanvas("coin_canvas", "coin_canvas", 1000, 600);
  coin_canvas ->Divide(1,2);
  coin_canvas->cd(1);
  coin_hist->SetTitle("HCal - SH time");
  coin_hist->Draw("E");
  coin_hist_proton->SetLineColor(kGreen);
  coin_hist_proton->Draw("same E");
  coin_hist_neutron->SetLineColor(kMagenta);
  coin_hist_neutron->Draw("same E");
  TLegend *coin_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  coin_leg ->AddEntry(coin_hist,"No Spot Cuts","l");
  coin_leg ->AddEntry(coin_hist_proton,"Proton Spot Cut","l");
  coin_leg ->AddEntry(coin_hist_neutron,"Neutron Spot Cut","l");
  coin_leg->Draw();
  coin_canvas->cd(2);
  // Replacing fit 
  utilityHandler.Fit_or_Replace_Pol0(np_coin_hist, coin_range);
  np_coin_hist->GetYaxis()->SetRangeUser(0,3);
  np_coin_hist->Draw("E");
  coin_canvas->Update();
  coin_canvas ->SaveAs(Form("%s/coin_np.pdf",outputfilelocation.Data() ) );

  // Drawing preshower energy
  // Adjusting histograms to make them look nice
  utilityHandler.AdjustHistLabelOffset(ps_e_hist);
  utilityHandler.AdjustHistLabelOffset(ps_e_hist_proton);
  utilityHandler.AdjustHistLabelOffset(ps_e_hist_neutron);
  utilityHandler.AdjustHistLabelOffset(np_ps_e_hist);
  TCanvas* ps_e_canvas = new TCanvas("ps_e_canvas", "ps_e_canvas", 1000, 600);
  ps_e_canvas ->Divide(1,2);
  ps_e_canvas->cd(1);
  ps_e_hist->SetTitle("Preshower Energy");
  ps_e_hist->Draw("E");
  ps_e_hist_proton->SetLineColor(kGreen);
  ps_e_hist_proton->Draw("same E");
  ps_e_hist_neutron->SetLineColor(kMagenta);
  ps_e_hist_neutron->Draw("same E");
  TLegend *ps_e_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  ps_e_leg ->AddEntry(ps_e_hist,"No Spot Cuts","l");
  ps_e_leg ->AddEntry(ps_e_hist_proton,"Proton Spot Cut","l");
  ps_e_leg ->AddEntry(ps_e_hist_neutron,"Neutron Spot Cut","l");
  ps_e_leg->Draw();
  ps_e_canvas->cd(2);
  // Replacing fit 
  utilityHandler.Fit_or_Replace_Pol0_xmin(np_ps_e_hist, ps_e_min);
  np_ps_e_hist->Draw("E");
  ps_e_canvas->Update();
  ps_e_canvas ->SaveAs(Form("%s/ps_e_np.pdf",outputfilelocation.Data() ) ); 
  
  // Drawing total energy
  // Adjusting histograms to make them look nice
  utilityHandler.AdjustHistLabelOffset(ps_sh_e_hist);
  utilityHandler.AdjustHistLabelOffset(ps_sh_e_hist_proton);
  utilityHandler.AdjustHistLabelOffset(ps_sh_e_hist_neutron);
  utilityHandler.AdjustHistLabelOffset(np_ps_sh_e_hist);
  TCanvas* ps_sh_e_canvas = new TCanvas("ps_sh_e_canvas", "ps_sh_e_canvas", 1000, 600);
  ps_sh_e_canvas ->Divide(1,2);
  ps_sh_e_canvas->cd(1);
  ps_sh_e_hist->SetTitle("Preshower + Shower Energy");
  ps_sh_e_hist->Draw("E");
  ps_sh_e_hist_proton->SetLineColor(kGreen);
  ps_sh_e_hist_proton->Draw("same E");
  ps_sh_e_hist_neutron->SetLineColor(kMagenta);
  ps_sh_e_hist_neutron->Draw("same E");
  TLegend *ps_sh_e_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  ps_sh_e_leg ->AddEntry(ps_sh_e_hist,"No Spot Cuts","l");
  ps_sh_e_leg ->AddEntry(ps_sh_e_hist_proton,"Proton Spot Cut","l");
  ps_sh_e_leg ->AddEntry(ps_sh_e_hist_neutron,"Neutron Spot Cut","l");
  ps_sh_e_leg->Draw();
  ps_sh_e_canvas->cd(2);
  // Replacing fit 
  utilityHandler.Fit_or_Replace_Pol0_xmin(np_ps_sh_e_hist, ps_sh_e_min);
  np_ps_sh_e_hist->Draw("E");
  ps_sh_e_canvas->Update();
  ps_sh_e_canvas ->SaveAs(Form("%s/ps_sh_e_np.pdf",outputfilelocation.Data() ) ); 		
		       

  // Drawing grinch cluster size
  // Adjusting histograms to make them look nice
  utilityHandler.AdjustHistLabelOffset(grinch_clus_size_hist);
  utilityHandler.AdjustHistLabelOffset(grinch_clus_size_hist_proton);
  utilityHandler.AdjustHistLabelOffset(grinch_clus_size_hist_neutron);
  utilityHandler.AdjustHistLabelOffset(np_grinch_clus_size_hist);
  TCanvas* grinch_clus_size_canvas = new TCanvas("grinch_clus_size_canvas", "grinch_clus_size_canvas", 1000, 600);
  grinch_clus_size_canvas ->Divide(1,2);
  grinch_clus_size_canvas->cd(1);
  grinch_clus_size_hist->SetTitle("GRINCH cluster size");
  grinch_clus_size_hist->Draw("E");
  grinch_clus_size_hist_proton->SetLineColor(kGreen);
  grinch_clus_size_hist_proton->Draw("same E");
  grinch_clus_size_hist_neutron->SetLineColor(kMagenta);
  grinch_clus_size_hist_neutron->Draw("same E");
  TLegend *grinch_clus_size_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  grinch_clus_size_leg ->AddEntry(grinch_clus_size_hist,"No Spot Cuts","l");
  grinch_clus_size_leg ->AddEntry(grinch_clus_size_hist_proton,"Proton Spot Cut","l");
  grinch_clus_size_leg ->AddEntry(grinch_clus_size_hist_neutron,"Neutron Spot Cut","l");
  grinch_clus_size_leg->Draw();
  grinch_clus_size_canvas->cd(2);
  // Replacing fit 
  utilityHandler.Fit_or_Replace_Pol0_xmin(np_grinch_clus_size_hist, grinch_clus_size_min);
  np_grinch_clus_size_hist->Draw("E");
  grinch_clus_size_canvas->Update();
  grinch_clus_size_canvas ->SaveAs(Form("%s/grinch_clus_size_np.pdf",outputfilelocation.Data() ) ); 	


  // Drawing GRINCH cluster adc (total time-over-threshold)
  // Adjusting histograms to make them look nice
  utilityHandler.AdjustHistLabelOffset(grinch_clus_adc_hist);
  utilityHandler.AdjustHistLabelOffset(grinch_clus_adc_hist_proton);
  utilityHandler.AdjustHistLabelOffset(grinch_clus_adc_hist_neutron);
  utilityHandler.AdjustHistLabelOffset(np_grinch_clus_adc_hist);
  TCanvas* grinch_clus_adc_canvas = new TCanvas("grinch_clus_adc_canvas", "grinch_clus_adc_canvas", 1000, 600);
  grinch_clus_adc_canvas ->Divide(1,2);
  grinch_clus_adc_canvas->cd(1);
  grinch_clus_adc_hist->SetTitle("GRINCH cluster ToT");
  grinch_clus_adc_hist->Draw("E");
  grinch_clus_adc_hist_proton->SetLineColor(kGreen);
  grinch_clus_adc_hist_proton->Draw("same E");
  grinch_clus_adc_hist_neutron->SetLineColor(kMagenta);
  grinch_clus_adc_hist_neutron->Draw("same E");
  TLegend *grinch_clus_adc_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  grinch_clus_adc_leg ->AddEntry(grinch_clus_adc_hist,"No Spot Cuts","l");
  grinch_clus_adc_leg ->AddEntry(grinch_clus_adc_hist_proton,"Proton Spot Cut","l");
  grinch_clus_adc_leg ->AddEntry(grinch_clus_adc_hist_neutron,"Neutron Spot Cut","l");
  grinch_clus_adc_leg->Draw();
  grinch_clus_adc_canvas->cd(2);
  // Replacing fit 
  utilityHandler.Fit_or_Replace_Pol0_xmin(np_grinch_clus_adc_hist, grinch_clus_adc_min);
  np_grinch_clus_adc_hist->Draw("E");
  grinch_clus_adc_canvas->Update();
  grinch_clus_adc_canvas ->SaveAs(Form("%s/grinch_clus_adc_np.pdf",outputfilelocation.Data() ) ); 	



  // Drawing total energy over momentum
  // Adjusting histograms to make them look nice
  utilityHandler.AdjustHistLabelOffset(e_over_p_hist);
  utilityHandler.AdjustHistLabelOffset(e_over_p_hist_proton);
  utilityHandler.AdjustHistLabelOffset(e_over_p_hist_neutron);
  utilityHandler.AdjustHistLabelOffset(np_e_over_p_hist);
  TCanvas* e_over_p_canvas = new TCanvas("e_over_p_canvas", "e_over_p_canvas", 1000, 600);
  e_over_p_canvas ->Divide(1,2);
  e_over_p_canvas->cd(1);
  e_over_p_hist->SetTitle("E over P");
  e_over_p_hist->Draw("E");
  e_over_p_hist_proton->SetLineColor(kGreen);
  e_over_p_hist_proton->Draw("same E");
  e_over_p_hist_neutron->SetLineColor(kMagenta);
  e_over_p_hist_neutron->Draw("same E");
  TLegend *e_over_p_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  e_over_p_leg ->AddEntry(e_over_p_hist,"No Spot Cuts","l");
  e_over_p_leg ->AddEntry(e_over_p_hist_proton,"Proton Spot Cut","l");
  e_over_p_leg ->AddEntry(e_over_p_hist_neutron,"Neutron Spot Cut","l");
  e_over_p_leg->Draw();
  e_over_p_canvas->cd(2);
  // Replacing fit 
  utilityHandler.Fit_or_Replace_Pol0(np_e_over_p_hist, e_over_p_range);
  np_e_over_p_hist->Draw("E");
  e_over_p_canvas->Update();
  e_over_p_canvas ->SaveAs(Form("%s/e_over_p_np.pdf",outputfilelocation.Data() ) );


  // Drawing vz 
  TCanvas* vz_canvas = new TCanvas("vz_canvas", "vz_canvas", 1000, 600);
  // Adjusting histograms to make them look nice
  utilityHandler.AdjustHistLabelOffset(vz_hist);
  utilityHandler.AdjustHistLabelOffset(vz_hist_proton);
  utilityHandler.AdjustHistLabelOffset(vz_hist_neutron);
  utilityHandler.AdjustHistLabelOffset(np_vz_hist);
  vz_canvas ->Divide(1,2);
  vz_canvas->cd(1);
  vz_hist->SetTitle("Track Vertex vz");
  vz_hist->Draw("E");
  vz_hist_proton->SetLineColor(kGreen);
  vz_hist_proton->Draw("same E");
  vz_hist_neutron->SetLineColor(kMagenta);
  vz_hist_neutron->Draw("same E");
  TLegend *vz_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  vz_leg ->AddEntry(vz_hist,"No Spot Cuts","l");
  vz_leg ->AddEntry(vz_hist_proton,"Proton Spot Cut","l");
  vz_leg ->AddEntry(vz_hist_neutron,"Neutron Spot Cut","l");
  vz_leg->Draw();
  vz_canvas->cd(2);
  // Replacing fit 
  utilityHandler.Fit_or_Replace_Pol0(np_vz_hist, vz_range);
  np_vz_hist->Draw("E");
  vz_canvas->Update();
  vz_canvas ->SaveAs(Form("%s/vz_np.pdf",outputfilelocation.Data() ) );

  // Drawing HCal Detection Eff across hcal y expected
  // Adjusting histograms to make them look nice
  utilityHandler.AdjustHistLabelOffset(expected_hcal_y_exp_hist);
  utilityHandler.AdjustHistLabelOffset(detected_hcal_y_exp_hist);
  utilityHandler.AdjustHistLabelOffset(HDE_hcal_y_exp_hist);
  TCanvas* HDE_hcal_y_exp_canvas = new TCanvas("HDE_hcal_y_exp_canvas", "HDE_hcal_y_exp_canvas", 1000, 600);
  HDE_hcal_y_exp_canvas ->Divide(1,2);
  HDE_hcal_y_exp_canvas->cd(1);
  gPad->SetLogy();
  expected_hcal_y_exp_hist->SetTitle("hcal y expected");
  expected_hcal_y_exp_hist->Draw("E");
  detected_hcal_y_exp_hist->SetLineColor(kViolet);
  detected_hcal_y_exp_hist->Draw("same E");
  TLegend *HDE_hcal_y_exp_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  HDE_hcal_y_exp_leg ->AddEntry(expected_hcal_y_exp_hist,"No HCal Cuts","l");
  HDE_hcal_y_exp_leg ->AddEntry(detected_hcal_y_exp_hist,"HCal Cuts and Spot Cuts","l");
  HDE_hcal_y_exp_leg->Draw();
  HDE_hcal_y_exp_canvas->cd(2);
  // Replacing fit 
  utilityHandler.Fit_or_Replace_Pol0(HDE_hcal_y_exp_hist, hcal_y_exp_range);
  HDE_hcal_y_exp_hist->Draw("E");
  HDE_hcal_y_exp_canvas->Update();
  HDE_hcal_y_exp_canvas ->SaveAs(Form("%s/HDE_hcal_y_exp.pdf",outputfilelocation.Data() ) );

  // Drawing HCal Detection Eff across hcal x expeted
  utilityHandler.AdjustHistLabelOffset(expected_hcal_x_exp_hist);
  utilityHandler.AdjustHistLabelOffset(detected_hcal_x_exp_hist);
  utilityHandler.AdjustHistLabelOffset(HDE_hcal_x_exp_hist);
  TCanvas* HDE_hcal_x_exp_canvas = new TCanvas("HDE_hcal_x_exp_canvas", "HDE_hcal_x_exp_canvas", 1000, 600);
  HDE_hcal_x_exp_canvas ->Divide(1,2);
  HDE_hcal_x_exp_canvas->cd(1);
  gPad->SetLogy();
  expected_hcal_x_exp_hist->SetTitle("hcal x expected");
  expected_hcal_x_exp_hist->Draw("E");
  detected_hcal_x_exp_hist->SetLineColor(kViolet);
  detected_hcal_x_exp_hist->Draw("same E");
  TLegend *HDE_hcal_x_exp_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  HDE_hcal_x_exp_leg ->AddEntry(expected_hcal_x_exp_hist,"No HCal Cuts","l");
  HDE_hcal_x_exp_leg ->AddEntry(detected_hcal_x_exp_hist,"HCal Cuts and Spot Cuts","l");
  HDE_hcal_x_exp_leg->Draw();
  HDE_hcal_x_exp_canvas->cd(2);
  // Replacing fit 
  utilityHandler.Fit_or_Replace_Pol0(HDE_hcal_x_exp_hist, hcal_x_exp_range);
  HDE_hcal_x_exp_hist->Draw("E");
  HDE_hcal_x_exp_canvas->Update();
  HDE_hcal_x_exp_canvas ->SaveAs(Form("%s/HDE_hcal_x_exp.pdf",outputfilelocation.Data() ) );


  fout->Write();
  //f1->Close();
}// end main
