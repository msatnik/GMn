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
#include "TClass.h"
#include "TLatex.h"
#include "TKey.h"
#include "TFile.h"
#include "TStyle.h"
#include "TEnv.h"


#include "/work/halla//sbs/msatnik/GMn/classes/FitHistogramGeneral.h"
#include "/work/halla//sbs/msatnik/GMn/classes/FitHistogramGeneral.cpp"

#include "/work/halla//sbs/msatnik/GMn/classes/Utility.h"
#include "/work/halla//sbs/msatnik/GMn/classes/Utility.cpp"


#include "/work/halla//sbs/msatnik/GMn/classes/FileNames.h"
#include "/work/halla//sbs/msatnik/GMn/classes/FileNames.cpp"

//// default values to be overwritten 
TString DataFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_2Dhistos.root";
TString ProtonFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_2Dhistos.root";
TString NeutronFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_2Dhistos.root";
TString InelasticFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_inel_2Dhistos.root";
TString HistogramName = "W2_1d_allcuts";
std::string AxisTitle ="W^{2}";
double xmin = 0; 
double xmax =3;

// SBS4 30p 
TString DataFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_2Dhistos_sept26.root";
TString ProtonFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_2Dhistos_sept26.root";
TString NeutronFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_2Dhistos_sept26.root";
TString InelasticFileString_sbs4_30p  = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_inel_2Dhistos.root";
double xmin_sbs4_30p = 0; // 
double xmax_sbs4_30p =3;//

// SBS4 50p 
TString DataFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_2Dhistos.root";
TString ProtonFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deep_2Dhistos.root";
TString NeutronFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deen_2Dhistos.root";
TString InelasticFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_inel_2Dhistos.root";
double xmin_sbs4_50p = 0; // 
double xmax_sbs4_50p = 3;//

// SBS8 70p
TString DataFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_2Dhistos_sept11.root";
TString ProtonFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deep_2Dhistos.root";
TString NeutronFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deen_2Dhistos.root";
TString InelasticFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_inel_2Dhistos.root";
double xmin_sbs8_70p = 0; // 
double xmax_sbs8_70p = 3;//

// SBS9 70p
TString DataFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_2Dhistos_sept13.root";
TString ProtonFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_deep_2Dhistos.root";
TString NeutronFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_deen_2Dhistos.root";
TString InelasticFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_inel_2Dhistos.root";
double xmin_sbs9_70p = 0; // need to adjust range
double xmax_sbs9_70p = 3;//


void w2_elas_inel_compare(TString KineString="sbs4_30p"){ // main
  gStyle->SetNumberContours(255); 
  gStyle->SetOptStat(0110);
  gStyle->SetCanvasPreferGL(1);

  FileNames fileNamesHandler;

  if (KineString == "sbs4_30p"){
    DataFileString =fileNamesHandler.DataFileString_sbs4_30p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs4_30p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs4_30p;
    InelasticFileString = fileNamesHandler.InelasticFileString_sbs4_30p;
    xmin = xmin_sbs4_30p;
    xmax =xmax_sbs4_30p;
  }else if (KineString == "sbs4_50p"){
    DataFileString = fileNamesHandler.DataFileString_sbs4_50p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs4_50p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs4_50p;
    InelasticFileString =fileNamesHandler.InelasticFileString_sbs4_50p;
    xmin = xmin_sbs4_50p;
    xmax = xmax_sbs4_50p;
  }else if (KineString == "sbs8_70p"){
    DataFileString = fileNamesHandler.DataFileString_sbs8_70p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs8_70p;
    InelasticFileString = fileNamesHandler.InelasticFileString_sbs8_70p;
    xmin = xmin_sbs8_70p;
    xmax = xmax_sbs8_70p;
  }else if (KineString == "sbs9_70p"){
    DataFileString = fileNamesHandler.DataFileString_sbs9_70p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs9_70p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs9_70p;
    InelasticFileString = fileNamesHandler.InelasticFileString_sbs9_70p;
    xmin = xmin_sbs9_70p;
    xmax = xmax_sbs9_70p;
  }  else {
    std::cout<<"Error with kinematic setting"<<std::endl;
    return;
  }

  cout<<"What we're running: "<<KineString<<endl;
  cout<<DataFileString<<endl;
  cout<<ProtonFileString<<endl;
  cout<<NeutronFileString<<endl;
  cout<<InelasticFileString<<endl;
  cout<<xmin<< ", "<<xmax<<endl;
  
  TFile *f1 = TFile::Open(DataFileString); // data
  TFile *f2 = TFile::Open(ProtonFileString); //proton
  TFile *f3 = TFile::Open(NeutronFileString); //neutron
  TFile *f4 = TFile::Open(InelasticFileString); //neutron

  //// set location and name for root output file 
  TString outputfilelocation="/volatile/halla/sbs/msatnik/output/w2compare/"+KineString;
  TString outputfilename = outputfilelocation +"/"+ KineString+".root";

  TFile *fout = new TFile(outputfilename,"RECREATE");
  cout<<"writing to file: "<< outputfilename <<endl;

  // Load Histograms
  TH1D *hist_data_orig = (TH1D*)f1->Get(HistogramName);
  TH1D *hist_proton_orig = (TH1D*)f2->Get(HistogramName);
  TH1D *hist_neutron_orig = (TH1D*)f3->Get(HistogramName);
  TH1D *hist_inelastic_orig = (TH1D*)f4->Get(HistogramName);

  // Make clones of the histograms
  TH1D *hist_data = (TH1D*)hist_data_orig->Clone("hist_data");
  TH1D *hist_data_dont_draw = (TH1D*)hist_data_orig->Clone("hist_data_dont_draw");
  TH1D *hist_proton = (TH1D*)hist_proton_orig->Clone("hist_proton");
  TH1D *hist_neutron = (TH1D*)hist_neutron_orig->Clone("hist_neutron");
  TH1D *hist_inelastic = (TH1D*)hist_inelastic_orig->Clone("hist_inelastic");



  Utility utilityHandler; // class that gives us access to various functions to use between programs.


  // testing incrementRange here quick:
  double min=-1.5;
  double max = 1;
  int nDecimals = 1;
  double plus_or_minus_xmin = 0.02;
  double plus_or_minus_xmax = 0;
  double stepsize = 0.005;


  // std::string teststring;
  // std::string teststring2;
  // teststring = utilityHandler.incrementRangeStringStepsize(min, max,plus_or_minus_xmin, plus_or_minus_xmax, stepsize, nDecimals);

  // cout<<"Testing writing as a string with stepsize"<<endl;
  // cout<<teststring<<endl;
   
  //teststring2 = utilityHandler.incrementRangeStringPercents( min,  max, plus_or_minus_percent_xmin,plus_or_minus_percent_xmax, percent_increments,  nDecimals);
   
  // cout<<"Testing writing as a string with percent of range"<<endl;
  // cout<<teststring2<<endl;
  
 
  
  /// Make a histogram that is a sum of the proton and neutron histograms
  /// (not 100% sure this is what we want to do, but I'll try it for now
  TH1D *hist_elastic =  (TH1D*)hist_proton->Clone("hist_elastic");
  hist_elastic->Add(hist_neutron, 1); 

  //// Class to help fit the sim histograms to the data histogram 
  FitHistogramGeneral fitHandlerGen(hist_inelastic,hist_elastic,xmin,xmax);
  //FitHistogramGeneral fitHandlerGen(hist_elastic,hist_inelastic,xmin,xmax);

  fitHandlerGen.fitData(hist_data_dont_draw);// bug in here that corrupts it and crashes if you try and draw

  /// Make a canvas called c1 to overwrite the default canvas and avoid crashes
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  hist_elastic->Draw("E");

  TF1 *fit_result = (TF1*)fitHandlerGen.fitFunc->Clone("fit_result");// bug in here that corrupts it and crashes if you try and draw

  double scale_inel = fitHandlerGen.scale1;
  double scale_inel_err = fitHandlerGen.scale1_err;
  double scale_elas = fitHandlerGen.scale2;
  double scale_elas_err = fitHandlerGen.scale2_err;
  double ChiSq = fitHandlerGen.ChiSq;
  double ndf = fitHandlerGen.NDF;

  cout<< "Testing Classes: "<<endl;
  cout<<"scale_inelastic = "<< scale_inel <<" +/- "<<scale_inel_err<<endl;
  cout<<"scale_elastic = "<< scale_elas <<" +/- "<<scale_elas_err<<endl;

  TH1D* hist_inelastic_result = utilityHandler.ScaleAndShiftHistogram(hist_inelastic,scale_inel ,0);
  TH1D* hist_elastic_result = utilityHandler.ScaleAndShiftHistogram(hist_elastic, scale_elas,0);
  
  TCanvas *testc = new TCanvas("testc","testc",800,600);
  //hist_data->Draw("E");
  //hist_inelastic_result->Draw("same E");
  hist_elastic_result->Draw("same E");
  fit_result->Draw("same");

  TLine *zero_line = new TLine(xmin, 0, xmax, 0);
  zero_line ->SetLineColor(kRed);
  zero_line->SetLineWidth(2);

  double markersize =1;

  // Convert histogram to TGraphErrors
  TGraphErrors *data_graph = utilityHandler.histogramToGraphErrors(hist_data); 
  utilityHandler.customizeGraph(data_graph, 20, kBlack, markersize); // assigns shape, color, and marker size to the points
  data_graph->GetXaxis() ->SetRangeUser(xmin, xmax);

  TH1D *residual_hist = utilityHandler.GetResidualHistogram(hist_data, fit_result);
  residual_hist->GetXaxis() ->SetRangeUser(xmin, xmax);

  // Convert TF1 fit to TGraphErrors 
  TGraphErrors *fit_graph = utilityHandler.createGraphFromFit(hist_data,fit_result );
  utilityHandler.customizeGraph(fit_graph, 33, kRed, markersize); // assigns shape, color and marker size to the points
  fit_graph->GetXaxis() ->SetRangeUser(xmin, xmax);


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

  utilityHandler.AdjustHistLabelOffset(hist_elastic_result);
  utilityHandler.AdjustHistLabelOffset(hist_inelastic_result);
  utilityHandler.AdjustGraphLabelOffset(data_graph);

  upperPad1->cd();
  data_graph ->SetTitle("W^{2}");
  data_graph->Draw("AP");
  hist_elastic_result->SetLineColorAlpha(kGreen,0.9);
  hist_elastic_result->SetLineWidth(1);
  hist_elastic_result  ->SetFillColorAlpha(kGreen,0.1);
  hist_elastic_result ->Draw("hist E same");
  hist_inelastic_result ->SetLineColorAlpha(kMagenta,0.9);
  hist_inelastic_result ->SetLineWidth(1);
  hist_inelastic_result ->SetFillColorAlpha(kMagenta,0.1);
  hist_inelastic_result ->Draw("hist E same");
  data_graph->Draw("P SAME");
  data_graph->GetXaxis()->SetTitle("W^{2}");
  fit_graph->Draw("P SAME");

  int nEntries_data = hist_data->GetEntries();
  
  graphcanvas->cd(1); 
  TLatex latex;
  latex.SetNDC(true);
  latex.DrawLatex(0.15, 0.85, Form("Scale Elastic = %.5f#pm %.5f", scale_elas, scale_elas_err));
  latex.DrawLatex(0.15, 0.8, Form("#chi^{2}/ndf = %.2f/%.0f", ChiSq,ndf));
  latex.DrawLatex(0.15, 0.75, Form("Data Events = %i", nEntries_data));
  graphcanvas->Update();
  upperPad1->Update();
  lowerPad1->Update();
  

  TLegend *graphlegend = new TLegend(0.6, 0.6, 0.9, 0.9); // Adjust the legend coordinates as needed
   graphlegend->SetTextSize(0.05);  // Adjust this value as needed
  /// Add entries to the legend for the fit parameters
  graphlegend->AddEntry(data_graph, "Data", "P");
  graphlegend->AddEntry(fit_graph,"Overall Fit","P");
  graphlegend->AddEntry(hist_elastic_result,"Elastic simc (scaled)","f l");
  graphlegend->AddEntry(hist_inelastic_result,"Inelastic g4sbs (scaled)","f l");
  graphlegend->Draw();

  // draw residual histogram on lower pad
  lowerPad1->cd();
  utilityHandler.AdjustHistLabelOffset(residual_hist);
  residual_hist->SetTitle("Residual");
  residual_hist ->Draw("E sames");
  zero_line ->Draw("same");

  graphcanvas->Update();
  graphcanvas->Draw();

  
  
}// end main
