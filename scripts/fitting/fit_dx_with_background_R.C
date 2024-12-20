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

#include "/work/halla//sbs/msatnik/GMn/classes/HelloWorld.h"
#include "/work/halla//sbs/msatnik/GMn/classes/HelloWorld.cpp"

#include "/work/halla//sbs/msatnik/GMn/classes/FitHistogram.h"
#include "/work/halla//sbs/msatnik/GMn/classes/FitHistogram.cpp"

#include "/work/halla//sbs/msatnik/GMn/classes/Utility.h"
#include "/work/halla//sbs/msatnik/GMn/classes/Utility.cpp"


#include "/work/halla//sbs/msatnik/GMn/classes/FileNames.h"
#include "/work/halla//sbs/msatnik/GMn/classes/FileNames.cpp"

// Global Params


//// default values to be overwritten 
TString DataFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_tight_2Dhistos.root";
TString ProtonFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_tight_2Dhistos.root";
TString NeutronFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_tight_2Dhistos.root";
TString InelasticFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_inel_2Dhistos.root";
TString BackgroundFileString = DataFileString;
double xmin = -2.1; // -2.1 for SBS4 30p, maybe -2.5 for 50p
double xmax = 1.4;//1.4 for SBS4 40p
double xminBG = -2.4;
double xmaxBG = 1.4;
int polyorder=5;
TString BackgroundHistogramName = "hcal_dx_1d_allcuts";
TString config_title = "Testing out";


// SBS4 30p 
// TString DataFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_2Dhistos_sept26.root";
// TString ProtonFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_2Dhistos_sept26.root";
// TString NeutronFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_2Dhistos_sept26.root";
// TString InelasticFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_inel_2Dhistos.root";
// double xmin_sbs4_30p = -2.15; 
// double xmax_sbs4_30p = 1.4;
// double xminBG_sbs4_30p = -2.15;
// double xmaxBG_sbs4_30p = 1.4;
TString config_title_sbs4_30p = "SBS4 30%";

// SBS4 50p 
// TString DataFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_2Dhistos_sept26.root";
// TString ProtonFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deep_2Dhistos_sept26.root";
// TString NeutronFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deen_2Dhistos_sept26.root";
// TString InelasticFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_inel_2Dhistos.root";
// double xmin_sbs4_50p = -2.63; 
// double xmax_sbs4_50p = 1.4;
// double xminBG_sbs4_50p = -2.63;
// double xmaxBG_sbs4_50p = 1.4;
TString config_title_sbs4_50p = "SBS4 50%";

// SBS8 70p
// TString DataFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_2Dhistos_sept26.root";
// TString ProtonFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deep_2Dhistos_sept26.root";
// TString NeutronFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deen_2Dhistos_sept26.root";
// TString InelasticFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_inel_2Dhistos.root";
double xmin_sbs8_70p = -1.9; 
double xmax_sbs8_70p = 1.0;
double xminBG_sbs8_70p = -1.9;
double xmaxBG_sbs8_70p = 1.0;
TString config_title_sbs8_70p = "SBS8 70%";


// // SBS8 50p
// TString DataFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_2Dhistos_sept26.root";
// TString ProtonFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_deep_2Dhistos_sept26.root";
// TString NeutronFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_deen_2Dhistos_sept26.root";
// TString InelasticFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_inel_2Dhistos.root";// NEED TO RUN INEL
// double xmin_sbs8_50p = -1.6; 
// double xmax_sbs8_50p = 1.0;
// double xminBG_sbs8_50p = -1.6;
// double xmaxBG_sbs8_50p = 1.0;
TString config_title_sbs8_50p = "SBS8 50%";

// // SBS8 100p
// TString DataFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_2Dhistos_sept26.root";
// TString ProtonFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_deep_2Dhistos_sept26_z.root";
// TString NeutronFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_deen_2Dhistos_sept26_z.root";
// TString InelasticFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_inel_2Dhistos.root";// NEED TO RUN INEL
// double xmin_sbs8_100p = -2.2; 
// double xmax_sbs8_100p = 1.0;
// double xminBG_sbs8_100p = -2.2;
// double xmaxBG_sbs8_100p = 1.0;
TString config_title_sbs8_100p = "SBS8 100%";


//// SBS9 70p
// TString DataFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_2Dhistos.root";
// TString ProtonFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_deep_2Dhistos_sept12_z.root";
// TString NeutronFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_deen_2Dhistos_sept12_z.root";
// TString InelasticFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_inel_2Dhistos.root";
// double xmin_sbs9_70p = -1.9; 
// double xmax_sbs9_70p = 1.0;
// double xminBG_sbs9_70p = -1.9;
// double xmaxBG_sbs9_70p = 1.0;
TString config_title_sbs9_70p = "SBS9 70%";


// Function Declarations


// MAIN
void fit_dx_with_background_R(TString KineString="sbs4_30p", TString bg_string="antidy", int PolyOrder_input=5){//main
  //gStyle->SetOptFit(11111);
  gStyle->SetOptFit(0);
  gStyle->SetCanvasPreferGL(1);
  gStyle -> SetOptStat(0);
  gStyle ->SetEndErrorSize(0);
  gStyle->SetOptTitle(1); // Ensure that titles are displayed

  FileNames fileNamesHandler;

  if (KineString == "sbs4_30p"){
    DataFileString = fileNamesHandler.DataFileString_sbs4_30p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs4_30p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs4_30p;
    InelasticFileString =  fileNamesHandler.InelasticFileString_sbs4_30p;
    xmin = fileNamesHandler.xmin_sbs4_30p;
    xmax =fileNamesHandler.xmax_sbs4_30p;
    config_title = config_title_sbs4_30p;
  }else if (KineString == "sbs4_50p"){
    DataFileString = fileNamesHandler.DataFileString_sbs4_50p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs4_50p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs4_50p;
    InelasticFileString =  fileNamesHandler.InelasticFileString_sbs4_50p;
    xmin = fileNamesHandler.xmin_sbs4_50p;
    xmax =fileNamesHandler.xmax_sbs4_50p;
    config_title = config_title_sbs4_50p;
  }else if (KineString == "sbs8_70p"){
    DataFileString = fileNamesHandler.DataFileString_sbs8_70p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs8_70p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs8_70p;
    InelasticFileString =  fileNamesHandler.InelasticFileString_sbs8_70p;
    xmin = fileNamesHandler.xmin_sbs8_70p;
    xmax =fileNamesHandler.xmax_sbs8_70p;
    config_title = config_title_sbs8_70p;
  }else if (KineString == "sbs8_50p"){
    DataFileString = fileNamesHandler.DataFileString_sbs8_50p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs8_50p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs8_50p;
    InelasticFileString =  fileNamesHandler.InelasticFileString_sbs8_50p;
    xmin = fileNamesHandler.xmin_sbs8_50p;
    xmax =fileNamesHandler.xmax_sbs8_50p;
    config_title = config_title_sbs8_50p;
  }else if (KineString == "sbs8_100p"){
    DataFileString = fileNamesHandler.DataFileString_sbs8_100p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs8_100p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs8_100p;
    InelasticFileString =  fileNamesHandler.InelasticFileString_sbs8_100p;
    xmin = fileNamesHandler.xmin_sbs8_100p;
    xmax =fileNamesHandler.xmax_sbs8_100p;;
    config_title = config_title_sbs8_100p;
  }else if (KineString == "sbs9_70p"){
    DataFileString = fileNamesHandler.DataFileString_sbs9_70p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs9_70p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs9_70p;
    InelasticFileString =  fileNamesHandler.InelasticFileString_sbs9_70p;
    xmin = fileNamesHandler.xmin_sbs9_70p;
    xmax =fileNamesHandler.xmax_sbs9_70p;
    config_title = config_title_sbs9_70p;
  }else {
    std::cout<<"Error with kinematic setting"<<std::endl;
    return;
  }
  
  if (bg_string == "antidy"){
    BackgroundFileString = DataFileString;
    BackgroundHistogramName = "hcal_dx_1d_antidy";
  }else if (bg_string == "anticoin"){
    BackgroundFileString = DataFileString;
    BackgroundHistogramName = "hcal_dx_1d_anticoin";
  }else if (bg_string == "mc"){
    BackgroundFileString = InelasticFileString;
    BackgroundHistogramName = "hcal_dx_1d_allcuts";
  }else {
    std::cout<<"Error with background string"<<std::endl;
    return;
  }

  xminBG = xmin - 0.2;
  xmaxBG = xmax + 0.2;
  

  if((PolyOrder_input < 0) && PolyOrder_input > 8){
    cout<<"error with poly order input. Setting to 5."<<endl;
    polyorder = 5;
  }else polyorder=PolyOrder_input; 

  cout<<endl;
  cout<<"What we're running: "<<KineString<<endl;
  cout<<DataFileString<<endl;
  cout<<ProtonFileString<<endl;
  cout<<NeutronFileString<<endl;
  cout<<InelasticFileString<<endl;
  cout<<"Background Histogram: "<<BackgroundHistogramName<<endl;
  cout<<"dx range: ("<<xmin<<", "<<xmax<<")"<<endl;
  cout<< "fitting background to poly order "<<polyorder<<" over range dx: ("<<xminBG<<", "<<xmaxBG<<" )"<<endl; 
  cout<<endl;
 
  // load files 
  TFile *f1 = TFile::Open(DataFileString); // data
  TFile *f2 = TFile::Open(ProtonFileString); //proton
  TFile *f3 = TFile::Open(NeutronFileString); //neutron
  TFile *f4 = TFile::Open(BackgroundFileString); // data file for anticuts, or inelastic file for mc shape 

  // set location and name for output file 
  TString outputfilelocation="../../output/final_Rsf_fits/"+KineString+"/" +bg_string;
  TString outputfilename = outputfilelocation +"/"+ KineString+"_"+bg_string+".root";
  // Declare outfile
  TFile *fout = new TFile(outputfilename,"RECREATE");
  cout<<"writing to file: "<< outputfilename <<endl;
  
  // Load Histograms
  TH1D *hist1 = (TH1D*)f1->Get("hcal_dx_1d_allcuts");
  TH1D *hist2 = (TH1D*)f2->Get("hcal_dx_1d_allcuts");
  TH1D *hist3 = (TH1D*)f3->Get("hcal_dx_1d_allcuts");
  TH1D *hist4 = (TH1D*)f4->Get(BackgroundHistogramName);

  TH1D *hist_data_dont_draw = (TH1D*)hist1->Clone("hist_data_dont_draw"); // fitting with this one
  TH1D *hist_data = (TH1D*)hist1->Clone("hist_data"); /// drawing with this one
  TH1D *hist_p =  (TH1D*)hist2->Clone("hist_p");
  TH1D *hist_n =  (TH1D*)hist3->Clone("hist_n");
  TH1D *histBG = (TH1D*)hist4->Clone("histBG");
  

  if (!hist_data || !hist_p || !hist_n || !histBG) {
    std::cerr << "Error: One of the histograms is null!" << std::endl;
    return;
  }


  
  Utility utilityHandler; // class that gives us access to various functions to use between programs.


  std::vector<double> polyresultBG;
  std::vector<double> polyresultBG_err;

  

  TF1 *FitBG = new TF1("FitBG",Form("pol%i",polyorder),xminBG,xmaxBG);
  // testing
  //polyresultBG.resize(polyorder+1,0.0);
  // FitBG->SetParameters(polyresultBG.data());


  FitBG->SetNpx(500);
  histBG->GetXaxis() ->SetRangeUser(xminBG, xmaxBG);
  histBG->Fit(FitBG, "Q R");
  for (int i = 0; i < polyorder +1; i++)
    {
      polyresultBG.push_back(FitBG->GetParameter(i));
      polyresultBG_err.push_back(FitBG->GetParError(i));			     
    }

 
  FitHistogram fitHandler(hist_p,hist_n,xmin,xmax);
  
  fitHandler.setPolynomialCoefficients(polyresultBG);
 

  fitHandler.fitDataWithBG(hist_data_dont_draw); // bug in here that corrupts it and crashes if you try and draw

  TF1 *fit_result = (TF1*)fitHandler.fitFunc->Clone("fit_result");// bug in here that corrupts it and crashes if you try and draw
  
  /// Make a canvas called c1 to overwrite the default canvas and avoid crashes
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  hist_p->Draw("E");
   

  double scale_p = fitHandler.scale_p;
  double scale_p_err = fitHandler.scale_p_err;
  double scale_n=  fitHandler.scale_n;
  double scale_n_err = fitHandler.scale_n_err;
  double Rsf= fitHandler.R;
  double  Rsf_err= fitHandler.R_err;
  double shift_p = fitHandler.shift_p;
  double shift_p_err = fitHandler.shift_p_err;
  double shift_n = fitHandler.shift_n;
  double shift_n_err = fitHandler.shift_n_err;
  double ChiSq = fitHandler.ChiSq;
  double ndf = fitHandler.NDF;
  double BGscale = fitHandler.BGscale;
  

  cout<<endl<<"--------------------------------------------------------------------------"<<endl;
  cout<<"R = "<<fitHandler.R<<" +/- "<<fitHandler.R_err <<endl;
  cout<<"Shift p = "<< shift_p <<" +/- "<<shift_p_err<<endl;
  cout<<"Shift n = "<<shift_n <<" +/- "<<shift_n_err<<endl;
  cout<<"chi2/ndf = "<< ChiSq<< "/"<<ndf<<" = "<<ChiSq/ndf<<endl;
  cout<<"--------------------------------------------------------------------------"<<endl; 
  
  

  // utility function that scales and shifts the provided histogram and returns a new histogram 
  TH1D *hist_p_result = utilityHandler.ScaleAndShiftHistogram(hist_p, fitHandler.scale_p,fitHandler.shift_p);
  TH1D *hist_n_result = utilityHandler.ScaleAndShiftHistogram(hist_n, fitHandler.scale_n,fitHandler.shift_n);

  
  TF1 *poly_result = new TF1("poly_result", Form("pol%d",polyorder), xmin, xmax);
  for (int i = 0 ; i < polyorder+1; i++)
    {
      double scaled_param = BGscale*polyresultBG[i];
      poly_result ->SetParameter(i,scaled_param);
    }


  // Get the mean and standard deviation from the MC histograms
  double proton_mean_mc = hist2->GetMean();
  double proton_StdDev_mc = hist2->GetStdDev();
  double neutron_mean_mc = hist3->GetMean();
  double neutron_StdDev_mc = hist3->GetStdDev();

  cout<<endl;
  cout<< "Proton: mean = "<<proton_mean_mc<<", StdDev = "<<proton_StdDev_mc<<endl;
  cout<< "Neutron: mean = "<<neutron_mean_mc<<", StdDev = "<<neutron_StdDev_mc<<endl;
  cout<<endl;

  TLine *zero_line = new TLine(xmin, 0, xmax, 0);
  zero_line ->SetLineColor(kRed);
  zero_line->SetLineWidth(2);

  // this will check to see if the background fit went below zero in the range (mosty need for sbs4). 
  bool check_fcn = utilityHandler.doesFunctionGoBelowZero(poly_result,xmin,xmax);

  TH1D *residual_hist = utilityHandler.GetResidualHistogram(hist_data, fit_result);
  residual_hist->GetXaxis() ->SetRangeUser(xmin, xmax);


  // testing out converting thigns to TGraphErrors 
  double markersize =1;

  // Convert histogram to TGraphErrors
  TGraphErrors *data_graph = utilityHandler.histogramToGraphErrors(hist_data); 
  utilityHandler.customizeGraph(data_graph, 20, kBlack, markersize); // assigns shape, color, and marker size to the points
  data_graph->GetXaxis() ->SetRangeUser(xmin, xmax);

  // Convert TF1 fit to TGraphErrors 
  TGraphErrors *fit_graph = utilityHandler.createGraphFromFit(hist_data,fit_result );
  utilityHandler.customizeGraph(fit_graph, 33, kRed, markersize); // assigns shape, color and marker size to the points
  fit_graph->GetXaxis() ->SetRangeUser(xmin, xmax);

   // Add the result histograms and background together for a sanity check
  TH1D* fit_result_hist=utilityHandler.sumHistogramsWithPolynomial(hist_p_result,hist_n_result, poly_result);


  TCanvas *fromMC_unscaled = new TCanvas("fromMC_unscaled","fromMC_unscaled",800,600);
  hist_p->SetLineColorAlpha(kViolet,0.9);
  hist_p->SetLineWidth(1);
  hist_p  ->SetFillColorAlpha(kViolet,0.1);
  hist_p ->Draw("hist E");
  hist_n ->SetLineColorAlpha(kOrange,0.9);
  hist_n->SetLineWidth(1);
  hist_n ->SetFillColorAlpha(kOrange,0.1);
  hist_n ->Draw("hist E same");
  TLegend *mc_unscaled_leg = new TLegend(0.6, 0.7, 0.9, 0.9); // Adjust the legend coordinates as needed
  /// Add entries to the legend for the fit parameters
  mc_unscaled_leg->AddEntry(hist_p,"Proton simc (unscaled unshifted)","f l");
  mc_unscaled_leg->AddEntry(hist_n_result,"Neutron simc (unscaled unshifted)","f l");
  mc_unscaled_leg->Draw();
  fromMC_unscaled->Update();
  fromMC_unscaled ->SaveAs(Form("%s/mc_unscaled_unshifted.pdf",outputfilelocation.Data() ));
  

  TCanvas *canvasBG = new TCanvas("canvasBG","canvasBG",800,600);
  histBG->Draw();
  // poly_result->Draw("same");
  int nEntries_bg = histBG->GetEntries();

  TLegend *bglegend = new TLegend(0.65, 0.7, 0.9, 0.9); // Adjust the legend coordinates as needed
  bglegend->SetTextSize(0.03);  // Adjust this value as needed
  /// Add entries to the legend for the fit parameters
  bglegend->AddEntry(histBG,"Background:" , "l");
  bglegend->AddEntry("", Form("%s", bg_string.Data() ), "");
  bglegend->AddEntry(FitBG,Form("poly%i fit",polyorder) ,"l");
  bglegend->AddEntry("", Form("Entries: %i", nEntries_bg ), "");
  bglegend->Draw();
  
  canvasBG ->SaveAs(Form("%s/background.pdf",outputfilelocation.Data() ));
  
  

  TCanvas *graphcanvas = new TCanvas("graphcanvas","graphcanvas",800,600);   
  //utilityHandler.adjustCanvas(graphcanvas);
  graphcanvas->Divide(1,2);
  // Upper pad
  graphcanvas->cd(1);
  TPad *upperPad1 = (TPad*)gPad;
  upperPad1->SetGrid();
  upperPad1->SetPad(0.01, 0.35, 0.99, 0.99); // Adjust the lower boundary of upper pad
  // utilityHandler.adjustPad(upperPad1);
  //upperPad1->SetFillColor(20);
  //upperPad1->SetFrameFillColor(20);
  // Lower pad
  graphcanvas->cd(2);
  TPad *lowerPad1 = (TPad*)gPad;
  lowerPad1->SetGrid();
  lowerPad1->SetPad(0.01, 0.01, 0.99, 0.34); // Adjust the upper boundary of lower pad
  //utilityHandler.adjustPad(lowerPad1);
  //  lowerPad->SetFillColor(18);
  //lowe1rPad->SetFrameFillColor(18);

  utilityHandler.AdjustHistLabelOffset(hist_n_result);
  utilityHandler.AdjustHistLabelOffset(hist_p_result);
  utilityHandler.AdjustGraphLabelOffset(data_graph);
  utilityHandler.AdjustGraphLabelOffset(fit_graph);
  utilityHandler.AdjustFitLabelOffset(poly_result);
  
  upperPad1->cd();
  data_graph->Draw("AP");
  //data_graph ->SetTitle("hcal dx");
  poly_result ->SetLineColorAlpha(kCyan,0.6);
  poly_result ->SetFillColorAlpha(kCyan,0.35);
  poly_result ->SetFillStyle(1001);
  poly_result ->Draw("same");
  hist_p_result->SetLineColorAlpha(kGreen,0.9);
  hist_p_result->SetLineWidth(1);
  hist_p_result  ->SetFillColorAlpha(kGreen,0.1);
  hist_p_result ->Draw("hist E same");
  hist_n_result ->SetLineColorAlpha(kMagenta,0.9);
  hist_n_result ->SetLineWidth(1);
  hist_n_result ->SetFillColorAlpha(kMagenta,0.1);
  hist_n_result ->Draw("hist E same");
  //data_graph->GetXaxis()->SetTitle("hcal_dx");
  //data_graph->Draw("P SAME");
  fit_graph->Draw("P SAME");

  int nEntries_data = hist_data->GetEntries();
  
  graphcanvas->cd(1); // Ensure you're in the right pad.
  
  TLatex latex;
  double upperleftX = 0.88;
  double upperleftY = 0.48;
  double spacing = 0.08;
  double textsize = 0.07;
  double allign  = 32;
  utilityHandler.DrawLatexRsfLabels(latex, Rsf, Rsf_err, ChiSq, ndf, nEntries_data, upperleftX, upperleftY,spacing,textsize, allign);
  
 
  hist_p_result->SetTitle(""); // Disable default title
  //hist_p_result->GetXaxis()->SetTitle("");/// Disable default title
  hist_n_result->SetTitle(""); // Disable default title
  //hist_n_result->GetXaxis()->SetTitle("");/// Disable default title
  data_graph->SetTitle(""); // Disable default title
  fit_graph->SetTitle(""); // Disable default title

  TLatex title_config;
  title_config.SetNDC();
  title_config.SetTextSize(0.06); // Set desired title size
  title_config.SetTextAlign(12); // left vertical, center horizontal
  title_config.DrawLatex(0.1, 0.94, Form("%s",config_title.Data() ) ); // Adjust x, y for positioning
  
  TLatex title_dx;
  title_dx.SetNDC();
  title_dx.SetTextSize(0.06); // Set desired title size
  title_dx.SetTextAlign(22); // 22 is centered 
  title_dx.DrawLatex(0.5, 0.94, "HCal #Delta x (m)"); // Adjust x, y for positioning
  
  graphcanvas->Update();
  upperPad1->Update();
  lowerPad1->Update();
   

  TLegend *graphlegend = new TLegend(0.6, 0.6, 0.9, 0.9); // Adjust the legend coordinates as needed
  graphlegend->SetTextSize(0.05);  // Adjust this value as needed
  /// Add entries to the legend for the fit parameters
  graphlegend->AddEntry(data_graph, "Data", "P");
  graphlegend->AddEntry(fit_graph,"Overall Fit","P");
  graphlegend->AddEntry(hist_p_result,"Proton simc (scaled)","f l");
  //graphlegend->AddEntry("", Form("   shifted by %.6f ", shift_p), "");
  graphlegend->AddEntry(hist_n_result,"Neutron simc (scaled)","f l");
  //graphlegend->AddEntry("", Form("   shifted by %.6f ", shift_n), "");
  graphlegend->AddEntry(poly_result,Form("Background: %s",bg_string.Data()),"f l");
  //graphlegend->AddEntry("", Form("Rsf = %.6f +/- %.6f ", Rsf,Rsf_err), "");
  //graphlegend->AddEntry("", Form("#chi^{2}/ndf = %.2f / %.0f  ", ChiSq ,ndf), "");
  graphlegend->Draw();

  // draw residual histogram on lower pad
  lowerPad1->cd();
  utilityHandler.AdjustHistPadLabelOffset(residual_hist, lowerPad1);
  residual_hist->SetTitle(""); // Disable default title
  residual_hist ->Draw("E sames");
  zero_line ->Draw("same");
  
  TLatex title_resid;
  title_resid.SetNDC();
  title_resid.SetTextSize(0.08); // Set desired title size
  // Set text alignment: 22 means center alignment both horizontally and vertically
  title_resid.SetTextAlign(22); 
  title_resid.DrawLatex(0.5, 0.94, "Residual Plot (counts)"); // Adjust x, y for positioning

  upperPad1->Update();
  lowerPad1->Update();
  graphcanvas->Update();
  graphcanvas->Draw();
  graphcanvas ->SaveAs(Form("%s/Rsf_graph.pdf",outputfilelocation.Data() ));



  // canvas
  TCanvas *histcanvas = new TCanvas("histcanvas","histcanvas",800,600);
  //divide canvas 
  histcanvas->Divide(1,2);
  // Upper pad
  histcanvas->cd(1);
  TPad *upperPad = (TPad*)gPad;
  upperPad->SetGrid();
  //Double_t xlow, Double_t ylow, Double_t xup, Double_t yup
  upperPad->SetPad(0.01, 0.35, 0.99, 0.99); 
  //utilityHandler.adjustPad(upperPad);
  //upperPad->SetFillColor(20);
  //upperPad->SetFrameFillColor(20);
    
  // Lower pad
  histcanvas->cd(2);
  TPad *lowerPad = (TPad*)gPad;
  lowerPad->SetGrid();
  lowerPad->SetPad(0.01, 0.01, 0.99, 0.34); // Adjust the upper boundary of lower pad
  // utilityHandler.adjustPad(lowerPad);
  //  lowerPad->SetFillColor(18);
  //lowerPad->SetFrameFillColor(18);

  utilityHandler.AdjustHistLabelOffset(hist_n_result);
  utilityHandler.AdjustHistLabelOffset(hist_p_result);
  utilityHandler.AdjustHistLabelOffset(hist_data);

  // Draw histogram and fit function on canvas
  upperPad->cd();
  hist_data->GetXaxis() ->SetRangeUser(xmin, xmax);
  hist_data->Draw("E");
  hist_data->SetLineWidth(2);
  //fit_result->Draw("same");
  poly_result ->SetLineColorAlpha(kCyan,0.9);
  poly_result ->SetFillColorAlpha(kCyan,0.35);
  poly_result ->SetFillStyle(1001);
  poly_result ->Draw("same");
  hist_p_result->SetLineColorAlpha(kGreen,0.9);
  hist_p_result ->SetLineWidth(1);
  hist_p_result ->SetFillColorAlpha(kGreen,0.1);
  hist_p_result ->Draw("hist E same");
  hist_n_result ->SetLineColorAlpha(kMagenta,0.9);
  hist_n_result ->SetLineWidth(1);
  hist_n_result ->SetFillColorAlpha(kMagenta,0.1);
  hist_n_result ->Draw("hist E same");
  fit_result_hist->SetLineColorAlpha(kRed,1);
  fit_result_hist->SetLineWidth(2);
  fit_result_hist->Draw("hist same");
  histcanvas->Update();
  upperPad->Update();

  
  // Create a legend.
  TLegend *legend = new TLegend(0.55, 0.5, 0.9, 0.9); // Adjust the legend coordinates as needed
  legend->SetTextSize(0.05);  // Adjust this value as needed
  /// Add entries to the legend for the fit parameters
  legend->AddEntry(hist_data, "Data", "l");
  legend->AddEntry(fit_result,"Overall Fit","l");
  legend->AddEntry(hist_p_result,"Proton simc","f l");
  legend->AddEntry("", Form("   shifted by %.6f ", shift_p), "");
  legend->AddEntry(hist_n_result,"Neutron simc","f l");
  legend->AddEntry("", Form("   shifted by %.6f ", shift_n), "");
  legend->AddEntry(poly_result,Form("Background: %s",bg_string.Data()),"f l");
  legend->AddEntry("", Form("Rsf = %.6f +/- %.6f ", Rsf,Rsf_err), "");
  legend->AddEntry("", Form("#chi^{2}/ndf = %.2f / %.0f = %.2f", ChiSq ,ndf,ChiSq/ndf), "");
  legend->Draw();

  // draw residual histogram on lower pad
  lowerPad->cd();
  //utilityHandler.AdjustHistLabelOffset(residual_hist);
  residual_hist ->Draw("E sames");
  zero_line->Draw("same");

  upperPad->Update(); // this is supposed to update the axis labels but it looks like crap
  lowerPad->Update();
  histcanvas->Update();
  histcanvas ->SaveAs(Form("%s/Rsf_hist.pdf",outputfilelocation.Data() ));

  TCanvas *residcanvas = new TCanvas("residcanvas","residcanvas",800,600);
  residcanvas->SetGrid();
  residual_hist ->Draw("E");
  zero_line->Draw("same");
  residcanvas ->SaveAs(Form("%s/resid.pdf",outputfilelocation.Data() ));


  std::string bg_title= histBG->GetTitle();
  utilityHandler.printParsedTitle(bg_title,outputfilelocation,"background");
   
  std::string title = hist_data->GetTitle();
  utilityHandler.printParsedTitle(title,outputfilelocation,"data");
 
  fout->Write();
  f1->Close();
  f2->Close();
  f3->Close();
  f4->Close();

}// end main

