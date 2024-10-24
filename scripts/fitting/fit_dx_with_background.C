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


/// This program will load in a histogram that is going to form the shape of the background.
/// For example, you can load in the anti-coin histogram.
/// It will fit the background shape to a high order polynomial.
/// Then use a scaled version of it as the overall fit. 

// Global Params
Double_t wide_x_range = 0.2;
Double_t nsigma = 0.6;
Double_t xmin = -2.1; // -2.1 for sbs4 30p
Double_t xmax = 1.4;//1.4 for sbs4 30p


const int polyorder = 7;
std::vector<double> polyresultBG;
std::vector<double> polyresultBG_err;
Double_t xminBG = -2.5; // -1.8
Double_t xmaxBG = 1.4; //1;




// Function Declarations
Double_t gausFit(Double_t *x, Double_t *par);
Double_t poly2(Double_t *x, Double_t *par);
Double_t poly3(Double_t *x, Double_t *par);
Double_t poly4(Double_t *x, Double_t *par);
Double_t poly2_2gaus(Double_t *x, Double_t *par);
Double_t poly4_2gaus(Double_t *x, Double_t *par);
Double_t mc_p_n_poly4_fit(Double_t *x, Double_t *par);
Double_t mc_p_n_poly2_fit(Double_t *x, Double_t *par);
Double_t mc_p_n_BG_fit(Double_t *x, Double_t *par);
Double_t mc_p_fit(Double_t *x, Double_t *par);
Double_t mc_n_fit(Double_t *x, Double_t *par);
TH1D* shiftHistogramX(TH1D* originalHist, double shiftValue);
TH1D* GetResidualHistogram(TH1D* hist, TF1* fit);
TGraphErrors*  histogramToGraphErrors(TH1D *hist);
void customizeGraph(TGraphErrors *graph, int markerStyle, int markerColor,double markersize);
TGraphErrors* createGraphFromFit(TH1D* hist, TF1* fitFunc);
void adjustCanvas(TCanvas* canvas,
                  double leftMargin = 0.15, double rightMargin = 0.05, 
                  double bottomMargin = 0.15, double topMargin = 0.10);
void adjustPad(TPad* pad,
                  double leftMargin = 0.15, double rightMargin = 0.05, 
		double bottomMargin = 0.15, double topMargin = 0.10);
void AdjustHistLabelOffset(TH1D* hist, double xoffset = 0.03, double yoffset= 0.02);
void printParsedTitle(const std::string& title);

/// making these histograms globals so they can be accessed by thte fit funcitons 
TH1D *histP;
TH1D *histN;
TH1D *histBG;

// MAIN
void fit_dx_with_background(){//main
  gStyle->SetOptFit(11111);
  gStyle->SetCanvasPreferGL(1);
  gStyle -> SetOptStat(0);
  gStyle ->SetEndErrorSize(0);


  // Load rootfiles for histograms. Should have this happen with a config file or something eventually. 
  //TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_sbs30p_sf0p359.root"); 
  // TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_feb22.root"); 
  //TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_mar14.root"); 
  // TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_march13.root"); 
  //  TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_may18.root"); // data
  //TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p.root"); // data
  // TFile *f2 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_simc_deep_may9.root"); //proton
  // TFile *f2 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_simc_deep.root"); // Load rootfile
  //TFile *f3 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_simc_deen_may9.root"); //neutron
  // TFile *f3 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_simc_deen.root"); 

  
  // Using histograms that come from the Make2Dhistos scripts.  

 TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_tight_2Dhistos.root"); // data
 TFile *f2 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_tight_2Dhistos.root"); // proton
 TFile *f3 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_tight_2Dhistos.root"); // neutron
 TFile *f4 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_inel_2Dhistos.root"); // inelastic sim
 

  // TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_2Dhistos.root"); //data
  // TFile *f2 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deep_2Dhistos.root"); ///proton
  // TFile *f3 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deen_2Dhistos.root"); //neutron
  // TFile *f4 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_inel_2Dhistos.root");/// inelastic sim


// TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_tight_BinStudy.root"); data
//  TFile *f2 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_tight_BinStudy.root"); proton
//  TFile *f3 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_tight_BinStudy.root"); neutron

 
  // Load Histograms
 TH1D *hist1 = (TH1D*)f1->Get("hcal_dx_1d_allcuts");// data
 TH1D *hist2 = (TH1D*)f2->Get("hcal_dx_1d_allcuts");// proton
 TH1D *hist3 = (TH1D*)f3->Get("hcal_dx_1d_allcuts");// neutron

 //TH1D *hist4=(TH1D*)f1->Get("hcal_dx_1d_antidy");// Anticut histogram from data that we will use as a background shape.
 TH1D *hist4=(TH1D*)f4->Get("hcal_dx_1d_allcuts");// dx histogram from the inelastic simulation 



  histP = (TH1D*)hist2->Clone("histP");
  histN = (TH1D*)hist3->Clone("histN");
  histBG=(TH1D*)hist4->Clone("histBG");


  TH1D *hist1_clone = (TH1D*)hist1->Clone("hist1_clone");
  hist1_clone->GetXaxis() ->SetRangeUser(xmin, xmax);

  //// Get max bin to find where the proton peak probaby is 
  double maxBin = hist1 -> GetMaximumBin();
  double maxBinValue = hist1 -> GetBinContent(maxBin);
  double maxBinXValue = hist1->GetXaxis()->GetBinCenter(maxBin);
  cout<<"maxBin: "<<maxBin<<endl;
  cout<<"maxBinValue: "<<maxBinValue<<endl;
  cout<<"maxBinXValue: "<<maxBinXValue<<endl;
  cout<<endl;


  // Get the mean and standard deviation from the MC histograms
  double proton_mean_mc = hist2->GetMean();
  double proton_StdDev_mc = hist2->GetStdDev();
  double neutron_mean_mc = hist3->GetMean();
  double neutron_StdDev_mc = hist3->GetStdDev();

  cout<< "Proton: mean = "<<proton_mean_mc<<", StdDev = "<<proton_StdDev_mc<<endl;
  cout<< "Neutron: mean = "<<neutron_mean_mc<<", StdDev = "<<neutron_StdDev_mc<<endl;
  cout<<endl;

  TLine *zero_line = new TLine(xmin, 0, xmax, 0);
  zero_line ->SetLineColor(kRed);
  zero_line->SetLineWidth(2);

  //// fit the background histogram to a high order polynomail
  /// polyorder is set up as a global param

  TF1 *FitBG = new TF1("FitBG",Form("pol%d",polyorder),xminBG,xmaxBG); // make sure this matches polyorder
  FitBG->SetNpx(500);
  histBG->GetXaxis() ->SetRangeUser(xminBG, xmaxBG);
  histBG->Fit(FitBG, "Q R");
  for (int i = 0; i < polyorder +1; i++)
    {
      polyresultBG.push_back(FitBG->GetParameter(i));
      polyresultBG_err.push_back(FitBG->GetParError(i));			     
    }

  
  //Double_t mc_p_n_BG_fit(Double_t *x, Double_t *par)

  double initialParameters[5] = {1,1,0,0,1};

   // Create custom fit function
    TF1 *fitFunc = new TF1("fitFunc", mc_p_n_BG_fit, xmin, xmax, 5); //

    // Set initial parameters for the fit function
    fitFunc->SetParameters(initialParameters); // Define initial parameters
    fitFunc ->SetNpx(500);
    fitFunc->SetParLimits(2, -0.10, 0.10); // setting the limit on parameter 2 (proton shift)
     fitFunc->SetParLimits(3, -0.10, 0.10); // setting the limit on parameter 3 (neutron shift)
     fitFunc->SetParLimits(4, 0, 10000); // setting the scale of the BG to be positive  
     // Fit combined histogram with custom fit function
     hist1->GetXaxis() ->SetRangeUser(xmin, xmax);
     hist1->Fit(fitFunc, "Q R");

    double normP_result, normN_result;
    double normP_result_error, normN_result_error;
    double Pshift_result, Nshift_result;
     double Pshift_result_error, Nshift_result_error;
    double BG_scale_result, BG_scale_result_error;
    double ChiSq, ndf;

    normP_result = fitFunc ->GetParameter(0);
    normP_result_error = fitFunc ->GetParError(0);
    normN_result = fitFunc ->GetParameter(1);
    normN_result_error = fitFunc ->GetParError(1);
    Pshift_result = fitFunc ->GetParameter(2);
    Pshift_result_error = fitFunc ->GetParError(2);
    Nshift_result = fitFunc ->GetParameter(3);
    Nshift_result_error = fitFunc ->GetParError(3);
    BG_scale_result = fitFunc->GetParameter(4);
    BG_scale_result_error = fitFunc->GetParError(4);
    
    
    ChiSq = fitFunc->GetChisquare();
    ndf = fitFunc->GetNDF();

    cout<<"Pshift " <<Pshift_result <<" +/- "<<Pshift_result_error<<endl;
    cout<<"Nshift " <<Nshift_result <<" +/- "<<Nshift_result_error<<endl;

    cout<<"BG scale " <<BG_scale_result<<" +/- "<< BG_scale_result_error<<endl;

      TH1D *histP_clone = shiftHistogramX(histP,  Pshift_result);
      TH1D *histN_clone =shiftHistogramX(histN, Nshift_result); 

      histP_clone ->Scale(normP_result); 
      histN_clone->Scale(normN_result);

      double Ratio = normN_result / normP_result;
      
      cout<< "ratio: scale_n / scale_p = "<< normN_result <<" / " <<normP_result <<" = "<<normN_result / normP_result <<endl;

      double Ratio_error = Ratio * sqrt( pow( (normN_result_error / normN_result), 2) + pow( (normP_result_error / normP_result),2) ); //just adding the uncert from the fit parameters in quadrature for now. 
      
      cout<<"ratio error = " << Ratio_error<<endl;

      TF1 *P_result =  new TF1("P_result", mc_p_fit,xmin,xmax,1);
      P_result ->SetParameter(0,normP_result);
      P_result ->SetParameter(1,Pshift_result);

      TF1 *N_result =  new TF1("N_result", mc_n_fit, xmin,xmax,1);
      N_result ->SetParameter(0,normN_result);
      N_result ->SetParameter(1,Nshift_result);

      TF1 *poly_result = new TF1("poly_result", Form("pol%d",polyorder), xmin, xmax);
      for (int i = 0 ; i < polyorder+1; i++)
	{
	  double scaled_param = BG_scale_result*polyresultBG[i];
	  poly_result ->SetParameter(i,scaled_param);
	}
      

      TH1D *residual_hist = GetResidualHistogram(hist1, fitFunc);
      residual_hist->GetXaxis() ->SetRangeUser(xmin, xmax);


      // testing out converting thigns to TGraphErrors 
      double markersize =1;

      // Convert histogram to TGraphErrors
      TGraphErrors *data_graph = histogramToGraphErrors(hist1_clone); 
      customizeGraph(data_graph, 20, kBlack, markersize); // assigns shape, color, and marker size to the points
      data_graph->GetXaxis() ->SetRangeUser(xmin, xmax);

      // Convert TF1 fit to TGraphErrors 
      TGraphErrors *fit_graph = createGraphFromFit(hist1,fitFunc );
      customizeGraph(fit_graph, 33, kRed, markersize); // assigns shape, color and marker size to the points
      fit_graph->GetXaxis() ->SetRangeUser(xmin, xmax);

      TCanvas *graphcanvas = new TCanvas("graphcanvas","graphcanvas",800,600);   
      adjustCanvas(graphcanvas);
      graphcanvas->Divide(1,2);
      // Upper pad
      graphcanvas->cd(1);
      TPad *upperPad1 = (TPad*)gPad;
      upperPad1->SetPad(0.01, 0.3, 0.99, 0.99); //Double_t xlow, Double_t ylow, Double_t xup, Double_t yup
      adjustPad(upperPad1);
      //upperPad1->SetFillColor(20);
      //upperPad1->SetFrameFillColor(20);
// Lower pad
      graphcanvas->cd(2);
      TPad *lowerPad1 = (TPad*)gPad;
      lowerPad1->SetPad(0.01, 0.01, 0.99, 0.29);
      adjustPad(lowerPad1);
      //  lowerPad->SetFillColor(18);
      //lowe1rPad->SetFrameFillColor(18);

      upperPad1->cd();
      data_graph->Draw("AP");
      poly_result ->SetLineColorAlpha(kCyan,0.6);
      poly_result ->SetFillColorAlpha(kCyan,0.35);
      poly_result ->SetFillStyle(1001);
      poly_result ->Draw("same");
      histP_clone->SetLineColorAlpha(kGreen,0.9);
      histP_clone ->SetLineWidth(1);
      histP_clone ->SetFillColorAlpha(kGreen,0.1);
      histP_clone ->Draw("hist same");
      histN_clone ->SetLineColorAlpha(kMagenta,0.9);
      histN_clone ->SetLineWidth(1);
      histN_clone ->SetFillColorAlpha(kMagenta,0.1);
      histN_clone ->Draw("hist same");
      data_graph->Draw("P SAME");
      fit_graph->Draw("P SAME");

   

      TLegend *graphlegend = new TLegend(0.6, 0.7, 0.9, 0.9); // Adjust the legend coordinates as needed
      /// Add entries to the legend for the fit parameters
      graphlegend->AddEntry(data_graph, "Data", "P");
      graphlegend->AddEntry(fit_graph,"Overall Fit","P");
      graphlegend->AddEntry(histP_clone,"Proton simc","f");
      graphlegend->AddEntry("", Form("   shifted by %.4f ", Pshift_result), "");
      graphlegend->AddEntry(histN_clone,"Neutron simc","f");
      graphlegend->AddEntry("", Form("   shifted by %.4f ", Nshift_result), "");
      graphlegend->AddEntry(poly_result,"background","f");
      graphlegend->AddEntry("", Form("R= %.4f +/- %.4f ", Ratio,Ratio_error), "");
      graphlegend->AddEntry("", Form("#chi^{2}/ndf = %.2f / %.0f  ", ChiSq ,ndf), "");
      graphlegend->Draw();

      // draw residual histogram on lower pad
      lowerPad1->cd();
      AdjustHistLabelOffset(residual_hist);
      residual_hist ->Draw("E sames");
      zero_line->Draw("same");

      graphcanvas->Update();
      graphcanvas->Draw();




      // canvas
      TCanvas *c2 = new TCanvas("c2","c2",800,600);
      //divide canvas 
      c2->Divide(1,2);
      // Upper pad
      c2->cd(1);
      TPad *upperPad = (TPad*)gPad;
      upperPad->SetPad(0.01, 0.3, 0.99, 0.99); //Double_t xlow, Double_t ylow, Double_t xup, Double_t yup
      adjustPad(upperPad);
      //upperPad->SetFillColor(20);
      //upperPad->SetFrameFillColor(20);
    
      // Lower pad
      c2->cd(2);
      TPad *lowerPad = (TPad*)gPad;
      lowerPad->SetPad(0.01, 0.01, 0.99, 0.29);
      adjustPad(lowerPad);
      //  lowerPad->SetFillColor(18);
      //lowerPad->SetFrameFillColor(18);

      // Draw histogram and fit function on canvas
      upperPad->cd();
      hist1->Draw("E");
      hist1->SetLineWidth(2);
      //fitFunc->Draw("same");
      poly_result ->SetLineColorAlpha(kCyan,0.9);
      poly_result ->SetFillColorAlpha(kCyan,0.35);
      poly_result ->SetFillStyle(1001);
      poly_result ->Draw("same");
      histP_clone->SetLineColorAlpha(kGreen,0.9);
      histP_clone ->SetLineWidth(1);
      histP_clone ->SetFillColorAlpha(kGreen,0.1);
      histP_clone ->Draw("hist same");
      histN_clone ->SetLineColorAlpha(kMagenta,0.9);
      histN_clone ->SetLineWidth(1);
      histN_clone ->SetFillColorAlpha(kMagenta,0.1);
      histN_clone ->Draw("hist same");

      //P_result ->SetLineColor(kBlack);
      // P_result ->Draw("same");
    
   
  // Create a legend.
  TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9); // Adjust the legend coordinates as needed
  /// Add entries to the legend for the fit parameters
  legend->AddEntry(hist1, "Data", "l");
  legend->AddEntry(fitFunc,"Overall Fit","l");
  legend->AddEntry(histP_clone,"Proton simc","f");
  legend->AddEntry("", Form("   shifted by %.4f ", Pshift_result), "");
  legend->AddEntry(histN_clone,"Neutron simc","f");
 legend->AddEntry("", Form("   shifted by %.4f ", Nshift_result), "");
  legend->AddEntry(poly_result,"Background","f");
  legend->AddEntry("", Form("R= %.4f +/- %.4f ", Ratio,Ratio_error), "");
   legend->AddEntry("", Form("#chi^{2}/ndf = %.2f / %.0f  ", ChiSq ,ndf), "");
  legend->Draw();

  // draw residual histogram on lower pad
  lowerPad->cd();
  AdjustHistLabelOffset(residual_hist);
  residual_hist ->Draw("E sames");
  zero_line->Draw("same");

  upperPad->Update(); // this is supposed to update the axis labels but it looks like crap
  lowerPad->Update();
  c2->Update();

  TCanvas *c3 = new TCanvas("c3","c3",800,600);
  residual_hist ->Draw("E hist");

   TCanvas *c4 = new TCanvas("c4","c4",800,600);
   // histP->SetLineColor(kRed);
   //histP->SetFillColor(kRed);
   // histP->SetFillStyle(3003);
   histP_clone->SetLineColorAlpha(kGreen,1);
   histP_clone->Draw("E");
   histN_clone->SetLineColorAlpha(kMagenta,1);
   histN_clone->Draw("sames E");


   TCanvas *canvasBG = new TCanvas("canvasBG","canvasBG",800,600);
   histBG->Draw();
   poly_result->Draw("same");

   std::string title = hist1->GetTitle();
   printParsedTitle(title);

    std::string title2 = hist4->GetTitle();
    cout<<endl;
    cout<<"background"<<endl;
    cout<<title2<<endl;



 //   /// set up a gaussian fit for the proton peak 
 //  TF1* gaussian_init =  new TF1("gaussian_init",gausFit , maxBinXValue -wide_x_range , maxBinXValue + wide_x_range,3);
 //  gaussian_init ->SetParameters(maxBinValue, maxBinXValue, 0.2); // amp, mean, sigma 
 //  /// fit histogram
 //  gaussian_init ->SetParNames("amp", "mean","sigma");
 //  hist1->Fit(gaussian_init, "Q", "", maxBinXValue - wide_x_range, maxBinXValue + wide_x_range);
 //  //get fit results
 //  double amp_init =  gaussian_init -> GetParameter(0);
 //  double mean_init = gaussian_init -> GetParameter(1);
 //  double sigma_init = gaussian_init ->GetParameter(2);
 //  cout<<"amp init: "<<amp_init<<endl;
 //  cout<<"mean init: "<<mean_init<<endl;
 //  cout<<"sigma init: "<<sigma_init<<endl;
 //  cout<<endl;


 //  /// set up a gaussian fit for the neutron peak. Mostly the same as proton, but around zero. 
 //  TF1* gaussian_init2 =  new TF1("gaussian_init2",gausFit , 0-wide_x_range , 0+wide_x_range,3);
 //  gaussian_init2 ->SetParameters(maxBinValue/2, 0, 0.2); // amp, mean, sigma 
 //  /// fit histogram
 //  gaussian_init2 ->SetParNames("amp", "mean","sigma");
 //  hist1->Fit(gaussian_init2, "Q", "", 0 - 0.25, 0 + 0.25);
 //  //get fit results
 //  double amp_init2 =  gaussian_init2 -> GetParameter(0);
 //  double mean_init2 = gaussian_init2 -> GetParameter(1);
 //  double sigma_init2 = gaussian_init2 ->GetParameter(2);
 //  cout<<"amp init2: "<<amp_init2<<endl;
 //  cout<<"mean init2: "<<mean_init2<<endl;
 //  cout<<"sigma init2: "<<sigma_init2<<endl;
 //  cout<<endl;



 //  /// fit again with a more contrained gaussian using the paramters from the looser gaussian for proton peak. 
 //  TF1* gaussian =  new TF1("gaussian",gausFit , mean_init - nsigma*sigma_init , mean_init + nsigma*sigma_init,3);
 //  gaussian ->SetParNames("amp", "mean","sigma");
 //  gaussian->SetParameters(amp_init, mean_init, sigma_init); // amp, mean, sigma 
 //  hist1->Fit(gaussian, "Q", "", mean_init - nsigma*sigma_init , mean_init + nsigma*sigma_init);
 //  //get fit results
 //  double amp =  gaussian -> GetParameter(0);
 //  double mean = gaussian -> GetParameter(1);
 //  double sigma = gaussian ->GetParameter(2);
 //  cout<<"amp: "<<amp<<endl;
 //  cout<<"mean: "<<mean<<endl;
 //  cout<<"sigma: "<<sigma<<endl;
 //  cout<<endl;


 //  /// fit again with a more contrained gaussian using the paramters from the looser gaussian for neutron peak. 
 //  TF1* gaussian2 =  new TF1("gaussian2",gausFit , mean_init2 - nsigma*sigma_init2 , mean_init2 + nsigma*sigma_init2,3);
 //  gaussian2 ->SetParNames("amp", "mean","sigma");
 //  gaussian2->SetParameters(amp_init2, mean_init2, sigma_init2); // amp, mean, sigma 
 //  hist1->Fit(gaussian2, "Q", "", mean_init2 - nsigma*sigma_init2 , mean_init2 + nsigma*sigma_init2);
 //  //get fit results
 //  double amp2 =  gaussian2 -> GetParameter(0);
 //  double mean2 = gaussian2 -> GetParameter(1);
 //  double sigma2= gaussian2 ->GetParameter(2);
 //  cout<<"amp2: "<<amp2<<endl;
 //  cout<<"mean2: "<<mean2<<endl;
 //  cout<<"sigma2: "<<sigma2<<endl;
 //  cout<<endl;

 //  // // Fit the whole distribution to two gaussians plus a polynomail background, using the gauss params we found as initial values. 
 //  // TF1* poly_2gaus = new TF1("poly_2gaus",poly2_2gaus,xmin,xmax,9);
 //  // //poly_2gaus ->SetParNames();
 //  // poly_2gaus ->SetParameters(100,-10,-10,amp,mean,sigma,amp2,mean2,sigma2);
 //  // poly_2gaus ->SetParLimits(7,mean2 - wide_x_range, mean2 + wide_x_range);
 //  // hist1->Fit(poly_2gaus, "Q", "", xmin, xmax);
 //  // double par0 = poly_2gaus -> GetParameter(0);
 //  // double par0_err = poly_2gaus ->GetParError(0);
 //  // double par1 = poly_2gaus ->GetParameter(1);
 //  // double par1_err = poly_2gaus ->GetParError(1);
 //  // double par2 = poly_2gaus ->GetParameter(2);
 //  // double par2_err = poly_2gaus ->GetParError(2);
 //  // double amp_p = poly_2gaus ->GetParameter(3);
 //  // double amp_p_err = poly_2gaus ->GetParError(3);
 //  // double mean_p = poly_2gaus ->GetParameter(4);
 //  // double mean_p_err = poly_2gaus ->GetParError(4);
 //  // double sigma_p = poly_2gaus ->GetParameter(5);
 //  // double sigma_p_err = poly_2gaus ->GetParError(5);
 //  // double amp_n = poly_2gaus ->GetParameter(6);
 //  // double amp_n_err = poly_2gaus ->GetParError(6);
 //  // double mean_n = poly_2gaus ->GetParameter(7);
 //  // double mean_n_err = poly_2gaus ->GetParError(7);
 //  // double sigma_n = poly_2gaus ->GetParameter(8);
 //  // double sigma_n_err = poly_2gaus ->GetParError(8);
 //  // TF1* poly_result = new TF1("poly_result",poly2,xmin,xmax,3);
 //  // poly_result ->SetParameters(par0,par1,par2);
 //  // TF1* gaus_result_p = new TF1("gaus_result_p",gausFit, xmin,xmax,3);
 //  // gaus_result_p ->SetParameters(amp_p, mean_p, sigma_p);
 //  // TF1* gaus_result_n = new TF1("gaus_result_n",gausFit,xmin,xmax,3);
 //  // gaus_result_n ->SetParameters(amp_n,mean_n, sigma_n);
 //  // cout << "Overall Fit"<<endl;
 //  // cout<<"par0: "<<par0<<" +- "<<par0_err<<endl;
 //  // cout<<"par1: "<<par1<<" +- "<<par1_err<<endl;
 //  // cout<<"par2: "<<par2<<" +- "<<par2_err<<endl;
 //  // cout<<"amp_p: "<<amp_p<<" +- "<<amp_p_err<<endl;
 //  // cout<<"mean_p: "<<mean_p <<" +- "<<mean_p_err<<endl;
 //  // cout<<"sigma_p: "<<sigma_p <<" +- "<<sigma_p_err<<endl;
 //  // cout<<"amp_n: "<<amp_n<<" +- "<<amp_n_err<<endl;
 //  // cout<<"mean_n: "<<mean_n <<" +- "<<mean_n_err<<endl;
 //  // cout<<"sigma_n: "<<sigma_n <<" +- "<<sigma_n_err<<endl;
 //  // cout<<endl;




 // // Fit the whole distribution to two gaussians plus a 4th orderpolynomail background, using the gauss params we found as initial values. 
 //  TF1* poly_2gaus = new TF1("poly_2gaus",poly4_2gaus,xmin,xmax,11);
 //  //poly_2gaus ->SetParNames();
 //  poly_2gaus ->SetParameters(100,-10,-10,-10,-10,amp,mean,sigma,amp2,mean2,sigma2);
 //  poly_2gaus ->SetParLimits(9,mean2 - wide_x_range, mean2 + wide_x_range);
 //  hist1->Fit(poly_2gaus, "Q", "", xmin, xmax);
 //  double par0 = poly_2gaus -> GetParameter(0);
 //  double par0_err = poly_2gaus ->GetParError(0);
 //  double par1 = poly_2gaus ->GetParameter(1);
 //  double par1_err = poly_2gaus ->GetParError(1);
 //  double par2 = poly_2gaus ->GetParameter(2);
 //  double par2_err = poly_2gaus ->GetParError(2);
 //  double par3 = poly_2gaus ->GetParameter(3);
 //  double par3_err = poly_2gaus ->GetParError(3);
 //  double par4 = poly_2gaus ->GetParameter(4);
 //  double par4_err = poly_2gaus ->GetParError(4);
 //  double amp_p = poly_2gaus ->GetParameter(5);
 //  double amp_p_err = poly_2gaus ->GetParError(5);
 //  double mean_p = poly_2gaus ->GetParameter(6);
 //  double mean_p_err = poly_2gaus ->GetParError(6);
 //  double sigma_p = poly_2gaus ->GetParameter(7);
 //  double sigma_p_err = poly_2gaus ->GetParError(7);
 //  double amp_n = poly_2gaus ->GetParameter(8);
 //  double amp_n_err = poly_2gaus ->GetParError(8);
 //  double mean_n = poly_2gaus ->GetParameter(9);
 //  double mean_n_err = poly_2gaus ->GetParError(9);
 //  double sigma_n = poly_2gaus ->GetParameter(10);
 //  double sigma_n_err = poly_2gaus ->GetParError(10);
 //  TF1* poly_result = new TF1("poly_result",poly4,xmin,xmax,5);
 //  poly_result ->SetParameters(par0,par1,par2,par3,par4);
 //  TF1* gaus_result_p = new TF1("gaus_result_p",gausFit, xmin,xmax,3);
 //  gaus_result_p ->SetParameters(amp_p, mean_p, sigma_p);
 //  TF1* gaus_result_n = new TF1("gaus_result_n",gausFit,xmin,xmax,3);
 //  gaus_result_n ->SetParameters(amp_n,mean_n, sigma_n);
 //  cout << "Overall Fit"<<endl;
 //  cout<<"par0: "<<par0<<" +- "<<par0_err<<endl;
 //  cout<<"par1: "<<par1<<" +- "<<par1_err<<endl;
 //  cout<<"par2: "<<par2<<" +- "<<par2_err<<endl;
 //  cout<<"par3: "<<par3<<" +- "<<par3_err<<endl;
 //  cout<<"par4: "<<par4<<" +- "<<par4_err<<endl;
 //  cout<<"amp_p: "<<amp_p<<" +- "<<amp_p_err<<endl;
 //  cout<<"mean_p: "<<mean_p <<" +- "<<mean_p_err<<endl;
 //  cout<<"sigma_p: "<<sigma_p <<" +- "<<sigma_p_err<<endl;
 //  cout<<"amp_n: "<<amp_n<<" +- "<<amp_n_err<<endl;
 //  cout<<"mean_n: "<<mean_n <<" +- "<<mean_n_err<<endl;
 //  cout<<"sigma_n: "<<sigma_n <<" +- "<<sigma_n_err<<endl;
 //  cout<<endl;


 //  // Set up canvas 
 //  TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
 //  hist1->SetLineColor(kBlue);
 //  hist1->Draw();
 //  //gaussian ->Draw("same");
 //  //gaussian2 ->Draw("same");
 //  poly_2gaus ->SetLineColor(kGreen);
 //  poly_2gaus ->Draw("same");
 //  poly_result ->SetLineColor(kRed);
 //  poly_result ->SetLineStyle(2);
 //  poly_result ->Draw("same");
 //  gaus_result_p ->SetLineColor(kMagenta); 
 //  gaus_result_p ->SetLineStyle(2);
 //  gaus_result_p ->Draw("same");
 //  gaus_result_n ->SetLineColor(kViolet);
 //  gaus_result_n ->SetLineStyle(2);
 //  gaus_result_n ->Draw("same");


 //  // Create a legend.
 //  TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9); // Adjust the legend coordinates as needed
 //  /// Add entries to the legend for the fit parameters
 //  legend->AddEntry(hist1, "Data", "l");
 //  legend->AddEntry(poly_2gaus,"Overall Fit","l");
 //  legend->AddEntry(gaus_result_p,"Proton Peak","l");
 //  legend->AddEntry("", Form("Mean= %.4f +- %.4f", mean_p,mean_p_err), "");
 //  legend->AddEntry("", Form("Sigma= %.4f +-%.4f",sigma_p, sigma_p_err), "");
 //  legend->AddEntry(gaus_result_n,"Neutron Peak","l");
 //  legend->AddEntry("", Form("Mean= %.4f +- %.4f", mean_n, mean_n_err), "");
 //  legend->AddEntry("", Form("Sigma= %.4f +- %.4f",sigma_n,sigma_n_err), "");
 //  legend->AddEntry(poly_result,"Poly Background","l");
 //  legend->Draw();


 //    c1->Update();


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


// Fit that is a combination of the scaled proton mc, scaled neutron mc, and 4th order polynomial. 
//// histP and histN are globals that need to be set to the proper histograms right before you call the fit function. 
Double_t mc_p_n_poly4_fit(Double_t *x, Double_t *par) {
    
    Double_t val = 0.0;

    // Get x value
    Double_t xx = x[0];

    // Retrieve parameters
    Double_t normP = par[0];
    Double_t normN = par[1];
    Double_t Pshift = par[2];
    Double_t Nshift = par[3];

    Double_t polyCoefficients[5]; 
    


    // Fill polynomial coefficients
    for (Int_t i = 0; i <=4; ++i) {
        polyCoefficients[i] = par[i + 4];
    }

    // Calculate value using combination of histograms and polynomial background
    val = normP * histP->Interpolate(xx-Pshift) + normN * histN->Interpolate(xx-Nshift);

    // Add polynomial background
    for (Int_t i = 0; i <= 4; ++i) {
        val += polyCoefficients[i] * TMath::Power(xx, i);
    }

    return val;
}

// Fit that is a combination of the scaled proton mc, scaled neutron mc, and 2nd order polynomial. 
//// histP and histN are globals that need to be set to the proper histograms right before you call the fit function. 
Double_t mc_p_n_poly2_fit(Double_t *x, Double_t *par) {
    
    Double_t val = 0.0;

    // Get x value
    Double_t xx = x[0];

    // Retrieve parameters
    Double_t normP = par[0];
    Double_t normN = par[1];
    Double_t Pshift = par[2];
    Double_t Nshift = par[3];

    Double_t polyCoefficients[3]; 
    


    // Fill polynomial coefficients
    for (Int_t i = 0; i <=2; ++i) {
        polyCoefficients[i] = par[i + 4];
    }

    // Calculate value using combination of histograms and polynomial background
    val = normP * histP->Interpolate(xx-Pshift) + normN * histN->Interpolate(xx-Nshift);

    // Add polynomial background
    for (Int_t i = 0; i <= 2; ++i) {
        val += polyCoefficients[i] * TMath::Power(xx, i);
    }

    return val;
}


// Fit that is a combination of the scaled proton mc, scaled neutron mc, and a scaleled version of the fit to a historam.
//// histP and histN and polyresultBG[] are globals that need to be set to the proper histograms right before you call the fit function. 
Double_t mc_p_n_BG_fit(Double_t *x, Double_t *par) {
    
    Double_t val = 0.0;

    // Get x value
    Double_t xx = x[0];

    // Retrieve parameters
    Double_t normP = par[0];
    Double_t normN = par[1];
    Double_t Pshift = par[2];
    Double_t Nshift = par[3];

    Double_t normBG = par[4];

    
    // Calculate value using combination of histP and histN histograms
    val = normP * histP->Interpolate(xx-Pshift) + normN * histN->Interpolate(xx-Nshift);


    // Add polynomial background
    for (Int_t i = 0; i <= polyorder+1; ++i) {
      val += normBG*polyresultBG[i]*TMath::Power(xx,i);
    }
	
    return val;
}




// // Custom fit function
// Double_t mc_p_n_poly4_fit(Double_t *x, Double_t *par) {
//     // Parameters:
//     // par[0] : normalization of hist1
//     // par[1] : normalization of hist2
   

//     Double_t val = 0.0;

//     // Get x value
//     Double_t xx = x[0];

//     // Retrieve parameters
//     Double_t normP = par[0];
//     Double_t normN = par[1];
//     Double_t polyCoefficients[5]; 
    


//     // Fill polynomial coefficients
//     for (Int_t i = 0; i <=4; ++i) {
//         polyCoefficients[i] = par[i + 2];
//     }

//     // Calculate value using combination of histograms and polynomial background
//     val = normP * histP->Interpolate(xx) + normN * histN->Interpolate(xx);

//     // Add polynomial background
//     for (Int_t i = 0; i <= 4; ++i) {
//         val += polyCoefficients[i] * TMath::Power(xx, i);
//     }

//     return val;
// }


Double_t mc_p_fit(Double_t *x, Double_t *par) {

  Double_t normP = par[0];
  Double_t shiftP = par[1];
  Double_t val = 0.0;
  Double_t xx = x[0];
  val = normP * histP->Interpolate(xx+shiftP);
  return val;
}

Double_t mc_n_fit(Double_t *x, Double_t *par) {

  Double_t normN = par[0];
  Double_t val = 0.0;
  Double_t xx = x[0];
  val = normN * histN->Interpolate(xx);// why no shift here?
  return val;
}


// Utility function to shift every bin of a TH1D along the x-axis
  TH1D* shiftHistogramX(TH1D* originalHist, double shiftValue) {
    if (!originalHist) return nullptr;
    // Preserve the total number of entries in the histogram for further analysis
    double totalEntries = originalHist->GetEntries();
    // Create a new histogram with the same binning as the original
    TH1D *shiftedHist = (TH1D*)(originalHist->Clone("shiftedHist"));
    // Clear the contents of the cloned histogram
    shiftedHist->Reset();
    // Shift each bin
    for (int i = 1; i <= originalHist->GetNbinsX(); ++i) {
      // Calculate new bin center
      double oldBinCenter = originalHist->GetBinCenter(i);
      double oldBinContent = originalHist->GetBinContent(i);
      double oldBinError = originalHist->GetBinError(i);
      // Find the bin in the new histogram that corresponds to the new bin center
      int newBin = shiftedHist->FindBin(oldBinCenter + shiftValue);
      // Add the content and error to the new bin
      // Note: If multiple old bins shift into the same new bin, their contents and errors are added
      double newBinContent = shiftedHist->GetBinContent(newBin) + oldBinContent;
      double newBinError = sqrt(pow(shiftedHist->GetBinError(newBin), 2) + pow(oldBinError, 2));
      shiftedHist->SetBinContent(newBin, newBinContent);
      shiftedHist->SetBinError(newBin, newBinError);
    }
    // Restore the total number of entries
    shiftedHist->SetEntries(totalEntries);
    return shiftedHist;
  }


TH1D* GetResidualHistogram(TH1D* hist, TF1* fit) {
  // Create a new histogram for residuals
  TH1D* h_residual = new TH1D(Form("%s_residual", hist->GetName()), "", hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

  // Loop over bins and calculate residuals
  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
    double binCenter = hist->GetBinCenter(i);
    double binContent = hist->GetBinContent(i);
    double fitValue = fit->Eval(binCenter);
    double residual = binContent - fitValue;
    h_residual->SetBinContent(i, residual);
    h_residual->SetBinError(i, sqrt(binContent));
  }

  return h_residual;
}


TGraphErrors*  histogramToGraphErrors(TH1D *hist) {
    int numBins = hist->GetNbinsX();

    double x[numBins];
    double y[numBins];
    double ex[numBins]; // No x errors for simplicity
    double ey[numBins];

    // Fill arrays with histogram data
    for (int i = 0; i < numBins; ++i) {
        x[i] = hist->GetBinCenter(i + 1);
        y[i] = hist->GetBinContent(i + 1);
        ex[i] = 0; // No x errors for simplicity
        ey[i] = hist->GetBinError(i + 1);
    }

    // Create a TGraphErrors
    TGraphErrors *graph = new TGraphErrors(numBins, x, y, ex, ey);

    // // Create a canvas
    // TCanvas *canvas = new TCanvas("canvas", "Graph from Histogram", 800, 600);

    // // Draw the graph
    // graph->Draw("AP");

    // canvas->Update();
    // canvas->Draw();

    graph->SetTitle("");
    return graph;
}

void customizeGraph(TGraphErrors *graph, int markerStyle, int markerColor, double markersize) {
    int numPoints = graph->GetN();

    // Set marker style and color for each point
    for (int i = 0; i < numPoints; ++i) {
        graph->SetMarkerStyle(markerStyle);
        graph->SetMarkerColor(markerColor);
	graph->SetMarkerSize(markersize);
    }
}


TGraphErrors* createGraphFromFit(TH1D* hist, TF1* fitFunc) {
    int numBins = hist->GetNbinsX();

    double x[numBins];
    double y[numBins];
    double ex[numBins]; // No x errors for simplicity
    double ey[numBins];

    // Fill arrays with fit function data
    for (int i = 0; i < numBins; ++i) {
        double binCenter = hist->GetBinCenter(i + 1);
        x[i] = binCenter;
        y[i] = fitFunc->Eval(binCenter);
        ex[i] = 0; // No x errors for simplicity
        ey[i] = 0; // No y errors for simplicity
    }

    // Create and return a TGraphErrors
    TGraphErrors* graph = new TGraphErrors(numBins, x, y, ex, ey);
    graph->SetTitle("");
    return graph;
}

void adjustCanvas(TCanvas* canvas,
                  double leftMargin = 0.15, double rightMargin = 0.05, 
                  double bottomMargin = 0.15, double topMargin = 0.10) {
  // cout<<"left "<< leftMargin<<endl;
  // cout<<"right "<< rightMargin<<endl;
  // cout<<"bottom "<< bottomMargin<<endl;
  // cout<<"top "<< bottomMargin<<endl;
    // Set canvas margins
    canvas->SetLeftMargin(leftMargin);
    canvas->SetRightMargin(rightMargin);
    canvas->SetBottomMargin(bottomMargin);
    canvas->SetTopMargin(topMargin);
}

void adjustPad(TPad* pad,
                  double leftMargin = 0.15, double rightMargin = 0.05, 
                  double bottomMargin = 0.3, double topMargin = 0.10) {
  // cout<<"left "<< leftMargin<<endl;
  // cout<<"right "<< rightMargin<<endl;
  // cout<<"bottom "<< bottomMargin<<endl;
  // cout<<"top "<< bottomMargin<<endl;
    // Set canvas margins
    pad->SetLeftMargin(leftMargin);
    pad->SetRightMargin(rightMargin);
    pad->SetBottomMargin(bottomMargin);
    pad->SetTopMargin(topMargin);
}

void AdjustHistLabelOffset(TH1D* hist, double xoffset = 0.03, double yoffset= 0.02){
  // Adjust axis label positions
  hist->GetXaxis()->SetLabelOffset(xoffset); // Adjust X-axis label offset
  hist->GetYaxis()->SetLabelOffset(yoffset); // Adjust Y-axis label offset
}


void printParsedTitle(const std::string& title) {
  
  cout<<title<<endl;
  
  // Find the position of the first '{' character
  size_t pos = title.find('{');
    
  // Extract the y_axis:x_axis part
  std::string axes = title.substr(0, pos);
  //cout<<"axis = "<<axes<<endl;
    
  // Extract the cuts part and remove '{' and '}'
  std::string cuts = title.substr(pos + 1, title.size() - pos - 2);
  //cout<<"cuts = "<<cuts<<endl;
   
  //cout<<"broken up cuts"<<endl;

  // Split the cuts into individual cut expressions
  std::vector<std::string> cutList;
  std::stringstream ss(cuts);
  std::string cut;
  while (std::getline(ss, cut, '&')) {
    // Remove leading and trailing whitespace
    cut.erase(0, cut.find_first_not_of(" \t"));
    cut.erase(cut.find_last_not_of(" \t") + 1);
        
    // Ensure the cut is not empty before processing
    if (!cut.empty() && cut.front() == '&') {
      cut.erase(cut.begin());
    }
        
    if (!cut.empty()) {
      cutList.push_back(cut);
      // cout<<cut<<endl;
    }
  }

    
  // Create a new canvas
  TCanvas* cuts_canvas = new TCanvas("cuts_canvas", "Parsed Histogram Title", 1000, 600);
    
  // Create a TLatex object to draw the text
  TLatex latex;
  latex.SetTextSize(0.03);  // Adjust text size
  latex.SetTextAlign(13);   // Align text to top left
    
  // Draw the axes part
  latex.DrawLatex(0.1, 0.9, axes.c_str());
    
  // Draw each cut expression on a new line
  double yPos = 0.8;  // Start position for the first cut
  for (const auto& cut : cutList) {
    latex.DrawLatex(0.1, yPos, cut.c_str());
    yPos -= 0.03;  // Move down for the next cut
  }
    
  // Update the canvas
  cuts_canvas->Update();
    
  // // Optionally, save the canvas as an image
  // // canvas->SaveAs(Form("%s/global_cuts_%s.pdf", output_dir.c_str(),outputname));
}
