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
Double_t xMin = 0;
Double_t xMax = 1;

// Function Declarations 
Double_t linear(Double_t *x, Double_t *par);


//MAIN 
void fit_tgraph(){

  // load files 

  // set up arrays for the TGraphErrors 

  // std::vector<double> x = {0, 0.3, 0.5, 0.7, 0.9};
  // std::vector<double> y = {0.01, -0.527, -0.878, -1.241, -1.585};
  // std::vector<double> x_err = {0, 0, 0, 0, 0};
  // std::vector<double> y_err = {0.16/2, 0.18/2, 0.2/2, 0.14/2, 0.21/2};


  std::vector<double> x = {0.41, 0.39, 0.37, 0.35, 0.33, 0.359};
  std::vector<double> y = {-0.72876, -0.697977, -0.656032, -0.624027, -0.585491, -0.635539};
  std::vector<double> x_err = {0, 0, 0, 0, 0,0};
  // std::vector<double> y_err = {0.183669/2, 0.190514/2, 0.190213/2, 0.184015/2, 0.176586/2, 0.18396/2}; // sigmas
  std::vector<double> y_err = {0.003, 0.005, 0.004, 0.003, 0.0038, 0.0027}; // mean uncert

  
  // make graph
  TGraphErrors *gr = new TGraphErrors(x.size(), &x[0], &y[0], &x_err[0], &y_err[0]);

  // make linear fit 
  TF1 *fit_linear = new TF1("fit_linear", linear, xMin, xMax, 2);
  TFitResultPtr r = gr ->Fit("fit_linear", "MQR","",xMin,xMax);
  Int_t fitStatus = r;
  //if (fitStatus !=0) {cout<<"fit error" <<endl;}
  Double_t linear_offset = fit_linear->GetParameter(0);
  Double_t linear_slope = fit_linear->GetParameter(1);
  Double_t linear_offset_error = fit_linear->GetParError(0);
  Double_t linear_slope_error = fit_linear->GetParError(1);
  cout<<endl;
  cout<<"linear fit: "<<endl;
  cout<< "offset "<<linear_offset<<endl; 
  cout<<"slope "<< linear_slope<<endl;
  cout<<endl; 

  // set up canvas 
  TCanvas *c1 = new TCanvas("c1", "c1", 200, 10, 700, 500);

  // Set the color and style of the points
    gr->SetMarkerColor(kBlack);
    gr->SetMarkerStyle(21);

    // Set the color and width of the error bars
    gr->SetLineColor(kBlue);
    gr->SetLineWidth(2);

    gr->SetTitle("SBS4");
    gr->GetXaxis()->SetTitle("SBS Scale Field Param");
    gr->GetYaxis()->SetTitle("Proton dx");
    gr->Draw("AP");

    TF1 *data_line = new TF1("data_line", linear, xMin, xMax, 2);
    data_line ->SetParameters(-0.639,0);
    data_line ->SetLineColor(kMagenta);
    data_line ->SetLineStyle(2);
    data_line ->Draw("same");
   

      // Create a legend
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9); // Adjust the legend coordinates as needed
   // Add entries to the legend for the fit parameters
   // legend->AddEntry(fit_linear, "Linear Fit", "l");
    legend->AddEntry("", Form("y-intercept = %.4f +- %.4f ", linear_offset, linear_offset_error), "");
    legend->AddEntry("", Form("slope = %.4f +- %.4f ", linear_slope, linear_slope_error) , "");
    legend->Draw();
    c1 ->Update();

}

Double_t linear(Double_t *x, Double_t *par)
{
  Double_t fit =  par[0] + par[1]*x[0];
  return fit;
}
