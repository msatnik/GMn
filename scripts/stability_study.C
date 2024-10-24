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
#include <TArrayD.h>


//// Vectors to store the sliced histograms in. Declaring them outside so all the functions can see them. 
  std::vector<TH1D*> hist_vector_data; 
  std::vector<TH1D*> hist_vector_p; 
  std::vector<TH1D*> hist_vector_n; 

std::vector<TH1D*> hist_result_vector_p; 
std::vector<TH1D*> hist_result_vector_n; 

std::vector<TH1D*> hist_residual_vector; 



// these global histograms are going to temporarily hold the proton and neutron slices for when I call the fit function. 
TH1D *hist_p;
TH1D *hist_n;


 double xmin = -2.1; // -1.8
 double xmax = 1.4;//1;


 //// default values to be overwritten 
  TString DataFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_tight_2Dhistos.root";
  TString ProtonFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_tight_2Dhistos.root";
  TString NeutronFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_tight_2Dhistos.root";
  TString HistogramName = "hcal_dx__hcal_nsigdy";
  std::string AxisTitle ="hcal_nsigdy";
  TString note = "";
  std::string xMin_xMax_string = "(-4.1,-2.1),(-3.1,-1.1),(-2,0.1),(-1,1),(0,2),(1.1,3),(2.1,4)";


//// ranges for slices. These are what we want to change.
////****** dy *************************************************************** 
//std::string xMin_xMax_string_dy = "(-1,-0.5),(-0.5,0),(0,0.5),(0.5,1)";
//std::string xMin_xMax_string_dy = "(-1,-0.5),(-0.75,-0.25),(-0.5,0),(-0.25,0.25),(0,0.5),(0.25,0.75),(0.5,1)";
std::string xMin_xMax_string_dy = "(-1,1),(-0.9,0.9),(-0.8,0.8),(-0.7,0.7),(-0.6,0.6),(-0.5,0.5),(-0.4,0.4),(-0.3,0.3),(-0.2,0.2),(-0.1,0.1),(-0.05,0.05)";
////*****nsigdy***************************************************************
std::string xMin_xMax_string_nsigdy = "(-4.1,-2.1),(-3.1,-1.1),(-2,0.1),(-1,1),(0,2),(1.1,3),(2.1,4)";
///********vz**********************************************************************
//std::string xMin_xMax_string_vz = "(-0.07,0.07),(-0.065,0.065),(-0.06,0.06),(-0.05,0.05),(-0.04,0.04),(-0.03,0.03)";
std::string xMin_xMax_string_vz = "(-0.08,-0.07),(-0.07,-0.06),(-0.06,-0.05),(-0.05,-0.04),(-0.04,-0.03),(-0.03,-0.02),(-0.02,-0.01),(-0.01,0),(0,0.01),(0.01,0.02),(0.02,0.03),(0.03,0.04),(0.04,0.05),(0.05,0.06),(0.06,0.07),(0.07,0.08)";
//std::string xMin_xMax_string_vz = "(-0.08,-0.06),(-0.07,-0.05),(-0.06,-0.04),(-0.05,-0.03),(-0.04,-0.02),(-0.03,-0.01),(-0.02,0),(-0.01,0.01),(0,0.02),(0.01,0.03),(0.02,0.04),(0.03,0.05),(0.04,0.06),(0.05,0.07),(0.06,0.08)";
//std::string xMin_xMax_string_vz = "(-0.065,0.065)";
//// ******ps_e**********************************************************************
//std::string xMin_xMax_string_ps =  "(0.1,0.2),(0.12,0.22),(0.14,0.24),(0.16,0.26),(0.18,0.28),(0.2,0.3),(0.22,0.32),(0.24,0.34),(0.26,0.36),(0.28,0.38)";
std::string xMin_xMax_string_ps =  "(0.1,2),(0.12,2),(0.14,2),(0.16,2),(0.18,2),(0.2,2),(0.22,2),(0.24,2),(0.26,2),(0.28,2)";
//std::string xMin_xMax_string_ps ="(0.2,2)";
////******* hcal energy **************************************************************
//std::string xMin_xMax_string_hcal_e = "(0.015,0.025), (0.020,0.030), (0.025,0.035),(0.030,0.040)";
//std::string xMin_xMax_string_hcal_e = "(0.015,1), (0.020,1), (0.025,1),(0.030,1)";
std::string xMin_xMax_string_hcal_e = "(0.025,1)";
//// fiducial x
//// fiducial y
//// ******w2****************************************************************
//std::string xMin_xMax_string_w2 = "(0.4,0.6),(0.5,0.7),(0.6,0.8),(0.7,0.9),(0.8,1.0),(0.9,1.1),(1.0,1.2),(1.1,1.3) ";
//std::string xMin_xMax_string_w2 = "(0.4,1.2),(0.5,1.1),(0.6,1.0),(0.7,0.9)";
std::string xMin_xMax_string_w2 = "(0.78,0.98),(0.68,1.08),(0.58,1.18)";
//std::string xMin_xMax_string_w2 = "(0.4,0.5),(0.5,0.6),(0.6,0.7),(0.7,0.8),(0.8,0.9),(0.9,1.0),(1.0,1.1),(1.1,1.2),(1.2,1.3),(1.3,1.4),(1.4,1.5)";
// e_over_p
//// *** x_expected*******************************************
//std::string xMin_xMax_string_x_exp = "(-1.5,-1),(-1.25,-0.75),(-1,-0.5),(-0.75,-0.25),(-0.5,0),(-0.25,0.25),(0,0.5),(0.25,0.75),(0.5,1),(0.75,1.25)";
std::string xMin_xMax_string_x_exp = "(-0.8,0.0),(-0.9,0.1),(-1.0,0.2),(-1.1,0.3 )(-1.2,0.4) ,(-1.3,0.5),(-1.4,0.6)";
///****** nsigx_fid *******************************************
//std::string xMin_xMax_string_nsigx_fid = "(0,1),(0.5,1.5),(1,2),(1.5,2.5),(2,3),(2.5,3.5),(3.5,4.5),(4,5),(5,6),(5.5,6.5),(6,7),(6.5,7.5)";
//std::string xMin_xMax_string_nsigx_fid = "(0,7),(0.5,7),(1,7),(1.5,7),(2,7),(2.5,7),(3,7),(3.5,7),(4,7),(4.5,7),(5,7),(5.5,7),(6,7)";
std::string xMin_xMax_string_nsigx_fid = "(1.5,7),(2,7),(2.5,7),(3,7),(3.5,7)";
///******* y expected*************************************************************
//std::string xMin_xMax_string_y_exp ="(-0.6,-0.4),(-0.4,-0.2),(-0.2,0),(0,0.2),(0.2,0.4),(0.4,0.6),(0.6,0.8)";
//std::string xMin_xMax_string_y_exp = "(-0.35,0.4)";
std::string xMin_xMax_string_y_exp = "(-0.6,0.65),(-0.55,0.6),(-0.5,0.55),(-0.45,0.5),(-0.4,0.45),(-0.35,0.4),(-0.3,0.35),(-0.25,0.3),(-0.2,0.25),(-0.15,0.2),(-0.1,0.15)";

//********* nsigy_fid ***************************************************************
std::string xMin_xMax_string_nsigy_fid = "(0,3),(0.5,3),(1,3),(1.5,3),(2,3)";

// SBS4 30p 
TString DataFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_2Dhistos.root";
TString ProtonFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_2Dhistos.root";
TString NeutronFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_2Dhistos.root";
double xmin_sbs4_30p = -2.1; // 
double xmax_sbs4_30p = 1.4;//


// SBS4 50p 
TString DataFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_2Dhistos.root";
TString ProtonFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deep_2Dhistos.root";
TString NeutronFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deen_2Dhistos.root";
 double xmin_sbs4_50p = -2.5; // 
 double xmax_sbs4_50p = 1.4;//

// SBS8 70p
TString DataFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_2Dhistos_sept3.root";
TString ProtonFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deep_2Dhistos.root";
TString NeutronFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deen_2Dhistos.root";
 double xmin_sbs8_70p = -2.5; // need to adjust range
double xmax_sbs8_70p = 1.4;//




// Functions 
void SliceAndProjectHistogram_xMaxFixed(TH2D* hist2D, const std::vector<double>& xSlices, std::vector<TH1D*>& histVector, double xMax, std::string xAxisName, std::string yAxisName, std::string type);
Double_t mc_p_n_poly4_slice_fit(Double_t *x, Double_t *par);
Double_t mc_p_n_poly2_slice_fit(Double_t *x, Double_t *par);
TH1D* shiftHistogramX(TH1D* originalHist, double shiftValue);
Double_t poly4(Double_t *x, Double_t *par);
Double_t poly2(Double_t *x, Double_t *par);
TH1D* GetResidualHistogram(TH1D* hist, TF1* fit);
void adjustCanvas(TCanvas* canvas, double leftMargin = 0.15, double rightMargin = 0.05, double bottomMargin = 0.15, double topMargin = 0.10);
void customizeGraph(TGraphErrors *graph, int markerStyle, int markerColor, double markerSize, 
                    const std::string& graphTitle ="", const std::string& xAxisLabel="", const std::string& yAxisLabel="",
		    double TitleOffsetX = 1.4, double TitleOffsetY = 2, 
		    double LabelOffsetX = 0.01, double LabelOffsetY = 0.01);
void printParsedTitle(const std::string& title, TString outputlocation);
TH1D* sumHistogramsWithPolynomial(TH1D* h1, TH1D* h2, TF1* poly);
void SliceAndProjectHistogram_AroundMean(TH2D* hist2D, const std::vector<double>& xSlices, std::vector<TH1D*>& histVector, double xMean, std::string xAxisName, std::string yAxisName, std::string type);
void SliceAndProjectHistogram_xMinxMax(TH2D* hist2D, const std::vector<double>& xMin,const std::vector<double>& xMax, std::vector<TH1D*>& histVector, std::string xAxisName, std::string yAxisName, std::string type);
void parseStringToVectors(const std::string& input, std::vector<double>& xMin, std::vector<double>& xMax);
double calculateMean(const std::vector<double>& data);
double calculateStandardDeviation(const std::vector<double>& data);
double calculateWeightedMean(const std::vector<double>& values, const std::vector<double>& uncertainties);
double calculateWeightedStandardDeviation(const std::vector<double>& values, const std::vector<double>& uncertainties);
double CalculateWeightedMean(const std::vector<double>& data, const std::vector<double>& uncert);
double CalculateWeightedStDev(const std::vector<double>& data, const std::vector<double>& uncert);
double CalculateMean(const std::vector<double>& data);
double CalculateStDev(const std::vector<double>& data);

void stability_study(TString KineString="sbs4_30p", TString VarString = "dy"){ // main
  // bit of a test for now. Will need to make this more sophisticated in the future. 

  gStyle->SetNumberContours(255); 
  gStyle->SetOptStat(0110);

  std::vector<int> color_vector = {kRed,kBlue,kGreen,kMagenta,kCyan,kYellow};
  

  // Taking user input on the kinematic and what variable we want to study. Going to try keeping it all here instead of config files. 

  if (VarString == "nsigdy"){
    HistogramName = "hcal_dx__hcal_nsigdy";
    AxisTitle ="hcal_nsigdy";
    xMin_xMax_string = xMin_xMax_string_nsigdy;
  }else if (VarString == "dy"){
    HistogramName = "hcal_dx__hcal_dy";
    AxisTitle ="hcal_dy";
    xMin_xMax_string = xMin_xMax_string_dy;
  } else if (VarString == "vz"){
    HistogramName = "hcal_dx__tr_vz";
    AxisTitle ="track vz";
    xMin_xMax_string = xMin_xMax_string_vz;
  }else if (VarString == "ps"){
    HistogramName = "hcal_dx__ps_e";
    AxisTitle ="preshower energy";
    xMin_xMax_string = xMin_xMax_string_ps;
  }else if (VarString == "hcal_e"){
    HistogramName = "hcal_dx__hcal_e";
    AxisTitle ="hcal energy";
    xMin_xMax_string = xMin_xMax_string_hcal_e;
  }else if (VarString == "w2"){
    HistogramName = "hcal_dx__W2";
    AxisTitle ="W^{2}";
    xMin_xMax_string = xMin_xMax_string_w2;
  }else if (VarString == "x_exp"){
    HistogramName = "hcal_dx__hcal_x_exp";
    AxisTitle ="x expected";
    xMin_xMax_string = xMin_xMax_string_x_exp;
  }else if (VarString == "nsigx_fid"){
    HistogramName = "hcal_dx__nsigx_fid";
    AxisTitle ="nsigx_fid";
    xMin_xMax_string = xMin_xMax_string_nsigx_fid;
  }else if (VarString == "y_exp"){
    HistogramName = "hcal_dx__hcal_y_exp";
    AxisTitle ="y expected";
    xMin_xMax_string = xMin_xMax_string_y_exp;
  }else if (VarString == "nsigy_fid"){
    HistogramName = "hcal_dx__nsigy_fid";
    AxisTitle ="nsigy_fid";
    xMin_xMax_string = xMin_xMax_string_nsigy_fid;
  } else {
    std::cout<<"Error in the Varibale String "<<std::endl;
    return;
  }
  

  if (KineString == "sbs4_30p"){
    DataFileString = DataFileString_sbs4_30p;
    ProtonFileString =ProtonFileString_sbs4_30p;
    NeutronFileString = NeutronFileString_sbs4_30p;
    xmin = xmin_sbs4_30p;
    xmax = xmax_sbs4_30p;
  }else if (KineString == "sbs4_50p"){
    DataFileString = DataFileString_sbs4_50p;
    ProtonFileString =ProtonFileString_sbs4_50p;
    NeutronFileString = NeutronFileString_sbs4_50p;
    xmin = xmin_sbs4_50p;
    xmax = xmax_sbs4_50p;
  }else if (KineString == "sbs8_70p"){
    DataFileString = DataFileString_sbs8_70p;
    ProtonFileString =ProtonFileString_sbs8_70p;
    NeutronFileString = NeutronFileString_sbs8_70p;
    xmin = xmin_sbs8_70p;
    xmax = xmax_sbs8_70p;
  }  else {
    std::cout<<"Error with kinematic setting"<<std::endl;
    return;
  }


  cout<<"What we're running: "<<KineString<< ", "<<VarString<<endl;
  cout<<DataFileString<<endl;
  cout<<ProtonFileString<<endl;
  cout<<NeutronFileString<<endl;
  cout<<xmin<< ", "<<xmax<<endl;
  cout<<HistogramName<<endl;
  cout<<AxisTitle<<endl;
  cout<<xMin_xMax_string<<endl;
 
      
      
  
  TFile *f1 = TFile::Open(DataFileString); // data
  TFile *f2 = TFile::Open(ProtonFileString); //proton
  TFile *f3 = TFile::Open(NeutronFileString); //neutron



  //// set location and name for root output file 
  TString outputfilelocation="../output/stability/"+KineString+"/"+VarString;
  TString outputfilename = outputfilelocation +"/"+ KineString+ "_" +VarString+".root";
  
  // TFile *fout = new TFile("../output/stability/testoutput.root","RECREATE");

  TFile *fout = new TFile(outputfilename,"RECREATE");
  cout<<"writing to file: "<< outputfilename <<endl;


  // Load Histograms
  TH2D *hist_data_orig = (TH2D*)f1->Get(HistogramName);
  TH2D *hist_proton_orig = (TH2D*)f2->Get(HistogramName);
  TH2D *hist_neutron_orig = (TH2D*)f3->Get(HistogramName);

  // Make clones of the histograms
  TH2D *hist_2D_data = (TH2D*)hist_data_orig->Clone("hist_2D_data");
  TH2D *hist_2D_proton = (TH2D*)hist_proton_orig->Clone("hist_2D_proton");
  TH2D *hist_2D_neutron = (TH2D*)hist_neutron_orig->Clone("hist_2D_neutron");

  

 // Set up the slices for the Y projections. 
  //std::vector<double> xSlices = {0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25,3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0};
  //std::vector<double> xSlices = {5,4.5,4,3.5,3,2.5,2,1.5,1};
  //double mean = 0;

  std::vector<double> xMinSlices;
 std:vector<double> xMaxSlices;
  parseStringToVectors(xMin_xMax_string, xMinSlices, xMaxSlices);
  
  // Output the results
  std::cout<<"xMin: ";
    for (double val : xMinSlices) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    std::cout << "xMax: ";
    for (double val : xMaxSlices) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

  //std::vector<double> xMinSlices ={-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0};
  // std::vector<double> xMaxSlices={-2.0,-1.0, 0.0, 1.0,2.0,3.0,4.0};


  //// Vectors to store the sliced histograms in. 
  std::vector<TH1D*> hist_vector_data; 
  std::vector<TH1D*> hist_vector_p; 
  std::vector<TH1D*> hist_vector_n;

 
  cout<<endl<<"Slicing Proton"<<endl;
  SliceAndProjectHistogram_xMinxMax(hist_2D_proton, xMinSlices, xMaxSlices, hist_vector_p, AxisTitle,"dx","p");
  cout<<endl<<"Slicing Neutron"<<endl;;
  SliceAndProjectHistogram_xMinxMax(hist_2D_neutron, xMinSlices, xMaxSlices, hist_vector_n, AxisTitle,"dx","n");
   cout<<endl<<"Slicing Data"<<endl;;
  SliceAndProjectHistogram_xMinxMax(hist_2D_data, xMinSlices,xMaxSlices, hist_vector_data, AxisTitle,"dx","data");

 // // use the function to make the slices and projections
 //  SliceAndProjectHistogram_AroundMean(hist_2D_proton, xSlices, hist_vector_p, mean, "hcal_nsigdy","dx","p");
 //  SliceAndProjectHistogram_AroundMean(hist_2D_neutron, xSlices, hist_vector_n, mean, "hcal_nsigdy","dx","n");
 //  SliceAndProjectHistogram_AroundMean(hist_2D_data, xSlices, hist_vector_data, mean, "hcal_nsigdy","dx","data");
 

  //   double initialParameters[9] = {1,1,0,0,1,1,1,1,1};
  //    TF1 *overall_fit = new TF1("overall_fit", mc_p_n_poly4_slice_fit, xmin, xmax, 9); //

  //    TH1D *hist_1d_test_data = hist_vector_data[9];
  //    //// global histograms that the fit function will use. Need to be set up properly before fit function is called. 
  //    hist_p = hist_vector_p[9];
  //    hist_n = hist_vector_n[9];

  //   // Set initial parameters for the fit function
  //      overall_fit->SetParameters(initialParameters); // Define initial parameters
  //      overall_fit->SetNpx(500);
    
  // // // Fit combined histogram with custom fit function
  //      hist_1d_test_data->GetXaxis() ->SetRangUser(xmin, xmax);
  //      hist_1d_test_data->Fit(overall_fit,"Q");

  std::vector<TF1*> fit_vector;

  std::vector<TH1D*> hist_result_p_vector;
  std::vector<TH1D*> hist_result_n_vector;

  std::vector<double> scale_p_vector;
  std::vector<double> scale_n_vector;
  std::vector<double> shift_p_vector;
  std::vector<double> shift_n_vector;

  std::vector<double> scale_p_err_vector;
  std::vector<double> scale_n_err_vector;
  std::vector<double> shift_p_err_vector;
  std::vector<double> shift_n_err_vector;
  

  std::vector<double> ChiSq_vector;
  std::vector<double> ndf_vector;

  std::vector<double> Rsf_vector;
  std::vector<double> Rsf_err_vector;
 
  std::vector<std::vector<double>> poly_result_vector_of_vectors;
  std::vector<std::vector<double>> poly_result_err_vector_of_vectors;

  fit_vector.clear();
  hist_result_p_vector.clear();
  hist_result_n_vector.clear();
  scale_p_vector.clear();
  scale_n_vector.clear();
  shift_p_vector.clear();
  shift_n_vector.clear();
  scale_p_err_vector.clear();
  scale_n_err_vector.clear();
  shift_p_err_vector.clear();
  shift_n_err_vector.clear();
  ChiSq_vector.clear();
  ndf_vector.clear();
  Rsf_vector.clear();
  Rsf_err_vector.clear();
  poly_result_vector_of_vectors.clear();
  poly_result_err_vector_of_vectors.clear();

  double initialParameters[7]={1,1,-0.02,-0.09,1,1,-1};
  // initialParameters = {0};
      
  TH1D *hist_data;
  /// Fit each slice 
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {
      /// set up the global histograms that the fit function is going to use
      // this uses clones of the histograms
      hist_p  =(TH1D*)hist_vector_p[sliceid]->Clone("hist_p");
      hist_n  =(TH1D*)hist_vector_n[sliceid]->Clone("hist_n");
      hist_data  =(TH1D*)hist_vector_data[sliceid]->Clone("hist_data");

      //// this uses the histograms themselves: you can't uncouple the fit from them easily 
      //  hist_p  =hist_vector_p[sliceid];
      //  hist_n  =hist_vector_n[sliceid];
      // hist_data  =hist_vector_data[sliceid];

	 
      TF1 *Fit = new TF1(Form("overall_fit_%i_",sliceid),mc_p_n_poly2_slice_fit, xmin, xmax, 7);
      Fit->SetParameters(initialParameters);
       // set parameter limits. SetParLimits -> (par#, min, max)
      Fit->SetParLimits(0, 0, 2000); // scale_p greater than 0
      Fit->SetParLimits(1, 0,2000); // scale_n greater than 0
      Fit->SetParLimits(2, -0.10,0.10); // shift_p less than +- 10cm
      Fit->SetParLimits(3, -0.10,0.10); // shift_n less than +- 10cm
      Fit->SetParLimits(6,-10000000,-0.0000000001); // x^2 term negative to force downward concavity 
      Fit->SetNpx(500);
      hist_data->GetXaxis() ->SetRangeUser(xmin, xmax);
      hist_data->Fit(Fit,"R Q");

      // retrieve fit results 
      double scale_p  = Fit ->GetParameter(0);
      double scale_p_err = Fit ->GetParError(0);
      double scale_n  = Fit ->GetParameter(1);
      double scale_n_err = Fit ->GetParError(1);

      double shift_p= Fit ->GetParameter(2);
      double shift_p_err= Fit ->GetParError(2); 
      double shift_n = Fit ->GetParameter(3);
      double shift_n_err = Fit ->GetParError(3);
	  
      double ChiSq= Fit->GetChisquare();
      double ndf = Fit->GetNDF();

      //cout<<"shift_p in the loop: "<<shift_p<< " +/- "<<shift_p_err<<endl;
      //cout<<"shift_n in the loop: "<<shift_n<< " +/- "<<shift_n_err<<endl;

      std::vector<double> poly_result;
      std::vector<double> poly_result_err;
      for (int i =0 ; i < 3; i++)
	{
	  poly_result.push_back( Fit->GetParameter(4+i) );
	  poly_result_err.push_back(Fit->GetParError(4+i) );
	}

      //compute results
      double Rsf = scale_n/scale_p;
      double Rsf_err = Rsf * sqrt( pow( (scale_n_err / scale_n), 2) + pow( (scale_p_err / scale_p),2) ); //just adding the uncert from the fit parameters in quadrature for now. 

      cout<<"sliceid: "<<sliceid<< " ratio: scale_n / scale_p = "<< scale_n <<" / " <<scale_p <<" = "<< Rsf <<" +/- "<< Rsf_err<<endl;
      cout<<"sliceid: "<<sliceid<<" shift p: "<<shift_p<<" +/- "<<shift_p_err<<" shift_n: "<<shift_n<<"+/-"<<shift_n_err<<endl;
      cout<<endl;

      //// save results into the vectors 
      //fit_vector.push_back(Fit);
      scale_p_vector.push_back(scale_p);
      scale_n_vector.push_back(scale_n);
      scale_p_err_vector.push_back(scale_p_err);
      scale_n_err_vector.push_back(scale_n_err);
      shift_p_vector.push_back(shift_p);
      shift_n_vector.push_back(shift_n);
      shift_p_err_vector.push_back(shift_p_err);
      shift_n_err_vector.push_back(shift_n_err);
      ChiSq_vector.push_back(ChiSq);
      ndf_vector.push_back(ndf);
      Rsf_vector.push_back(Rsf);
      Rsf_err_vector.push_back(Rsf_err);
      // poly_result_vector_of_arrays.push_back(poly_result);
      poly_result_vector_of_vectors.push_back(poly_result);
      poly_result_err_vector_of_vectors.push_back(poly_result_err);
	  
      delete Fit; //
    }// end loop over slices


     // // Print the contents of the vector of vectors
     // for (const auto& vec : poly_result_vector_of_vectors) {
     //   for (double val : vec) {
     // 	 std::cout << val << " ";
     //   }
     //   std::cout << std::endl;
     // }


  cout<<"---------------------------------------------------------------------------------------------"<<endl;
  double Rsf_mean = CalculateMean(Rsf_vector);
  double Rsf_stdev = CalculateStDev(Rsf_vector);
  double Rsf_mean_w = CalculateWeightedMean(Rsf_vector,Rsf_err_vector);
  double Rsf_stdev_w = CalculateWeightedStDev(Rsf_vector,Rsf_err_vector);
  cout<<"Rsf mean = " << Rsf_mean<<endl;
  cout<<"Rsf StDev = " << Rsf_stdev<<endl;
  cout<<"weighted Rsf mean = "<<Rsf_mean_w<<endl;
  cout<<"weighted Rsf StDev = "<<Rsf_stdev_w<<endl;
  cout<<"---------------------------------------------------------------------------------------------"<<endl;
  cout<<endl;


 
  
  
  TH1D *hist_temp_p;
  TH1D *hist_temp_n;

  // loop over the slices to create shifted and scaled versions of the mc histograms to plot
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {	
      /// shift the histograms. Function returns a clone. 
      hist_temp_p = shiftHistogramX(hist_vector_p[sliceid], shift_p_vector[sliceid] );
      hist_temp_n = shiftHistogramX(hist_vector_n[sliceid], shift_n_vector[sliceid] );

      /// scale the histograms 
      hist_temp_p -> Scale(scale_p_vector[sliceid]);
      hist_temp_n -> Scale(scale_n_vector[sliceid]);
	 
      // save histograms to the global vectors 
      hist_result_vector_p.push_back(hist_temp_p);
      hist_result_vector_n.push_back(hist_temp_n);	
    }// end loop over slices


  std::vector<TF1*> poly_fit_result;
  // loop over slices and make fits to plot the polynomial background result
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {
      // get the fit results from the vector of vectors
      std::vector<double> poly_params_vector;
      poly_params_vector = poly_result_vector_of_vectors[sliceid];

      // convert the vector into an array that TF1 can use
      int size = poly_params_vector.size();
      double poly_params_array[size]; 
      for (int i  = 0 ; i < size; i++){
	poly_params_array[i] = poly_params_vector[i];
      }

      TF1 *fit = new TF1(Form("poly2_%i_",sliceid),poly2, xmin, xmax, 3);
      fit->SetParameters(poly_params_array); 
      fit->SetNpx(500);
      poly_fit_result.push_back(fit);
    }// end loop over slices


     //// Plot Rsf
     //// Make arrays that TGraphErrors can use 
  double x[xMinSlices.size()];
  double y[xMinSlices.size()];
  double x_err[xMinSlices.size()];
  double y_err[xMinSlices.size()];
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {
      double xCenter = (xMaxSlices[sliceid] + xMinSlices[sliceid] )/2;
      double xWidth = xMaxSlices[sliceid] - xMinSlices[sliceid];
     
      x[sliceid] = xCenter;
      y[sliceid] = Rsf_vector[sliceid];
      x_err[sliceid] = xWidth/2;
      y_err[sliceid] = Rsf_err_vector[sliceid];

      //cout<< "sliceid = "<<sliceid<<", xCenter = "<<x[sliceid]<<", x= "<<x[sliceid]<<", y = "<<y[sliceid]<<", x_err= "<<x_err[sliceid]<<", y_err = "<<y_err[sliceid]<<endl;
    }    
  TGraphErrors *Rsf_graph = new TGraphErrors(xMinSlices.size(), x, y, x_err, y_err);
  customizeGraph(Rsf_graph, 33, kBlue, 3,"",AxisTitle,"Rsf");

  
  
  //// canvas
  TCanvas *graphcanvas = new TCanvas("graphcanvas","graphcanvas",800,600);
  graphcanvas->SetGrid();
  adjustCanvas(graphcanvas);
  Rsf_graph->Draw("AP");
  graphcanvas->Update();
  graphcanvas->SaveAs(Form("%s/Rsf_xCenter.pdf",outputfilelocation.Data() ) );


  //// Plot nEntries
     //// Make arrays that TGraphErrors can use 
  double x_n[xMinSlices.size()];
  double y_n[xMinSlices.size()];
  double x_n_err[xMinSlices.size()];
  double y_n_err[xMinSlices.size()];
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {
      double xCenter_n = (xMaxSlices[sliceid] + xMinSlices[sliceid] )/2;
      double xWidth_n = (xMaxSlices[sliceid] - xMinSlices[sliceid]);
     
      x_n[sliceid] = xCenter_n;
      y_n[sliceid] = hist_vector_data[sliceid]->GetEntries();
      x_n_err[sliceid] = xWidth_n/2;
      y_n_err[sliceid] = 0;
    }    
  TGraphErrors *nEntries_graph = new TGraphErrors(xMinSlices.size(), x_n, y_n, x_n_err, y_n_err);
  customizeGraph(nEntries_graph, 33, kBlue, 3,"",AxisTitle,"nEntries");
  //// canvas
  TCanvas *nEntriesCanvas = new TCanvas("nEntriesCanvas","nEntriesCanvas",800,600);  nEntriesCanvas->SetGrid();
  adjustCanvas(nEntriesCanvas);
  nEntries_graph->Draw("AP");
  nEntriesCanvas->Update();
  nEntriesCanvas->SaveAs(Form("%s/nEntries.pdf",outputfilelocation.Data()));


   double x_max[xMinSlices.size()];
  double y_max[xMinSlices.size()];
  double x_max_err[xMinSlices.size()];
  double y_max_err[xMinSlices.size()];
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {
      x_max[sliceid] = xMaxSlices[sliceid];
      y_max[sliceid] = Rsf_vector[sliceid];
      x_max_err[sliceid] = 0;
      y_max_err[sliceid] = Rsf_err_vector[sliceid];
    }    
  TGraphErrors *Rsf_xMax_graph = new TGraphErrors(xMinSlices.size(), x_max, y_max, x_max_err, y_max_err);
  customizeGraph(Rsf_xMax_graph, 33, kBlue, 3,"",AxisTitle,"Rsf");
  
 //// canvas
  TCanvas *Rsf_xMaxCanvas = new TCanvas("Rsf_xMaxCanvas","Rsf_xMaxCanvas",800,600);
  Rsf_xMaxCanvas->SetGrid();
  adjustCanvas(Rsf_xMaxCanvas);
  Rsf_xMax_graph->Draw("AP");
  Rsf_xMaxCanvas->Update();
  Rsf_xMaxCanvas->SaveAs(Form("%s/Rsf_xMax.pdf",outputfilelocation.Data()));
  

   double x_m[xMinSlices.size()];
  double y_m[xMinSlices.size()];
  double x_m_err[xMinSlices.size()];
  double y_m_err[xMinSlices.size()];
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {
      x_m[sliceid] = xMinSlices[sliceid];
      y_m[sliceid] = Rsf_vector[sliceid];
      x_m_err[sliceid] = 0;
      y_m_err[sliceid] = Rsf_err_vector[sliceid];
    }    
  TGraphErrors *Rsf_xMin_graph = new TGraphErrors(xMinSlices.size(), x_m, y_m, x_m_err, y_m_err);
  customizeGraph(Rsf_xMin_graph, 33, kBlue, 3,"",AxisTitle,"Rsf");
  
 //// canvas
  TCanvas *Rsf_xMinCanvas = new TCanvas("Rsf_xMinCanvas","Rsf_xMinCanvas",800,600);
  Rsf_xMinCanvas->SetGrid();
  adjustCanvas(Rsf_xMinCanvas);
  Rsf_xMin_graph->Draw("AP");
  Rsf_xMinCanvas->Update();
  Rsf_xMinCanvas->SaveAs(Form("%s/Rsf_xMin.pdf",outputfilelocation.Data()));

  

  int nHist = hist_vector_data.size();
  int nCols = 4;
  int nRows = (nHist + nCols - 1) / nCols;


  std::vector<TH1D*> overall_fit_as_histogram;

  TCanvas* fits_canvas = new TCanvas("fits_canvas", "fits_canvas", 1000, 600);
  if (hist_vector_data.size()!=1)
    {
      fits_canvas->Divide(nCols, nRows);
    }
  for (int i = 0; i < nHist; ++i) {
    fits_canvas->cd(i + 1);
    fits_canvas->SetGrid();
    hist_vector_data[i]->GetXaxis() ->SetRangeUser(xmin, xmax);
    hist_vector_data[i]->Draw();
      
    //// make a histogram that is the sum of the poly fit and the scaled mc histograms. 
    TH1D* sum_histo = sumHistogramsWithPolynomial(hist_result_vector_p[i],hist_result_vector_n[i] , poly_fit_result[i]);
    overall_fit_as_histogram.push_back(sum_histo);

    // // Create and customize the legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    // legend->AddEntry(hist_vector_data[i], hist_vector_data[i]->GetName(), "l");
    legend->AddEntry("", Form("R= %.4f +/- %.4f ", Rsf_vector[i],Rsf_err_vector[i]), "");
    legend->AddEntry("", Form("#chi^{2}/ndf = %.2f / %.0f  ", ChiSq_vector[i] ,ndf_vector[i]), "");
    legend->Draw();

    poly_fit_result[i]->SetLineColor(kCyan);
    poly_fit_result[i]->Draw("same");
    hist_result_vector_p[i] ->SetLineColor(kGreen);
    hist_result_vector_p[i]->Draw("same");
    hist_result_vector_n[i] ->SetLineColor(kMagenta);
    hist_result_vector_n[i]->Draw("same");
    sum_histo ->SetLineColor(kRed);
    sum_histo->Draw("same");
     
  }// end loop over slices
  fits_canvas->Update();
   fits_canvas->SaveAs(Form("%s/fitted_slices.pdf",outputfilelocation.Data()));
  

  std::vector<TH1D*> hist_residual_vector;
 
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {
      TH1D* hist_residual = (TH1D*)hist_vector_data[sliceid]->Clone("hist_residual");  
      hist_residual  ->Add(overall_fit_as_histogram[sliceid], -1);
      hist_residual->GetXaxis() ->SetRangeUser(xmin, xmax);
      hist_residual_vector.push_back(hist_residual);
    }// end loop over slices

   
  TCanvas* resid_canvas = new TCanvas("resid_canvas", "resid_canvas", 1000, 600);
  if (hist_vector_data.size()!=1)
    {
      resid_canvas->Divide(nCols, nRows);
    }
  for (int i = 0; i < nHist; ++i) {
    resid_canvas->cd(i + 1);
    resid_canvas->SetGrid();
    hist_residual_vector[i]->Draw("E sames");
  }
  resid_canvas->Update();
  resid_canvas->SaveAs(Form("%s/residuals.pdf",outputfilelocation.Data() ));

     
  // Create a canvas to draw the histogram
  TCanvas* cutcanvas = new TCanvas("cutcanvas", "visulaized cuts", 800, 600);

  // Draw the TH2D histogram
  hist_2D_data->Draw("COLZ");

  // Draw vertical lines at each x value in xSlices
  for (size_t i = 0; i < xMinSlices.size(); ++i) {
    int color = color_vector[i%(color_vector.size()-1)];
    // cout<<" color vec index ="<<i%(color_vector.size()-1)<<endl;
    double xLeft = xMinSlices[i];
    double xRight = xMaxSlices[i];
    TLine* lineLeft = new TLine(xLeft, hist_2D_data->GetYaxis()->GetXmin(), xLeft, hist_2D_data->GetYaxis()->GetXmax());
    TLine* lineRight = new TLine(xRight, hist_2D_data->GetYaxis()->GetXmin(),xRight, hist_2D_data->GetYaxis()->GetXmax());
    lineLeft->SetLineColor(color);  // Set line color 
    lineLeft->SetLineWidth(2);     // Set line width
    lineLeft->Draw("SAME");
    lineRight->SetLineColor(color);  // Set line color 
    lineRight->SetLineWidth(2);     // Set line width
    lineRight->Draw("SAME");
    //cout<<"sliceid = "<<i<<" xLeft = "<<xLeft<<" xRight = "<<xRight <<" color = "<<color<<endl;
  }

  // // Save the canvas to a file or display it
  cutcanvas->SaveAs(Form("%s/visualized_cuts.pdf",outputfilelocation.Data() ));


  
  // Create a projection of the TH2D histogram onto the x-axis
  TH1D* hist1D_data = hist_2D_data->ProjectionX("dataProjX");
  TH1D* hist1D_p = hist_2D_proton->ProjectionX("protonProjX");
  TH1D* hist1D_n = hist_2D_neutron->ProjectionX("neutronProjX");

  cout<<endl;
  cout<<"Projection onto x-axis"<<endl;
  cout<<"data:"<<endl;
  cout<< "Mean = "<<hist1D_data->GetMean()<<", StdDev = "<<hist1D_data->GetStdDev()<<endl;
  cout<<"proton:"<<endl;
  cout<< "Mean = "<<hist1D_p->GetMean()<<", StdDev = "<<hist1D_p->GetStdDev()<<endl;
  cout<<"neutron:"<<endl;
  cout<< "Mean = "<<hist1D_n->GetMean()<<", StdDev = "<<hist1D_n->GetStdDev()<<endl;
  cout<<endl;

  // Create a projection of the TH2D histogram onto the y-axis
  TH1D* hist1D_data_y = hist_2D_data->ProjectionY("dataProjY");
  TH1D* hist1D_p_y = hist_2D_proton->ProjectionY("protonProjY");
  TH1D* hist1D_n_y = hist_2D_neutron->ProjectionY("neutronProjY");

  cout<<"Projection onto Y-axis"<<endl;
  cout<<"data:"<<endl;
  cout<< "Mean = "<<hist1D_data_y->GetMean()<<", StdDev = "<<hist1D_data_y->GetStdDev()<<endl;
  cout<<"proton:"<<endl;
  cout<< "Mean = "<<hist1D_p_y->GetMean()<<", StdDev = "<<hist1D_p_y->GetStdDev()<<endl;
  cout<<"neutron:"<<endl;
  cout<< "Mean = "<<hist1D_n_y->GetMean()<<", StdDev = "<<hist1D_n_y->GetStdDev()<<endl;
  cout<<endl;

  

  // Create a canvas to draw the histogram
  TCanvas* cutcanvas1D = new TCanvas("cutcanvas1D", "X Projection with Vertical Lines", 800, 600);
  cutcanvas1D->Divide(1,3);

  // Draw the TH1D histogram
  cutcanvas1D->cd(1);
  //hist1D_data->GetXaxis() ->SetRangeUser(0, 0.1);
  hist1D_data->Draw();
  // Draw vertical lines at each x value in xSlices
  for (size_t i = 0; i < xMinSlices.size(); ++i) {
    int color = color_vector[i%(color_vector.size()-1)];
    double xLeft = xMinSlices[i];
    double xRight = xMaxSlices[i];
    TLine* lineLeft = new TLine(xLeft, hist1D_data->GetMinimum(), xLeft, hist1D_data->GetMaximum());
    TLine* lineRight = new TLine(xRight, hist1D_data->GetMinimum(), xRight, hist1D_data->GetMaximum());
    lineLeft->SetLineColor(color);  // Set line color (e.g., red)
    lineLeft->SetLineWidth(2);     // Set line width
    lineLeft->Draw("SAME");
    lineRight->SetLineColor(color);  // Set line color (e.g., red)
    lineRight->SetLineWidth(2);     // Set line width
    lineRight->Draw("SAME");
  }

  // Draw the TH1D histogram
  cutcanvas1D->cd(2);
  //hist1D_p->GetXaxis() ->SetRangeUser(0, 0.1);
  hist1D_p->Draw("hist");
  // Draw vertical lines at each x value in xSlices
  for (size_t i = 0; i < xMinSlices.size(); ++i) {
 int color = color_vector[i%(color_vector.size()-1)];
    double xLeft = xMinSlices[i];
    double xRight = xMaxSlices[i];
    TLine* lineLeft = new TLine(xLeft, hist1D_p->GetMinimum(), xLeft, hist1D_p->GetMaximum());
    TLine* lineRight = new TLine(xRight, hist1D_p->GetMinimum(), xRight, hist1D_p->GetMaximum());
    lineLeft->SetLineColor(color);  // Set line color (e.g., red)
    lineLeft->SetLineWidth(2);     // Set line width
    lineLeft->Draw("SAME");
    lineRight->SetLineColor(color);  // Set line color (e.g., red)
    lineRight->SetLineWidth(2);     // Set line width
    lineRight->Draw("SAME"); 
  }
  
  // Draw the TH1D histogram
  cutcanvas1D->cd(3);
  //hist1D_n->GetXaxis() ->SetRangeUser(0, 0.1);
  hist1D_n->Draw("hist");
  // Draw vertical lines at each x value in xSlices
  for (size_t i = 0; i < xMinSlices.size(); ++i) { 
     int color = color_vector[i%(color_vector.size()-1)];
    double xLeft = xMinSlices[i];
    double xRight = xMaxSlices[i];
    TLine* lineLeft = new TLine(xLeft, hist1D_n->GetMinimum(), xLeft, hist1D_n->GetMaximum());
    TLine* lineRight = new TLine(xRight, hist1D_n->GetMinimum(), xRight, hist1D_n->GetMaximum());
    lineLeft->SetLineColor(color);  // Set line color (e.g., red)
    lineLeft->SetLineWidth(2);     // Set line width
    lineLeft->Draw("SAME");
    lineRight->SetLineColor(color);  // Set line color (e.g., red)
    lineRight->SetLineWidth(2);     // Set line width
    lineRight->Draw("SAME"); 
  }

  cutcanvas1D->SaveAs(Form("%s/visulaized_cuts_1D.pdf",outputfilelocation.Data() ));
  

  // // TCanvas* testcanvas1  = new TCanvas("testcanvas1", "testcanvas1", 800, 600);
  // // TH1D* sumHistoTest = sumHistogramsWithPolynomial(hist_result_vector_p[0],hist_result_vector_n[0] , poly_fit_result[0]);
  // // hist_vector_data[0]->Draw();
  // // sumHistoTest ->Draw("same");

  
 
  // //// extract the histogram title and print it
  std::string title = hist_2D_data->GetTitle();
  cout<<"data title"<<endl;
  printParsedTitle(title,outputfilelocation);
  std::string proton_title = hist_2D_proton->GetTitle();
  cout<<"proton title: "<<endl;
  cout<< proton_title<<endl;
  std::string neutron_title =hist_2D_neutron->GetTitle();
  cout<<"neutron title: "<<endl;
  cout<< neutron_title<<endl;


  fout->Write(); 

  f1->Close();
  f2->Close();
  f3->Close();

}// End Main



TH1D* sumHistogramsWithPolynomial(TH1D* h1, TH1D* h2, TF1* poly) {
    if (h1->GetNbinsX() != h2->GetNbinsX() || 
        h1->GetXaxis()->GetXmin() != h2->GetXaxis()->GetXmin() ||
        h1->GetXaxis()->GetXmax() != h2->GetXaxis()->GetXmax()) {
        std::cerr << "Histograms must have the same binning and range!" << std::endl;
        return nullptr;
    }

    // Create a new histogram for the sum
    TH1D *h_sum = (TH1D*)h1->Clone("h_sum");
    h_sum->SetTitle("Sum of Histograms and Polynomial");
    h_sum->Reset();


    // Loop over bins and add the content of h1, h2, and the polynomial
    for (int i = 1; i <= h_sum->GetNbinsX(); ++i) {
        double bin_center = h_sum->GetBinCenter(i);
        double content = h1->GetBinContent(i) + h2->GetBinContent(i) + poly->Eval(bin_center);
        h_sum->SetBinContent(i, content);
    }

    return h_sum;
}// end sumHistogramsWithPolynomial


void SliceAndProjectHistogram_xMinxMax(TH2D* hist2D, const std::vector<double>& xMinimum,const std::vector<double>& xMaximum, std::vector<TH1D*>& histVector, std::string xAxisName, std::string yAxisName, std::string type) {
    // Clear the vector to ensure it's empty before filling
    histVector.clear();

    if (xMinimum.size()!=xMaximum.size())
      {
	cout<<"error in the xMin and xMax ranges. Not same size."<<endl;
      }

    cout<< "xMinimum.size() "<<xMinimum.size() <<" , xMaximum.size() "<<xMaximum.size()<<endl;
  
    // loop over the slices 
    for (size_t i = 0; i < xMinimum.size() ; ++i) { //-1
      double xMin = xMinimum[i] ;
      int binMin = hist2D ->GetXaxis()->FindBin(xMin)+1;

      double xMax = xMaximum[i];
      int binMax = hist2D ->GetXaxis()->FindBin(xMax)-1;

      // Define the name and title for the TH1D histogram
      // Format the histogram name to display 4 decimal places
      std::ostringstream stream;
      stream << std::fixed << std::setprecision(4) << "_"<<xMin << "_to_" << xMax;
      std::string histName = yAxisName+"__"+xAxisName + stream.str() + "_" +type;
      TH1D *projY = hist2D->ProjectionY(histName.c_str(), binMin, binMax);
      histVector.push_back(projY); //

      cout<<"xMin = "<<xMin<<" , xMax = "<<xMax<<endl;
      cout<< "FindBin(xMin) = "<<hist2D ->GetXaxis()->FindBin(xMin)<<", using binMin= "<<binMin<< " to be exclusive"<<endl;
       cout<< "FindBin(xMax) = "<<hist2D ->GetXaxis()->FindBin(xMax)<<", using binMin= "<<binMax<< " to be exclusive"<<endl;
      

    }// end loop over slices
}// end SliceAndProjectHistogram_AroundMean



void SliceAndProjectHistogram_AroundMean(TH2D* hist2D, const std::vector<double>& xSlices, std::vector<TH1D*>& histVector, double xMean, std::string xAxisName, std::string yAxisName, std::string type) {
    // Clear the vector to ensure it's empty before filling
    histVector.clear();
  
    // loop over the slices 
    for (size_t i = 0; i < xSlices.size() ; ++i) { //-1
      double xMin = xMean - xSlices[i];
      int binMin = hist2D ->GetXaxis()->FindBin(xMin);

      double xMax = xMean + xSlices[i];
      int binMax = hist2D ->GetXaxis()->FindBin(xMax);

      // Define the name and title for the TH1D histogram
      // Format the histogram name to display only two decimal places
      std::ostringstream stream;
      stream << std::fixed << std::setprecision(4) << "_"<<xMin << "_to_" << xMax;
      std::string histName = yAxisName+"__"+xAxisName + stream.str() + "_" +type;
      TH1D *projY = hist2D->ProjectionY(histName.c_str(), binMin, binMax);
      histVector.push_back(projY); // 

    }// end loop over slices
}// end SliceAndProjectHistogram_AroundMean



void SliceAndProjectHistogram_xMaxFixed(TH2D* hist2D, const std::vector<double>& xSlices, std::vector<TH1D*>& histVector, double xMax, std::string xAxisName, std::string yAxisName, std::string type) {
    // Clear the vector to ensure it's empty before filling
    histVector.clear();
    // Find the bin number of the x endpoint 
    int binMax = hist2D->GetXaxis()->FindBin(xMax);

    // loop over the slices 
    for (size_t i = 0; i < xSlices.size() ; ++i) { //-1
      double xMin = xSlices[i];
      int binMin = hist2D ->GetXaxis()->FindBin(xMin);

      // Define the name and title for the TH1D histogram
      // Format the histogram name to display only two decimal places
      std::ostringstream stream;
      stream << std::fixed << std::setprecision(4) << "_"<<xMin << "_to_" << xMax;
      std::string histName = yAxisName+"__"+xAxisName + stream.str() + "_" +type;
      TH1D *projY = hist2D->ProjectionY(histName.c_str(), binMin, binMax);
      histVector.push_back(projY); // 

    }// end loop over slices
}// end SliceAndProjectHistogram_xMaxFixed


// Fit that is a combination of the scaled proton mc, scaled neutron mc, and 2nd order polynomial. 
//// hist_p and hist_n are globals that need to be set to the proper histograms right before you call the fit function. 
Double_t mc_p_n_poly2_slice_fit(Double_t *x, Double_t *par) {
    
    Double_t val = 0.0;

    // Get x value
    Double_t xx = x[0];

    // Retrieve parameters
    Double_t scale_p = par[0];
    Double_t scale_n = par[1];
    Double_t shift_p = par[2];
    Double_t shift_n = par[3];

    Double_t polyCoefficients[3]; 
    
    // Fill polynomial coefficients 
    for (Int_t i = 0; i <=2; ++i) {
        polyCoefficients[i] = par[i + 4];
    }

    //// hist_p and hist_n are globals that need to be set to the proper histograms right before you call the fit function. 

    // Calculate value using combination of histograms and polynomial background
    val = scale_p *hist_p ->Interpolate(xx-shift_p) + scale_n *hist_n ->Interpolate(xx-shift_n);

    // Add polynomial background
    for (Int_t i = 0; i <= 2; ++i) {
        val += polyCoefficients[i] * TMath::Power(xx, i);
    }

    return val;
}




// Fit that is a combination of the scaled proton mc, scaled neutron mc, and 4th order polynomial. 
//// hist_p and hist_n are globals that need to be set to the proper histograms right before you call the fit function. 
Double_t mc_p_n_poly4_slice_fit(Double_t *x, Double_t *par) {
    
    Double_t val = 0.0;

    // Get x value
    Double_t xx = x[0];

    // Retrieve parameters
    Double_t scale_p = par[0];
    Double_t scale_n = par[1];
    Double_t shift_p = par[2];
    Double_t shift_n = par[3];

    Double_t polyCoefficients[5]; 
    
    // Fill polynomial coefficients 
    for (Int_t i = 0; i <=4; ++i) {
        polyCoefficients[i] = par[i + 4];
    }

    //// hist_p and hist_n are globals that need to be set to the proper histograms right before you call the fit function. 

    // Calculate value using combination of histograms and polynomial background
    val = scale_p *hist_p ->Interpolate(xx-shift_p) + scale_n *hist_n ->Interpolate(xx-shift_n);

    // Add polynomial background
    for (Int_t i = 0; i <= 4; ++i) {
        val += polyCoefficients[i] * TMath::Power(xx, i);
    }

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

Double_t poly2(Double_t *x, Double_t *par)
{
  Double_t fit = par[0] + par[1] * x[0] + par[2] * pow(x[0],2) ;
  return fit;
}

Double_t poly4(Double_t *x, Double_t *par)
{
  Double_t fit = par[0] + par[1] * x[0] + par[2] * pow(x[0],2) + par[3] * pow(x[0],3) + par[4] * pow(x[0],4) ;
  return fit;
}


TH1D* GetResidualHistogram(TH1D* hist, TF1* fit) {
 // Create a new histogram with the same binning as the original
  TH1D *h_residual = (TH1D*)(hist->Clone(Form("%s_residual",hist->GetName())));
  TF1 *fit_clone = (TF1*)(fit->Clone(Form("%s_clone",fit->GetName())));

  // Loop over bins and calculate residuals
  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
    double binCenter = hist->GetBinCenter(i);
    double binContent = hist->GetBinContent(i);
    double fitValue = fit_clone->Eval(binCenter);
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
    double ex[numBins]; // No x errors 
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


    graph->SetTitle("");
    return graph;
}

void customizeGraph(TGraphErrors *graph, int markerStyle, int markerColor, double markerSize, 
                    const std::string& graphTitle ="", const std::string& xAxisLabel="", const std::string& yAxisLabel="",
		    double TitleOffsetX = 1.4, double TitleOffsetY = 2, 
		    double LabelOffsetX = 0.01, double LabelOffsetY = 0.01) {
  // Set marker style, color, and size
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerColor(markerColor);
  graph->SetMarkerSize(markerSize);

  // Set graph title and axis labels
  graph->SetTitle(graphTitle.c_str());
  graph->GetXaxis()->SetTitle(xAxisLabel.c_str());
  graph->GetYaxis()->SetTitle(yAxisLabel.c_str());

// Adjust axis title offsets to provide more space
  graph->GetXaxis()->SetTitleOffset(TitleOffsetX); // Adjust as needed
  graph->GetYaxis()->SetTitleOffset(TitleOffsetY); // Adjust as needed


// Adjust axis label offsets to provide more space
  graph->GetXaxis()->SetLabelOffset(LabelOffsetX); // Adjust as needed
  graph->GetYaxis()->SetLabelOffset(LabelOffsetY); // Adjust as needed
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




//// Get the title from the histogram and display it on a canvas. 
/// This expects the title to be in the form:
///  y_axis:x_axis {cut1&&cut2&&cut3....}
void printParsedTitle(const std::string& title, TString outputlocation) {

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
    
    // Optionally, save the canvas as an image
   cuts_canvas->SaveAs(Form("%s/global_cuts.pdf", outputlocation.Data()));
}


double CalculateMean(const std::vector<double>& data){
  double sum = 0;
  for (double value : data) {
    sum += value;
    }
  return sum/data.size();
}

double CalculateStDev(const std::vector<double>& data){
  double mean = CalculateMean(data);
  double numerator_sum = 0;
  for (double value : data){
    numerator_sum += pow( (value-mean),2);
  }
  double variance = numerator_sum/data.size();
  return sqrt(variance);
}


double CalculateWeightedMean(const std::vector<double>& data, const std::vector<double>& uncert){
  double numerator_sum=0;
  double sum_of_weights=0;
  double value =0;
  double weight = 0;
  if (data.size()!=uncert.size()){
    cout<<"data and uncert vectors not the same size in Calculate Weighted Mean"<<endl;
  }
  for (int i = 0;i<data.size();i++) {
    value = data[i];
    weight = 1/( pow(uncert[i],2));
    numerator_sum += value*weight;
    sum_of_weights += weight;
    }
  return numerator_sum/sum_of_weights;
}

double CalculateWeightedStDev(const std::vector<double>& data, const std::vector<double>& uncert){
  double numerator_sum=0;
  double sum_of_weights=0;
  double value =0;
  double weight = 0;
  double mean_w = CalculateWeightedMean(data, uncert);
  if (data.size()!=uncert.size()){
    cout<<"data and uncert vectors not the same size in Calculate Weighted Mean"<<endl;
  }
  for (int i = 0;i<data.size();i++) {
    value = data[i];
    weight = 1/( pow(uncert[i],2));
    numerator_sum += weight*pow((value-mean_w),2);
    sum_of_weights += weight;
    }
  double variance = numerator_sum/sum_of_weights;
  return sqrt(variance);
}





void parseStringToVectors(const std::string& input, std::vector<double>& xMin, std::vector<double>& xMax) {
    std::stringstream ss(input);
    std::string token;

    // Clear the input vectors in case they contain previous data
    xMin.clear();
    xMax.clear();

    while (std::getline(ss, token, ')')) {
        // Find the opening parenthesis
        std::size_t openParen = token.find('(');
        if (openParen == std::string::npos) continue;

        // Remove parentheses from the token
        std::string cleanToken = token.substr(openParen + 1);
        cleanToken.erase(remove(cleanToken.begin(), cleanToken.end(), ')'), cleanToken.end());

        // Split the token into xMin and xMax values
        std::stringstream pairStream(cleanToken);
        std::string xMinStr, xMaxStr;

        std::getline(pairStream, xMinStr, ',');
        std::getline(pairStream, xMaxStr, ',');

        try {
            // Convert strings to doubles and add to respective vectors
            xMin.push_back(std::stod(xMinStr));
            xMax.push_back(std::stod(xMaxStr));
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid argument: " << e.what() << std::endl;
        } catch (const std::out_of_range& e) {
            std::cerr << "Out of range: " << e.what() << std::endl;
        }
    }
}
