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


#include "/work/halla//sbs/msatnik/GMn/classes/FitHistogram.h"
#include "/work/halla//sbs/msatnik/GMn/classes/FitHistogram.cpp"

#include "/work/halla//sbs/msatnik/GMn/classes/Utility.h"
#include "/work/halla//sbs/msatnik/GMn/classes/Utility.cpp"

#include "/work/halla//sbs/msatnik/GMn/classes/FileNames.h"
#include "/work/halla//sbs/msatnik/GMn/classes/FileNames.cpp"

//// Vectors to store the sliced histograms in. Declaring them outside so all the functions can see them. 
std::vector<TH1D*> hist_vector_data; 
std::vector<TH1D*> hist_vector_p; 
std::vector<TH1D*> hist_vector_n; 

std::vector<TH1D*> hist_result_vector_p; 
std::vector<TH1D*> hist_result_vector_n; 

std::vector<TH1D*> hist_residual_vector; 



//// default values to be overwritten 
TString DataFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_tight_2Dhistos.root";
TString ProtonFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_tight_2Dhistos.root";
TString NeutronFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_tight_2Dhistos.root";
TString HistogramName = "hcal_dx_1d_allcuts";
std::string xMin_xMax_string = "(-4.1,-2.1),(-3.1,-1.1),(-2,0.1),(-1,1),(0,2),(1.1,3),(2.1,4)";
int PolyOrder = 2;
double xmin_start = -0.74 - 0.4;
double xmin_step = -0.10;
double xmax_start= 0 + 0.4;
double xmax_step = 0.10;
int nSteps = 16;


// SBS4 30p 
// TString DataFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_2Dhistos_sept26.root";
// TString ProtonFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_2Dhistos_sept26.root";
// TString NeutronFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_2Dhistos_sept26.root";
std::string xMin_xMax_string_sbs4_30p = "(-2.75,2.0),(-2.65,1.9),(-2.55,1.8),(-2.45,1.7),(-2.35,1.6),(-2.25,1.5),(-2.15,1.4),(-2.05,1.3),(-1.95,1.2),(-1.85,1.1),(-1.75,1.0),(-1.65,0.9),(-1.55,0.8),(-1.45,0.7),(-1.35,0.6),(-1.25,0.5)";
double xmin_start_sbs4_30p = -0.74 - 0.4;
double xmin_step_sbs4_30p = -0.10;
double xmax_start_sbs4_30p= 0 + 0.4;
double xmax_step_sbs4_30p = 0.10;
int nSteps_sbs4_30p = 16;

int PolyOrder_sbs4_30p = 2;

// SBS4 50p 
// TString DataFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_2Dhistos.root";
// TString ProtonFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deep_2Dhistos.root";
// TString NeutronFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deen_2Dhistos.root";
std::string xMin_xMax_string_sbs4_50p = "(-3.2,2.0),(-3.1,1.9),(-3.0,1.8),(-2.9,1.7),(-2.8,1.6),(-2.7,1.5),(-2.6,1.4),(-2.5,1.3),(-2.4,1.2),(-2.3,1.1),(-2.2,1.0),(-2.1,0.9)";
int PolyOrder_sbs4_50p = 2;
double xmin_start_sbs4_50p = -1.22 - 0.4;
double xmin_step_sbs4_50p = -0.10;
double xmax_start_sbs4_50p= 0 + 0.4;
double xmax_step_sbs4_50p = 0.10;
int nSteps_sbs4_50p = 16;


// SBS8 70p
// TString DataFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_2Dhistos_sept11.root";
// TString ProtonFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deep_2Dhistos.root";
// TString NeutronFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deen_2Dhistos.root";
std::string xMin_xMax_string_sbs8_70p = "(-3.2,2.0),(-3.1,1.9),(-3.0,1.8),(-2.9,1.7),(-2.8,1.6),(-2.7,1.5),(-2.6,1.4),(-2.5,1.3),(-2.4,1.2),(-2.3,1.1),(-2.2,1.0),(-2.1,0.9),(-2,0.8),(-1.9,0.7),(-1.8,0.6),(-1.7,0.5)";
int PolyOrder_sbs8_70p = 4;
double xmin_start_sbs8_70p = -0.87 - 0.4;
double xmin_step_sbs8_70p = -0.10;
double xmax_start_sbs8_70p= 0 + 0.4;
double xmax_step_sbs8_70p = 0.10;
int nSteps_sbs8_70p = 16;

// SBS8 50p
double xmin_start_sbs8_50p = -0.61 - 0.4;
double xmin_step_sbs8_50p = -0.10;
double xmax_start_sbs8_50p= 0 + 0.4;
double xmax_step_sbs8_50p = 0.10;
int nSteps_sbs8_50p = 16;

// SBS8 100p
double xmin_start_sbs8_100p = -1.23 - 0.4;
double xmin_step_sbs8_100p = -0.10;
double xmax_start_sbs8_100p= 0 + 0.4;
double xmax_step_sbs8_100p = 0.10;
int nSteps_sbs8_100p = 16;

// SBS9 70p
double xmin_start_sbs9_70p = -0.91 - 0.4;
double xmin_step_sbs9_70p = -0.10;
double xmax_start_sbs9_70p= 0 + 0.4;
double xmax_step_sbs9_70p = 0.10;
int nSteps_sbs9_70p = 16;

// Functions 



void fit_range_stability(TString KineString="sbs4_30p", TString IterationString="both", int PolyOrder_input = 0){ // main
 
  gStyle->SetNumberContours(255); 
  // gStyle->SetOptStat(0110);
  gStyle->SetOptStat(0);
  
  Utility utilityHandler; // class that gives us access to various functions to use between programs.
  
  FileNames fileNamesHandler;
 

  //// Histogram to store the Rsf distriubtion
  TH1D *h_Rsf = new TH1D("h_Rsf","h_Rsf",300, 0.8,1.2);

  /// Histogram for the pull distribution
  TH1D *h_pull = new TH1D("h_pull","h_pull",10,-1,1);

  
  if ( !(IterationString == "both" || IterationString == "left"|| IterationString=="right") ){
    std::cout<<"Error with Itteration String "<<std::endl;
    return;
  }

										       
   
  if (KineString == "sbs4_30p"){
    DataFileString =fileNamesHandler.DataFileString_sbs4_30p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs4_30p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs4_30p;
    if (IterationString == "both"){
      xMin_xMax_string = utilityHandler.incrementRangeStringStepsizeNotNested(xmin_start_sbs4_30p, xmax_start_sbs4_30p,xmin_step_sbs4_30p, xmax_step_sbs4_30p, nSteps_sbs4_30p, 2);
    }else if (IterationString == "left"){
      xMin_xMax_string = utilityHandler.incrementRangeStringStepsizeNotNested(xmin_start_sbs4_30p, xmax_start_sbs4_30p,xmin_step_sbs4_30p, 0, nSteps_sbs4_30p, 2);
    } else if (IterationString == "right"){
      xMin_xMax_string = utilityHandler.incrementRangeStringStepsizeNotNested(xmin_start_sbs4_30p, xmax_start_sbs4_30p,0, xmax_step_sbs4_30p, nSteps_sbs4_30p, 2);
    }
  }else if (KineString == "sbs4_50p"){
    DataFileString = fileNamesHandler.DataFileString_sbs4_50p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs4_50p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs4_50p;
    if (IterationString == "both"){
      xMin_xMax_string = utilityHandler.incrementRangeStringStepsizeNotNested(xmin_start_sbs4_50p, xmax_start_sbs4_50p,xmin_step_sbs4_50p, xmax_step_sbs4_50p, nSteps_sbs4_50p, 2);
    }else if (IterationString == "left"){
      xMin_xMax_string = utilityHandler.incrementRangeStringStepsizeNotNested(xmin_start_sbs4_50p, xmax_start_sbs4_50p,xmin_step_sbs4_50p, 0, nSteps_sbs4_50p, 2);
    } else if (IterationString == "right"){
      xMin_xMax_string = utilityHandler.incrementRangeStringStepsizeNotNested(xmin_start_sbs4_50p, xmax_start_sbs4_50p,0, xmax_step_sbs4_50p, nSteps_sbs4_50p, 2);
    }
  }else if (KineString == "sbs8_70p"){
    DataFileString = fileNamesHandler.DataFileString_sbs8_70p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs8_70p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs8_70p;
    if (IterationString == "both"){
      xMin_xMax_string = utilityHandler.incrementRangeStringStepsizeNotNested(xmin_start_sbs8_70p, xmax_start_sbs8_70p,xmin_step_sbs8_70p, xmax_step_sbs8_70p, nSteps_sbs8_70p, 2);
    }else if (IterationString == "left"){
      xMin_xMax_string = utilityHandler.incrementRangeStringStepsizeNotNested(xmin_start_sbs8_70p, xmax_start_sbs8_70p,xmin_step_sbs8_70p, 0, nSteps_sbs8_70p, 2);
    } else if (IterationString == "right"){
      xMin_xMax_string = utilityHandler.incrementRangeStringStepsizeNotNested(xmin_start_sbs8_70p, xmax_start_sbs8_70p,0, xmax_step_sbs8_70p, nSteps_sbs8_70p, 2);
    }
  }else if (KineString == "sbs8_50p"){
    DataFileString = fileNamesHandler.DataFileString_sbs8_50p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs8_50p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs8_50p;
    if (IterationString == "both"){
      xMin_xMax_string = utilityHandler.incrementRangeStringStepsizeNotNested(xmin_start_sbs8_50p, xmax_start_sbs8_50p,xmin_step_sbs8_50p, xmax_step_sbs8_50p, nSteps_sbs8_50p, 2);
    }else if (IterationString == "left"){
      xMin_xMax_string = utilityHandler.incrementRangeStringStepsizeNotNested(xmin_start_sbs8_50p, xmax_start_sbs8_50p,xmin_step_sbs8_50p, 0, nSteps_sbs8_50p, 2);
    } else if (IterationString == "right"){
      xMin_xMax_string = utilityHandler.incrementRangeStringStepsizeNotNested(xmin_start_sbs8_50p, xmax_start_sbs8_50p,0, xmax_step_sbs8_50p, nSteps_sbs8_50p, 2);
    }
  }else if (KineString == "sbs8_100p"){
    DataFileString = fileNamesHandler.DataFileString_sbs8_100p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs8_100p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs8_100p;
    if (IterationString == "both"){
      xMin_xMax_string = utilityHandler.incrementRangeStringStepsizeNotNested(xmin_start_sbs8_100p, xmax_start_sbs8_100p,xmin_step_sbs8_100p, xmax_step_sbs8_100p, nSteps_sbs8_100p, 2);
    }else if (IterationString == "left"){
      xMin_xMax_string = utilityHandler.incrementRangeStringStepsizeNotNested(xmin_start_sbs8_100p, xmax_start_sbs8_100p,xmin_step_sbs8_100p, 0, nSteps_sbs8_100p, 2);
    } else if (IterationString == "right"){
      xMin_xMax_string = utilityHandler.incrementRangeStringStepsizeNotNested(xmin_start_sbs8_100p, xmax_start_sbs8_100p,0, xmax_step_sbs8_100p, nSteps_sbs8_100p, 2);
    }
  } else if (KineString == "sbs9_70p"){
    DataFileString = fileNamesHandler.DataFileString_sbs9_70p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs9_70p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs9_70p;
    if (IterationString == "both"){
      xMin_xMax_string = utilityHandler.incrementRangeStringStepsizeNotNested(xmin_start_sbs9_70p, xmax_start_sbs9_70p,xmin_step_sbs9_70p, xmax_step_sbs9_70p, nSteps_sbs9_70p, 2);
    }else if (IterationString == "left"){
      xMin_xMax_string = utilityHandler.incrementRangeStringStepsizeNotNested(xmin_start_sbs9_70p, xmax_start_sbs9_70p,xmin_step_sbs9_70p, 0, nSteps_sbs9_70p, 2);
    } else if (IterationString == "right"){
      xMin_xMax_string = utilityHandler.incrementRangeStringStepsizeNotNested(xmin_start_sbs9_70p, xmax_start_sbs9_70p,0, xmax_step_sbs9_70p, nSteps_sbs9_70p, 2);
    }
  }else {
    std::cout<<"Error with kinematic setting"<<std::endl;
    return;
  }


  if((PolyOrder_input < 0) && PolyOrder_input > 6){
    cout<<"error with poly order input. Setting to 2."<<endl;
    PolyOrder = 2;
  }else PolyOrder=PolyOrder_input; 
 
  cout<<"What we're running: "<<KineString<< ", "<< IterationString<<endl;
  cout<<DataFileString<<endl;
  cout<<ProtonFileString<<endl;
  cout<<NeutronFileString<<endl;
  cout<<HistogramName<<endl;
  cout<<"poly order: "<<PolyOrder <<endl;
  cout<<xMin_xMax_string<<endl;
 
      
      
  
  TFile *f1 = TFile::Open(DataFileString); // data
  TFile *f2 = TFile::Open(ProtonFileString); //proton
  TFile *f3 = TFile::Open(NeutronFileString); //neutron


  // set location and name for output file 
  TString outputfilelocation="../output/stability/"+KineString+"/fitRange";
  TString outputfilename = outputfilelocation +"/"+ KineString+ "_.root";

  // Declare outfile
  TFile *fout = new TFile(outputfilename,"RECREATE");
  cout<<"writing to file: "<< outputfilename <<endl;


  // Load Histograms
  TH1D *hist_data_orig = (TH1D*)f1->Get(HistogramName.Data());
  TH1D *hist_proton_orig_1d  = (TH1D*)f2->Get(HistogramName.Data());
  TH1D *hist_neutron_orig_1d  =  (TH1D*)f3->Get(HistogramName.Data());

  // Make clones of the histograms
  TH1D *hist_1D_data = (TH1D*)hist_data_orig->Clone("hist_1D_data");
  TH1D *hist_p = (TH1D*)hist_proton_orig_1d->Clone("hist_p"); 
  TH1D *hist_n = (TH1D*)hist_neutron_orig_1d->Clone("hist_n");

 

  // Set up the slices for the Y projections. 
  // std::vector<double> xSlices = {-8,-5,-3,-2.5,-2,-1.5,-1,-0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3,5,8};
  
  // std::vector<double> xMin_range =  {-3.2,-2.7,-2.2,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9};
  // std::vector<double> xMax_range = { 2.5, 2.0, 1.5, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2};

  // std::vector<double> xMin_range =  {-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6};
  // std::vector<double> xMax_range = { 2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9};


  std::vector<double> xMin_range;
 std:vector<double> xMax_range;
  utilityHandler.parseStringToVectors(xMin_xMax_string, xMin_range, xMax_range);
  
  if (xMin_range.size() != xMax_range.size())
    {
      cout<<"error with the fit ranges"<<endl;
      // exit;
    }
  int nRanges = xMin_range.size();
  
    
  std::vector<FitHistogram> fitHandler_vector; // vector of objects fo the FitHistogram class that will help us do data-mc comparion. .
  

  //// Vectors to store the sliced histograms in. 
  std::vector<TH1D*> hist_vector_data; 
 

  for (int i = 0; i<nRanges ;  i++)
    {
      TH1D* cloneHist = (TH1D*)hist_1D_data->Clone(Form("clone_%d", i));
      hist_vector_data.push_back(cloneHist);
    }
  


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

 
      
  TH1D *hist_data;
  TH1D *hist_data_dont_draw;
  /// Fit each slice 
  for (int rangeID = 0 ; rangeID < nRanges ; rangeID++)
    {
      /// set up the global histograms that the fit function is going to use
      // this uses clones of the histograms
      hist_data  =(TH1D*)hist_vector_data[rangeID]->Clone("hist_data");
      hist_data_dont_draw = (TH1D*)hist_vector_data[rangeID]->Clone("hist_data_dont_draw");
    

      if (!hist_data_dont_draw || !hist_p || !hist_n) {
	std::cerr << "Error: One of the histograms is null!" << std::endl;
	return;
      }
    

      FitHistogram fitHandler(hist_p,hist_n, xMin_range[rangeID],xMax_range[rangeID]);
      fitHandler.setPolyOrder(PolyOrder);// setting the order for the polynomial background
      fitHandler.fitDataPoly(hist_data_dont_draw);

      TF1 *fit_result = (TF1*)fitHandler.fitFunc->Clone("fit_result");// bug in here that corrupts it and crashes if you try and draw
	

      // get fit results 
      double scale_p = fitHandler.scale_p;
      double scale_n=  fitHandler.scale_n;
      double scale_p_err = fitHandler.scale_p_err;
      double scale_n_err = fitHandler.scale_n_err;
      double Rsf= fitHandler.R;
      double Rsf_err= fitHandler.R_err;
      double shift_p = fitHandler.shift_p;
      double shift_p_err = fitHandler.shift_p_err;
      double shift_n = fitHandler.shift_n;
      double shift_n_err = fitHandler.shift_n_err;
      double ChiSq = fitHandler.ChiSq;
      double ndf = fitHandler.NDF;

      std::vector<double> poly_result;
      std::vector<double> poly_result_err;
      for (int i =0 ; i <= PolyOrder; i++)
	{
	  poly_result.push_back( fitHandler.poly_result[i]);
	  poly_result_err.push_back(fitHandler.poly_result[i] );
	}


      cout<<"rangeID: "<<rangeID<<endl;
      cout<<"scale_p = "<< scale_p <<endl;
      cout<<"Rsf = "<< Rsf <<" +/- "<< Rsf_err<<endl;
      cout<<"shift p: "<<shift_p<<" +/- "<<shift_p_err<<", shift_n: "<<shift_n<<"+/-"<<shift_n_err<<endl;
      cout<<endl;

      h_Rsf->Fill(Rsf);
      
      //// save results into the vectors 
      fit_vector.push_back(fit_result);
      scale_p_vector.push_back(scale_p);
      scale_n_vector.push_back(scale_n);
      scale_p_err_vector.push_back(scale_p_err);
      scale_n_err_vector.push_back(scale_n_err);
      shift_p_vector.push_back(shift_p);
      shift_n_vector.push_back(shift_n);
      shift_p_vector.push_back(shift_p);
      shift_n_vector.push_back(shift_n);
      ChiSq_vector.push_back(ChiSq);
      ndf_vector.push_back(ndf);
      Rsf_vector.push_back(Rsf);
      Rsf_err_vector.push_back(Rsf_err);
      // poly_result_vector_of_arrays.push_back(poly_result);
      poly_result_vector_of_vectors.push_back(poly_result);
      poly_result_err_vector_of_vectors.push_back(poly_result_err);
	  
      //delete Fit; //
    }// end loop over slices



     // Print the contents of the vector of vectors
  for (const auto& vec : poly_result_vector_of_vectors) {
    for (double val : vec) {
      // std::cout << val << " ";
    }
    std::cout << std::endl;
  }

  
  
  cout<<"---------------------------------------------------------------------------------------------"<<endl;
  double Rsf_mean = utilityHandler.CalculateMean(Rsf_vector);
  double Rsf_stdev = utilityHandler.CalculateStDev(Rsf_vector);
  double Rsf_mean_w = utilityHandler.CalculateWeightedMean(Rsf_vector,Rsf_err_vector);
  double Rsf_stdev_w = utilityHandler.CalculateWeightedStDev(Rsf_vector,Rsf_err_vector);
  cout<<"Rsf mean = " << Rsf_mean<<endl;
  cout<<"Rsf StDev = " << Rsf_stdev<<endl;
  cout<<"weighted Rsf mean = "<<Rsf_mean_w<<endl;
  cout<<"weighted Rsf StDev = "<<Rsf_stdev_w<<endl;
  cout<<"---------------------------------------------------------------------------------------------"<<endl;
  cout<<endl;

  // calculate the Pull for each point and put in the histogram 
  for (int i = 0; i < Rsf_vector.size();i++)
    {
      double pull = (Rsf_mean - Rsf_vector[i]) / Rsf_err_vector[i];
      h_pull ->Fill(pull);
    }
  

  TH1D *hist_temp_p;
  TH1D *hist_temp_n;

  // loop over the slices to create shifted and scaled versions of the mc histograms to plot
  for (int rangeID = 0 ; rangeID < nRanges ; rangeID++)
    {	
      // utility function that scales and shifts the provided histogram and returns a new histogram 
      hist_temp_p = utilityHandler.ScaleAndShiftHistogram(hist_p, scale_p_vector[rangeID], shift_p_vector[rangeID]);
      hist_temp_n = utilityHandler.ScaleAndShiftHistogram(hist_n, scale_n_vector[rangeID], shift_n_vector[rangeID]);
	 
      // save histograms to the global vectors 
      hist_result_vector_p.push_back(hist_temp_p);
      hist_result_vector_n.push_back(hist_temp_n);
    }// end loop over slices

  /// Make a canvas called c1 to overwrite the default canvas and avoid crashes
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  hist_result_vector_p[0] ->Draw("E");
  

  std::vector<TF1*> poly_fit_result;
  // loop over slices and make fits to plot the polynomial background result
  for (int rangeID = 0 ; rangeID < nRanges ; rangeID++)
    {
      // get the fit results from the vector of vectors
      std::vector<double> poly_params_vector;
      poly_params_vector = poly_result_vector_of_vectors[rangeID];

      // convert the vector into an array that TF1 can use
      int size = poly_params_vector.size();
      double poly_params_array[size]; 
      for (int i  = 0 ; i < size; i++){
  	poly_params_array[i] = poly_params_vector[i];
      }

      TF1 *fit = new TF1(Form("poly_%i_",rangeID),Form("pol%i",PolyOrder),xMin_range[rangeID],xMax_range[rangeID] );
      fit->SetParameters(poly_params_array); 
      fit->SetNpx(500);
      poly_fit_result.push_back(fit);

      if (utilityHandler.doesFunctionGoBelowZero(fit,xMin_range[rangeID],xMax_range[rangeID]))
	{
	  // background function went below zero, which isn't physical
	  std::cout<<"On slice ID: "<<rangeID<<std::endl;
	  std::cout<<"---------------------------"<<std::endl;
	}
      
    }// end loop over slices


     //// Plot Rsf
     //// Make arrays that TGraphErrors can use 
  double x[nRanges];
  double y[nRanges];
  double x_err[nRanges];
  double y_err[nRanges];
  for (int rangeID = 0 ; rangeID < nRanges ; rangeID++)
    {
      x[rangeID] = xMax_range[rangeID];
      y[rangeID] = Rsf_vector[rangeID];
      x_err[rangeID] = 0;
      y_err[rangeID] = Rsf_err_vector[rangeID];
    }    
  TGraphErrors *Rsf_graph = new TGraphErrors(nRanges, x, y, x_err, y_err);
  utilityHandler.customizeGraphMore(Rsf_graph, 33, kBlue, 3,"","xMax","Rsf");

  // Fit a straight line to the graph
  TF1* fit_Rsf_graph = new TF1("fit_Rsf_graph", "[0]",Rsf_graph ->GetX()[0], Rsf_graph->GetX()[Rsf_graph->GetN()-1]);
  Rsf_graph->Fit(fit_Rsf_graph, "Q R");
  
  //// canvas
  TCanvas *graphcanvas = new TCanvas("graphcanvas","graphcanvas",800,600);  
  utilityHandler.adjustCanvas(graphcanvas);
  Rsf_graph->Draw("AP");
  graphcanvas->Update();
  // Get fit parameters
  double constant = fit_Rsf_graph->GetParameter(0);       // The constant value
  double constantError = fit_Rsf_graph->GetParError(0);    // The error on the constant
  double chi2 = fit_Rsf_graph->GetChisquare();             // The chi-squared value
  int ndf = fit_Rsf_graph->GetNDF();                       // The number of degrees of freedom
  double chi2_ndf = chi2 / ndf;                          // chi2/ndf

  // Draw the graph with the fit result again
  Rsf_graph->Draw("AP");

  // Use TLatex to add the fit result and chi²/ndf to the canvas
  TLatex latex;
  latex.SetNDC();  // Use normalized coordinates (0 to 1)
  latex.SetTextSize(0.04);
  latex.DrawLatex(0.2, 0.85, Form("Fit: y = %.5f #pm %.5f", constant, constantError));
  latex.DrawLatex(0.2, 0.80, Form("#chi^{2}/ndf = %.2f/%i", chi2,ndf));
  latex.DrawLatex(0.2, 0.2, Form("Rsf Mean = %.5f #pm StDev = %.5f", Rsf_mean, Rsf_stdev));
  
  graphcanvas->Update();
  graphcanvas->SaveAs(Form("%s/Rsf.pdf",outputfilelocation.Data() ));


  //// Plot Rsf
  //// Make arrays that TGraphErrors can use 
  double x_min[nRanges];
  double y_min[nRanges];
  double x_min_err[nRanges];
  double y_min_err[nRanges];
  for (int rangeID = 0 ; rangeID < nRanges ; rangeID++)
    {
      x_min[rangeID] = xMin_range[rangeID];
      y_min[rangeID] = Rsf_vector[rangeID];
      x_min_err[rangeID] = 0;
      y_min_err[rangeID] = Rsf_err_vector[rangeID];
    }    
  TGraphErrors *Rsf_graph_xmin = new TGraphErrors(nRanges, x_min, y_min, x_min_err, y_min_err);
  utilityHandler.customizeGraphMore(Rsf_graph_xmin, 33, kBlue, 3,"","xMin","Rsf");

  // Fit a straight line to the graph
  TF1* fit_Rsf_graph_xmin = new TF1("fit_Rsf_graph_xmin", "[0]",Rsf_graph_xmin ->GetX()[0], Rsf_graph_xmin->GetX()[Rsf_graph_xmin->GetN()-1]);
  Rsf_graph_xmin->Fit(fit_Rsf_graph_xmin, "Q R");
  
  //// canvas
  TCanvas *xMin_canvas = new TCanvas("xMin_canvas","xMin_canvas",800,600);  
  utilityHandler.adjustCanvas(xMin_canvas);
  Rsf_graph_xmin->Draw("AP");
  xMin_canvas->Update();
  // Get fit parameters
  double constant_min = fit_Rsf_graph_xmin->GetParameter(0);       // The constant value
  double constantError_min = fit_Rsf_graph_xmin->GetParError(0);    // The error on the constant
  double chi2_min = fit_Rsf_graph_xmin->GetChisquare();             // The chi-squared value
  int ndf_min = fit_Rsf_graph_xmin->GetNDF();                       // The number of degrees of freedom
  double chi2_ndf_min = chi2_min / ndf_min;                          // chi2/ndf

  // Draw the graph with the fit result again
  Rsf_graph_xmin->Draw("AP");

  // Use TLatex to add the fit result and chi²/ndf to the canvas
  TLatex latex_min;
  latex_min.SetNDC();  // Use normalized coordinates (0 to 1)
  latex_min.SetTextSize(0.04);
  latex_min.DrawLatex(0.2, 0.85, Form("Fit: y = %.5f #pm %.5f", constant_min, constantError_min));
  latex_min.DrawLatex(0.2, 0.80, Form("#chi^{2}/ndf = %.2f/%i", chi2_min,ndf_min));
  latex.DrawLatex(0.2, 0.2, Form("Rsf Mean = %.5f #pm StDev = %.5f", Rsf_mean, Rsf_stdev));
  
  xMin_canvas->Update();
  xMin_canvas->SaveAs(Form("%s/Rsf_xMin.pdf",outputfilelocation.Data() ));

  
  //// Plot Chi2/ndf
  //// Make arrays that TGraphErrors can use
  double x_ch[nRanges];
  double y_ch[nRanges];
  double x_ch_err[nRanges];
  double y_ch_err[nRanges];
  for (int rangeID = 0 ; rangeID < nRanges ; rangeID++)
    {
      x_ch[rangeID] = xMax_range[rangeID];
      y_ch[rangeID] = ChiSq_vector[rangeID] /ndf_vector[rangeID];
      x_ch_err[rangeID] = 0;
      y_ch_err[rangeID] = 0;
    }    
  TGraphErrors *Chi2_ndf_graph = new TGraphErrors(nRanges, x_ch, y_ch, x_ch_err, y_ch_err);
  utilityHandler.customizeGraphMore(Chi2_ndf_graph, 33, kBlue, 3,"","xMax","Chi^{2}/ndf");

  //// canvas
  TCanvas *chi2_ndf_canvas = new TCanvas("chi2_ndf_canvas","chi2_ndf_canvas",800,600);  
  utilityHandler.adjustCanvas(chi2_ndf_canvas);
  Chi2_ndf_graph ->Draw("AP");
  chi2_ndf_canvas->Update();
  chi2_ndf_canvas->SaveAs(Form("%s/Chi2_ndf.pdf", outputfilelocation.Data() ) );


  
  //// Plot Chi2/ndf
  //// Make arrays that TGraphErrors can use
  double x_ch_min[nRanges];
  double y_ch_min[nRanges];
  double x_ch_min_err[nRanges];
  double y_ch_min_err[nRanges];
  for (int rangeID = 0 ; rangeID < nRanges ; rangeID++)
    {
      x_ch_min[rangeID] = xMin_range[rangeID];
      y_ch_min[rangeID] = ChiSq_vector[rangeID] /ndf_vector[rangeID];
      x_ch_min_err[rangeID] = 0;
      y_ch_min_err[rangeID] = 0;
    }    
  TGraphErrors *Chi2_ndf_graph_xmin = new TGraphErrors(nRanges, x_ch_min, y_ch_min, x_ch_min_err, y_ch_min_err);
  utilityHandler.customizeGraphMore(Chi2_ndf_graph_xmin, 33, kBlue, 3,"","xMin","Chi^{2}/ndf");

  //// canvas
  TCanvas *xmin_chi2_ndf_canvas = new TCanvas("xmin_chi2_ndf_canvas","xmin_chi2_ndf_canvas",800,600);  
  utilityHandler.adjustCanvas(xmin_chi2_ndf_canvas);
  Chi2_ndf_graph_xmin ->Draw("AP");
  xmin_chi2_ndf_canvas->Update();
  xmin_chi2_ndf_canvas->SaveAs(Form("%s/Chi2_ndf_xmin.pdf", outputfilelocation.Data() ) );

  

  TCanvas *Rsf_hist_canvas = new TCanvas("Rsf_hist_canvas","Rsf_hist_canvas",800,600);
  Rsf_hist_canvas->SetGrid();
  utilityHandler.adjustCanvas(Rsf_hist_canvas);
  h_Rsf->Draw("hist");

  TCanvas *pull_canvas = new TCanvas("pull_canvas","pull_canvas",800,600);
  pull_canvas->SetGrid();
  utilityHandler.adjustCanvas(pull_canvas);
  h_pull->Draw("hist");


  int nHist = hist_vector_data.size();
  int nCols = 4;
  int nRows = (nHist + nCols - 1) / nCols;

  std::vector<TH1D*> overall_fit_as_histogram;

  TCanvas* fits_canvas = new TCanvas("fits_canvas", "fits_canvas", 1000, 600);
  fits_canvas->Divide(nCols, nRows); 
  for (int i = 0; i < nHist; ++i) {
    fits_canvas->cd(i + 1);
    
    hist_vector_data[i]->GetXaxis() ->SetRangeUser(xMin_range[i], xMax_range[i]);
    hist_vector_data[i]->Draw();

    int nEntries = hist_vector_data[i]->GetEntries();
    
    //// make a histogram that is the sum of the poly fit and the scaled mc histograms. 
    TH1D* sum_histo = utilityHandler.sumHistogramsWithPolynomial(hist_result_vector_p[i],hist_result_vector_n[i] , poly_fit_result[i]);
    overall_fit_as_histogram.push_back(sum_histo);

    double padHeight = fits_canvas->GetWh() / nRows;
    double legendTextSize = 0.015 * (600.0 / padHeight);  // Adjust based on canvas height

    // // Create and customize the legend
    TLegend* legend = new TLegend(0.5, 0.6, 0.9, 0.9);
    legend->SetTextSize(legendTextSize);  // Adjust size dynamically
    legend->SetMargin(0.10);  // Adjust margin to reduce space (default is around 0.25)
    // legend->AddEntry(hist_vector_data[i], hist_vector_data[i]->GetName(), "l");
    legend->AddEntry("", Form("R= %.4f +/- %.4f ", Rsf_vector[i],Rsf_err_vector[i]), "");
    legend->AddEntry("", Form("#chi^{2}/ndf = %.2f / %.0f  ", ChiSq_vector[i] ,ndf_vector[i]), "");
    legend->AddEntry("", Form("Entries: %i", nEntries ), "");
    legend->Draw();


    sum_histo ->SetLineColor(kRed);
    sum_histo->Draw("same");
    poly_fit_result[i]->SetLineColor(kCyan);
    poly_fit_result[i]->Draw("same");
    hist_result_vector_p[i] ->SetLineColor(kGreen);
    hist_result_vector_p[i]->Draw("same");
    hist_result_vector_n[i] ->SetLineColor(kMagenta);
    hist_result_vector_n[i]->Draw("same");
     
  }// end loop over slices
  fits_canvas->Update();
  fits_canvas->SaveAs(Form("%s/fitted_slices.pdf",outputfilelocation.Data() ) );

  

  std::vector<TH1D*> hist_residual_vector;
 
  for (int rangeID = 0 ; rangeID < nRanges ; rangeID++)
    {
      TH1D* hist_residual = (TH1D*)hist_vector_data[rangeID]->Clone(Form("hist_residual_%i",rangeID));  
      hist_residual  ->Add(overall_fit_as_histogram[rangeID], -1);
      hist_residual->GetXaxis() ->SetRangeUser(xMin_range[rangeID],xMax_range[rangeID]);
      hist_residual_vector.push_back(hist_residual);
    }// end loop over slices

   
  TCanvas* resid_canvas = new TCanvas("resid_canvas", "resid_canvas", 1000, 600);
  resid_canvas->Divide(nCols, nRows);
  for (int i = 0; i < nHist; ++i) {
    resid_canvas->cd(i + 1);
    hist_residual_vector[i]->Draw();
  }
  resid_canvas->Update();
  resid_canvas->SaveAs(Form("%s/residuals.pdf",outputfilelocation.Data() ) );

  
 
  //// extract the histogram title and print it
  std::string title = hist_1D_data->GetTitle();
  utilityHandler.printParsedTitle(title,outputfilelocation,"data");

  fout->Write();

  f1->Close();
  f2->Close();
  f3->Close();

}// End Main





