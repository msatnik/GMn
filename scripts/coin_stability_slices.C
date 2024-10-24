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
//std::vector<TH1D*> hist_vector_p; 
//std::vector<TH1D*> hist_vector_n; 

TH1D* hist_result_p;
TH1D* hist_resut_n;

std::vector<TH1D*> hist_result_vector_p; 
std::vector<TH1D*> hist_result_vector_n; 

std::vector<TH1D*> hist_residual_vector; 


//// default values to be overwritten 
TString DataFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_tight_2Dhistos.root";
TString ProtonFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_tight_2Dhistos.root";
TString NeutronFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_tight_2Dhistos.root";
TString HistogramName = "hcal_dx__hcal_sh_atime_diff";
std::string AxisTitle ="hcal time - sh time";
TString note = "";
std::string xMin_xMax_string ="(-15,-13),(-13,-11),(-11,-9),(-9,-7),(-7,-5),(-5,-3),(-3,-1),(-1,1),(1,3),(3,5),(5,7),(7,9),(9,11),(11,13),(13,15)"; 
int PolyOrder = 2;
double xmin = -2.1; // -1.8
double xmax = 1.4;//1;
// for if we are going to use the range incrementers. Default to be overwriten. 
double left=-0.055;
double right = 0.055;
double plus_or_minus_xmin = 0.01;
double plus_or_minus_xmax = 0;
double stepsize = 0.001;
int nDecimals = 3;


///******* coin time *************************************************************
//std::string xMin_xMax_string_coin ="(-5,-4),(-4,-3),(-3,-2),(-2,-1),(-1,0),(0,1),(1,2),(2,3),(3,4),(4,5)";
std::string xMin_xMax_string_coin ="(-5,-4.5),(-4.5,-4),(-4,-3.5),(-3.5,-3),(-3,-2.5),(-2.5,-2),(-2,-1.5),(-1.5,-1),(-1,-0.5),(-0.5,0),(0,0.5),(0.5,1),(1,1.5),(1.5,2),(2,2.5),(2.5,3),(3,3.5),(3.5,4),(4,4.5),(4.5,5)";

//*********Coin Time Range Study ***************************************************************
double left_coin=-10;
double right_coin = 10;
double plus_or_minus_coin =1;
double stepsize_coin = 0.1;
int nDecimals_coin = 1;


void coin_stability_slices(TString KineString="sbs4_30p", TString VarString = "coin", int PolyOrder_input=2){ // main
  // bit of a test for now. Will need to make this more sophisticated in the future. 

  gStyle->SetNumberContours(255); 
  gStyle->SetOptStat(0110);
  
  FileNames fileNamesHandler;

 
  //std::vector<int> color_vector = {kRed,kBlue,kGreen,kMagenta,kCyan,kYellow};
  std::vector<int> color_vector = {kRed,kRed};


  //// Histogram to store the Rsf distriubtion
  TH1D *h_Rsf = new TH1D("h_Rsf","h_Rsf",300, 0.8,1.2);

  /// Histogram for the pull distribution
  TH1D *h_pull = new TH1D("h_pull","h_pull",10,-1,1);

  Utility utilityHandler; // class that gives us access to various functions to use between programs.

  if (VarString == "coin"){
    HistogramName = "hcal_dx__hcal_sh_atime_diff";
    AxisTitle ="hcal - sh time ";
    xMin_xMax_string = xMin_xMax_string_coin;
  }else if (VarString == "coin_range_left"){
    HistogramName = "hcal_dx__hcal_sh_atime_diff";
    AxisTitle ="hcal -  sh time (ns)";
    xMin_xMax_string  = utilityHandler.incrementRangeStringStepsize(left_coin, right_coin,plus_or_minus_coin, 0,stepsize_coin, nDecimals_coin);
  }else if (VarString == "coin_range_right"){
    HistogramName = "hcal_dx__hcal_sh_atime_diff";
    AxisTitle ="hcal - sh time (ns)";
    xMin_xMax_string = utilityHandler.incrementRangeStringStepsize(left_coin,right_coin,0, plus_or_minus_coin,stepsize_coin, nDecimals_coin);
  }else {
    std::cout<<"Error in the Varibale String "<<std::endl;
    return;
  }

  if (KineString == "sbs4_30p"){
    DataFileString = fileNamesHandler.DataFileString_sbs4_30p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs4_30p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs4_30p;
    xmin = fileNamesHandler.xmin_sbs4_30p;
    xmax = fileNamesHandler.xmax_sbs4_30p;
    //PolyOrder = PolyOrder_sbs4_30p;
  }else if (KineString == "sbs4_50p"){
    DataFileString = fileNamesHandler.DataFileString_sbs4_50p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs4_50p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs4_50p;
    xmin = fileNamesHandler.xmin_sbs4_50p;
    xmax = fileNamesHandler.xmax_sbs4_50p;
    //PolyOrder = PolyOrder_sbs4_50p;
  }else if (KineString == "sbs8_70p"){
    DataFileString = fileNamesHandler.DataFileString_sbs8_70p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs8_70p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs8_70p;
    xmin = fileNamesHandler.xmin_sbs8_70p;
    xmax = fileNamesHandler.xmax_sbs8_70p;
    //PolyOrder = PolyOrder_sbs8_70p;
  }else if (KineString == "sbs8_50p"){
    DataFileString = fileNamesHandler.DataFileString_sbs8_50p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs8_50p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs8_50p;
    xmin = fileNamesHandler.xmin_sbs8_50p;
    xmax = fileNamesHandler.xmax_sbs8_50p;
    //PolyOrder = PolyOrder_sbs8_50p;
  }else if (KineString == "sbs8_100p"){
    DataFileString = fileNamesHandler.DataFileString_sbs8_100p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs8_100p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs8_100p;
    xmin = fileNamesHandler.xmin_sbs8_100p;
    xmax = fileNamesHandler.xmax_sbs8_100p;
    //PolyOrder = PolyOrder_sbs8_100p;
  }else if (KineString == "sbs9_70p"){
    DataFileString = fileNamesHandler.DataFileString_sbs9_70p;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs9_70p;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs9_70p;
    xmin = fileNamesHandler.xmin_sbs9_70p;
    xmax = fileNamesHandler.xmax_sbs9_70p;
    // PolyOrder = PolyOrder_sbs9_70p;
  }  else {
    std::cout<<"Error with kinematic setting"<<std::endl;
    return;
  }

  if((PolyOrder_input < 0) && PolyOrder_input > 6){
    cout<<"error with poly order input. Setting to 2."<<endl;
    PolyOrder = 2;
  }else PolyOrder=PolyOrder_input;

  cout<<"What we're running: "<<KineString<< ", "<<VarString<<endl;
  cout<<DataFileString<<endl;
  cout<<ProtonFileString<<endl;
  cout<<NeutronFileString<<endl;
  cout<<xmin<< ", "<<xmax<<endl;
  cout<<HistogramName<<endl;
  cout<<AxisTitle<<endl;
  cout<<"poly order: "<<PolyOrder <<endl;
  // cout<<xMin_xMax_string<<endl;

  TFile *f1 = TFile::Open(DataFileString); // data
  TFile *f2 = TFile::Open(ProtonFileString); //proton
  TFile *f3 = TFile::Open(NeutronFileString); //neutron

  


  //// set location and name for root output file 
  TString outputfilelocation="../output/stability/"+KineString+"/"+VarString;
  TString outputfilename = outputfilelocation +"/"+ KineString+ "_" +VarString+".root";
  

  
  TFile *fout = new TFile(outputfilename,"RECREATE");
  cout<<"writing to file: "<< outputfilename <<endl;


  // Load Histograms
  TH2D *hist_data_orig = (TH2D*)f1->Get(HistogramName);
  TH1D *hist_proton_orig_1d  = (TH1D*)f2->Get("hcal_dx_1d_allcuts");
  TH1D *hist_neutron_orig_1d  =  (TH1D*)f3->Get("hcal_dx_1d_allcuts");

  // Make clones of the histograms
  TH2D *hist_2D_data = (TH2D*)hist_data_orig->Clone("hist_2D_data");
  TH1D* hist_p = (TH1D*)hist_proton_orig_1d->Clone("hist_p"); // global histogram
  TH1D* hist_n = (TH1D*)hist_neutron_orig_1d->Clone("hist_n"); // global histogram

 

  // Set up the slices for the Y projections. 
  // std::vector<double> xSlices = {-8,-5,-3,-2.5,-2,-1.5,-1,-0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3,5,8};


  std::vector<double> xMinSlices;
 std:vector<double> xMaxSlices;
  utilityHandler.parseStringToVectors(xMin_xMax_string, xMinSlices, xMaxSlices);
  
  cout<<endl;
  cout<< "nSlices: "<<xMinSlices.size()<<endl;
  cout<<endl;

  if (xMinSlices.size() != xMaxSlices.size())
    {
      std::cout<<"Error: max and min are different sizes"<<std::endl;
      return;
    }

  std::cout<<"Slices:"<<std::endl;
  for (int i = 0;i<xMinSlices.size();i++)
    {
      std::cout<<"("<<xMinSlices[i]<<", "<<xMaxSlices[i]<<")"<<", ";
    }
  std::cout<<std::endl;

  
  std::vector<FitHistogram> fitHandler_vector; // vector of objects fo the FitHistogram class that will help us do data-mc comparion.

 

  cout<<endl<<"Slicing Data"<<endl;
  utilityHandler.SliceAndProjectHistogram_xMinxMax_inclusiveMin(hist_2D_data, xMinSlices,xMaxSlices, hist_vector_data, AxisTitle,"dx","data");

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

 
  double initialParameters[7]={1,1,0,0,1,1,-1};

      
  
  /// Fit each slice 
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {
      /// set up the global histograms that the fit function is going to use
      // this uses clones of the histograms
      TH1D *hist_data_dont_draw  =(TH1D*)hist_vector_data[sliceid]->Clone("hist_data_dont_draw");
      TH1D *hist_data  =(TH1D*)hist_vector_data[sliceid]->Clone("hist_data");

      if (!hist_data_dont_draw || !hist_data) {
	std::cerr << "Error: One of the histograms is null!" << std::endl;
	return;
      }

      FitHistogram fitHandler(hist_p,hist_n,xmin,xmax);

      fitHandler.setPolyOrder(PolyOrder);// setting the order for the polynomial background
      
      fitHandler.fitDataPoly(hist_data_dont_draw);

      TF1 *fit_result = (TF1*)fitHandler.fitFunc->Clone("fit_result");// bug in here that corrupts it and crashes if you try and draw

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
      for (int i =0 ; i <=PolyOrder; i++)
	{
	  poly_result.push_back(fitHandler.poly_result[i]);
	  poly_result_err.push_back(fitHandler.poly_result_err[i]);
	}

      h_Rsf->Fill(Rsf);

      //cout<<"sliceid: "<<sliceid<<endl;
      // cout<<"scale_p = "<< scale_p <<endl;
      // cout<<"Rsf = "<< Rsf <<" +/- "<< Rsf_err<<endl;
      // cout<<"shift p: "<<shift_p<<" +/- "<<shift_p_err<<", shift_n: "<<shift_n<<"+/-"<<shift_n_err<<endl;
      // cout<<endl;


      //// save results into the vectors 
      fit_vector.push_back(fit_result);
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
   
      poly_result_vector_of_vectors.push_back(poly_result);
      poly_result_err_vector_of_vectors.push_back(poly_result_err);
	  
    }// end loop over slices


     // // Print the contents of the vector of vectors
     // for (const auto& vec : poly_result_vector_of_vectors) {
     //   for (double val : vec) {
     // 	 std::cout << val << " ";
     //   }
     //   std::cout << std::endl;
     // }

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
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {	    
      // utility function that scales and shifts the provided histogram and returns a new histogram 
      hist_temp_p = utilityHandler.ScaleAndShiftHistogram(hist_p, scale_p_vector[sliceid], shift_p_vector[sliceid]);
      hist_temp_n = utilityHandler.ScaleAndShiftHistogram(hist_n, scale_n_vector[sliceid], shift_n_vector[sliceid]);
	 
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

      TF1 *fit = new TF1(Form("poly_%i_",sliceid),Form("pol%i",PolyOrder), xmin, xmax);
      fit->SetParameters(poly_params_array);
      fit->SetNpx(500);
      poly_fit_result.push_back(fit);
      
      if (utilityHandler.doesFunctionGoBelowZero(fit,xmin,xmax))
	{
	  // background function went below zero, which isn't physical
	  std::cout<<"On slice ID: "<<sliceid<<std::endl;
	  std::cout<<"---------------------------"<<std::endl;
	}

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
  utilityHandler.customizeGraphMore(Rsf_graph, 33, kBlue, 3,"",AxisTitle,"Rsf");

  // Fit a straight line to the graph
  TF1* fit_Rsf_graph = new TF1("fit_Rsf_graph", "[0]",Rsf_graph ->GetX()[0], Rsf_graph->GetX()[Rsf_graph->GetN()-1]);
  Rsf_graph->Fit(fit_Rsf_graph, "Q R");
  
  //// canvas
  TCanvas *graphcanvas = new TCanvas("graphcanvas","graphcanvas",800,600);
  graphcanvas->SetGrid();
  utilityHandler.adjustCanvas(graphcanvas);
  Rsf_graph->Draw("AP");

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
  latex.DrawLatex(0.2, 0.80, Form("#chi^{2}/ndf = %.2f/%i = %.2f", chi2,ndf,chi2_ndf));
  latex.DrawLatex(0.8, 0.2, Form("Rsf Mean = %.5f #pm StDev = %.5f", Rsf_mean, Rsf_stdev));
  
  graphcanvas->Update();
  graphcanvas->SaveAs(Form("%s/Rsf_xCenter.pdf",outputfilelocation.Data() ) );

  //// Plot Chi2/ndf
  //// Make arrays that TGraphErrors can use
  
  double x_ch[xMinSlices.size()];
  double y_ch[xMinSlices.size()];
  double x_ch_err[xMinSlices.size()];
  double y_ch_err[xMinSlices.size()];
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {
      double xCenter = (xMaxSlices[sliceid] + xMinSlices[sliceid] )/2;
      double xWidth = xMaxSlices[sliceid] - xMinSlices[sliceid];
      
      x_ch[sliceid] = xCenter;
      y_ch[sliceid] = ChiSq_vector[sliceid] /ndf_vector[sliceid];
      x_ch_err[sliceid] = xWidth;
      y_ch_err[sliceid] = 0;
    }    
  TGraphErrors *Chi2_ndf_graph = new TGraphErrors(xMinSlices.size(), x_ch, y_ch, x_ch_err, y_ch_err);
  utilityHandler.customizeGraphMore(Chi2_ndf_graph, 33, kBlue, 3,"","xMax","Chi^{2}/ndf");


  //// canvas
  TCanvas *chi2_ndf_canvas = new TCanvas("chi2_ndf_canvas","chi2_ndf_canvas",800,600);  
  utilityHandler.adjustCanvas(chi2_ndf_canvas);
  Chi2_ndf_graph ->Draw("AP");
  chi2_ndf_canvas->Update();
  chi2_ndf_canvas->SaveAs(Form("%s/Chi2_ndf.pdf",outputfilelocation.Data()));
  
  TCanvas *Rsf_hist_canvas = new TCanvas("Rsf_hist_canvas","Rsf_hist_canvas",800,600);
  Rsf_hist_canvas->SetGrid();
  utilityHandler.adjustCanvas(Rsf_hist_canvas);
  h_Rsf->Draw("hist");

  TCanvas *pull_canvas = new TCanvas("pull_canvas","pull_canvas",800,600);
  pull_canvas->SetGrid();
  utilityHandler.adjustCanvas(pull_canvas);
  h_pull->Draw("hist");

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
  utilityHandler.customizeGraphMore(nEntries_graph, 33, kBlue, 3,"",AxisTitle,"nEntries");
  //// canvas
  TCanvas *nEntriesCanvas = new TCanvas("nEntriesCanvas","nEntriesCanvas",800,600);  nEntriesCanvas->SetGrid();
  utilityHandler.adjustCanvas(nEntriesCanvas);
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
  utilityHandler.customizeGraphMore(Rsf_xMax_graph, 33, kBlue, 3,"",AxisTitle,"Rsf");
  
  // Fit a straight line to the graph
  TF1* fit_Rsf_xMax_graph = new TF1("fit_Rsf_xMax_graph", "[0]",Rsf_xMax_graph ->GetX()[0], Rsf_xMax_graph->GetX()[Rsf_xMax_graph->GetN()-1]);
  Rsf_xMax_graph->Fit(fit_Rsf_xMax_graph, "Q R");
  
  //// canvas
  TCanvas *Rsf_xMaxCanvas = new TCanvas("Rsf_xMaxCanvas","Rsf_xMaxCanvas",800,600);
  Rsf_xMaxCanvas->SetGrid();
  utilityHandler.adjustCanvas(Rsf_xMaxCanvas);
  Rsf_xMax_graph->Draw("AP");
  Rsf_xMaxCanvas->Update();

  // Get fit parameters
  double constant_xMax = fit_Rsf_xMax_graph->GetParameter(0);       // The constant value
  double constantError_xMax = fit_Rsf_xMax_graph->GetParError(0);    // The error on the constant
  double chi2_xMax = fit_Rsf_xMax_graph->GetChisquare();             // The chi-squared value
  int ndf_xMax = fit_Rsf_xMax_graph->GetNDF();                       // The number of degrees of freedom
  double chi2_ndf_xMax = chi2 / ndf;                          // chi2/ndf

  // Draw the graph with the fit result again
  Rsf_xMax_graph->Draw("AP");

  // Use TLatex to add the fit result and chi²/ndf to the canvas
  TLatex latex_xMax;
  latex_xMax.SetNDC();  // Use normalized coordinates (0 to 1)
  latex_xMax.SetTextSize(0.04);
  latex_xMax.DrawLatex(0.2, 0.85, Form("Fit: y = %.5f #pm %.5f", constant_xMax, constantError_xMax));
  latex_xMax.DrawLatex(0.2, 0.80, Form("#chi^{2}/ndf = %.2f/%i = %.2f", chi2_xMax,ndf_xMax,chi2_ndf_xMax));
  latex_xMax.DrawLatex(0.2, 0.2, Form("Rsf Mean = %.5f #pm StDev = %.5f", Rsf_mean, Rsf_stdev));

  // Update the canvas
  Rsf_xMaxCanvas ->Update();
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
  utilityHandler.customizeGraphMore(Rsf_xMin_graph, 33, kBlue, 3,"",AxisTitle,"Rsf");

  // Fit a straight line to the graph
  TF1* fit_Rsf_xMin_graph = new TF1("fit_Rsf_xMin_graph", "[0]",Rsf_xMin_graph ->GetX()[0], Rsf_xMin_graph->GetX()[Rsf_xMin_graph->GetN()-1]);
  Rsf_xMin_graph->Fit(fit_Rsf_xMin_graph, "Q R");
  
  //// canvas
  TCanvas *Rsf_xMinCanvas = new TCanvas("Rsf_xMinCanvas","Rsf_xMinCanvas",800,600);
  Rsf_xMinCanvas->SetGrid();
  utilityHandler.adjustCanvas(Rsf_xMinCanvas);
  Rsf_xMin_graph->Draw("AP");
  Rsf_xMinCanvas->Update();
  Rsf_xMinCanvas->SaveAs(Form("%s/Rsf_xMin.pdf",outputfilelocation.Data()));
  // Get fit parameters
  double constant_xMin = fit_Rsf_xMin_graph->GetParameter(0);       // The constant value
  double constantError_xMin = fit_Rsf_xMin_graph->GetParError(0);    // The error on the constant
  double chi2_xMin = fit_Rsf_xMin_graph->GetChisquare();             // The chi-squared value
  int ndf_xMin = fit_Rsf_xMin_graph->GetNDF();                       // The number of degrees of freedom
  double chi2_ndf_xMin = chi2 / ndf;                          // chi2/ndf

  // Draw the graph with the fit result again
  Rsf_xMin_graph->Draw("AP");

  // Use TLatex to add the fit result and chi²/ndf to the canvas
  TLatex latex_xMin;
  latex_xMin.SetNDC();  // Use normalized coordinates (0 to 1)
  latex_xMin.SetTextSize(0.04);
  latex_xMin.DrawLatex(0.2, 0.85, Form("Fit: y = %.5f #pm %.5f", constant_xMin, constantError_xMin));
  latex_xMin.DrawLatex(0.2, 0.80, Form("#chi^{2}/ndf = %.2f/%i = %.2f", chi2_xMin,ndf_xMin,chi2_ndf_xMin));
  latex_xMin.DrawLatex(0.2, 0.2, Form("Rsf Mean = %.5f #pm StDev = %.5f", Rsf_mean, Rsf_stdev));
  

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

    int nEntries = hist_vector_data[i]->GetEntries();
    //// make a histogram that is the sum of the poly fit and the scaled mc histograms. 
    TH1D* sum_histo = utilityHandler.sumHistogramsWithPolynomial(hist_result_vector_p[i],hist_result_vector_n[i] , poly_fit_result[i]);
    overall_fit_as_histogram.push_back(sum_histo);

    double padHeight = fits_canvas->GetWh() / nRows;
    double legendTextSize = 0.015 * (600.0 / padHeight);  // Adjust based on canvas height

    // // Create and customize the legend
    TLegend* legend = new TLegend(0.5, 0.5, 0.9, 0.9);
    legend->SetTextSize(legendTextSize);  // Adjust size dynamically
    legend->SetMargin(0.10);  // Adjust margin to reduce space (default is around 0.25)
    // legend->AddEntry(hist_vector_data[i], hist_vector_data[i]->GetName(), "l");
    legend->AddEntry("", Form("R= %.4f +/- %.4f ", Rsf_vector[i],Rsf_err_vector[i]), "");
    legend->AddEntry("", Form("#chi^{2}/ndf = %.2f / %.0f  ", ChiSq_vector[i] ,ndf_vector[i]), "");
    legend->AddEntry("", Form("Entries: %i", nEntries ), "");
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
      TH1D* hist_residual = (TH1D*)hist_vector_data[sliceid]->Clone(Form("hist_residual_%d",sliceid));  
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
 
  

  // Create a canvas to draw the histogram
  TCanvas* datacanvas1d = new TCanvas("datacanvas1d", "datacanvas1d", 800, 600);
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

  datacanvas1d ->SaveAs(Form("%s/visulaized_cuts_1D_data.pdf",outputfilelocation.Data() ));
 
  // //// extract the histogram title and print it
  std::string title = hist_2D_data->GetTitle();
  cout<<"data title"<<endl;
  utilityHandler.printParsedTitle(title,outputfilelocation,"data");
  std::string proton_title = hist_p->GetTitle();
  cout<<"proton title: "<<endl;
  cout<< proton_title<<endl;
  std::string neutron_title =hist_n->GetTitle();
  cout<<"neutron title: "<<endl;
  cout<< neutron_title<<endl;


  /// Make a canvas called c1 to overwrite the default canvas and avoid crashes
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  hist_2D_data->Draw("colz");
  
  
  
  fout->Write(); 

  f1->Close();
  f2->Close();
  f3->Close();
  
}// End Main


