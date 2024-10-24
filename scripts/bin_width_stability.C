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
double xmin = -2;
double xmax = -1.4;
std::string AxisTitle ="Bin Width (m)";
std::vector<double> initialParams ={1,1,0,0,1,1,-1};



std::vector<double> binwidth_vector = {0.001,0.005,0.01,0.02,0.03,0.04,0.05};

std::vector<TString> hist_subtitle = {"0p001","0p005","0p01","0p02","0p03","0p04","0p05"}; 


void bin_width_stability(TString KineString="sbs4_30p", int PolyOrder_input = 2){ // main
 
  gStyle->SetNumberContours(255); 
  // gStyle->SetOptStat(0110);
  gStyle->SetOptStat(0);
  
  Utility utilityHandler; // class that gives us access to various functions to use between programs.
  
  FileNames fileNamesHandler;
 

  //// Histogram to store the Rsf distriubtion
  TH1D *h_Rsf = new TH1D("h_Rsf","h_Rsf",300, 0.8,1.2);

  /// Histogram for the pull distribution
  TH1D *h_pull = new TH1D("h_pull","h_pull",10,-1,1);

										       
   
  if (KineString == "sbs4_30p"){
    DataFileString =fileNamesHandler.DataFileString_sbs4_30p_BinStudy;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs4_30p_BinStudy;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs4_30p_BinStudy;
    xmin = fileNamesHandler.xmin_sbs4_30p;
    xmax = fileNamesHandler.xmax_sbs4_30p;
    initialParams =fileNamesHandler.initialParams_sbs4_30p;
  }else if (KineString == "sbs4_50p"){
    DataFileString = fileNamesHandler.DataFileString_sbs4_50p_BinStudy;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs4_50p_BinStudy;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs4_50p_BinStudy;
    xmin = fileNamesHandler.xmin_sbs4_50p;
    xmax = fileNamesHandler.xmax_sbs4_50p;
    initialParams =fileNamesHandler.initialParams_sbs4_50p;
  }
  else if (KineString == "sbs8_70p"){
    DataFileString = fileNamesHandler.DataFileString_sbs8_70p_BinStudy;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs8_70p_BinStudy;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs8_70p_BinStudy;
    xmin = fileNamesHandler.xmin_sbs8_70p;
    xmax = fileNamesHandler.xmax_sbs8_70p;
    initialParams =fileNamesHandler.initialParams_sbs8_70p;
  }
else if (KineString == "sbs8_50p"){
    DataFileString = fileNamesHandler.DataFileString_sbs8_50p_BinStudy;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs8_50p_BinStudy;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs8_50p_BinStudy;
    xmin = fileNamesHandler.xmin_sbs8_50p;
    xmax = fileNamesHandler.xmax_sbs8_50p;
    initialParams =fileNamesHandler.initialParams_sbs8_50p;
  }
  else if (KineString == "sbs8_100p"){
    DataFileString = fileNamesHandler.DataFileString_sbs8_100p_BinStudy;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs8_100p_BinStudy;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs8_100p_BinStudy;
    xmin = fileNamesHandler.xmin_sbs8_100p;
    xmax = fileNamesHandler.xmax_sbs8_100p;
    initialParams =fileNamesHandler.initialParams_sbs8_100p;
  }else if (KineString == "sbs9_70p"){
    DataFileString = fileNamesHandler.DataFileString_sbs9_70p_BinStudy;
    ProtonFileString =fileNamesHandler.ProtonFileString_sbs9_70p_BinStudy;
    NeutronFileString = fileNamesHandler.NeutronFileString_sbs9_70p_BinStudy;
    xmin = fileNamesHandler.xmin_sbs9_70p;
    xmax = fileNamesHandler.xmax_sbs9_70p;
    initialParams =fileNamesHandler.initialParams_sbs9_70p;
  } else {
    std::cout<<"Error with kinematic setting"<<std::endl;
    return;
  }


  if((PolyOrder_input < 0) && PolyOrder_input > 6){
    cout<<"error with poly order input. Setting to 2."<<endl;
    PolyOrder = 2;
  }else PolyOrder=PolyOrder_input; 
 
  cout<<"What we're running: "<<KineString<<endl;
  cout<<DataFileString<<endl;
  cout<<ProtonFileString<<endl;
  cout<<NeutronFileString<<endl;
  cout<<HistogramName<<endl;
  cout<<"poly order: "<<PolyOrder <<endl;
 
      
      
  
  TFile *f1 = TFile::Open(DataFileString); // data
  TFile *f2 = TFile::Open(ProtonFileString); //proton
  TFile *f3 = TFile::Open(NeutronFileString); //neutron


  // set location and name for output file 
  TString outputfilelocation="../output/stability/"+KineString+"/binWidth";
  TString outputfilename = outputfilelocation +"/"+ KineString+ ".root";

  // Declare outfile
  TFile *fout = new TFile(outputfilename,"RECREATE");
  cout<<"writing to file: "<< outputfilename <<endl;


  std::vector<TH1D*> hist_vector_data; 
  std::vector<TH1D*> hist_vector_p; 
  std::vector<TH1D*> hist_vector_n;

  TH1D *hist_data_orig;
  TH1D *hist_proton_orig;
  TH1D *hist_neutron_orig;

  for(int i = 0; i< binwidth_vector.size();i++){
    TString HistogramFullName = HistogramName +"_"+  hist_subtitle[i];
   hist_data_orig = (TH1D*)f1->Get(HistogramFullName.Data());
    hist_proton_orig = (TH1D*)f2->Get(HistogramFullName.Data());
    hist_neutron_orig  =  (TH1D*)f3->Get(HistogramFullName.Data());
    // Make clones of the histograms
    TH1D *hist_data = (TH1D*)hist_data_orig->Clone("hist_data");
    TH1D *hist_p = (TH1D*)hist_proton_orig->Clone("hist_p"); 
    TH1D *hist_n = (TH1D*)hist_neutron_orig->Clone("hist_n");
    hist_vector_data.push_back(hist_data);
    hist_vector_p.push_back(hist_p);
    hist_vector_n.push_back(hist_n);
  }
  


    
  std::vector<FitHistogram> fitHandler_vector; // vector of objects fo the FitHistogram class that will help us do data-mc comparion. .
  



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

 
 for (int histid = 0 ; histid < binwidth_vector.size() ; histid++)
    {
      // this uses clones of the histograms
      TH1D *hist_p  =(TH1D*)hist_vector_p[histid]->Clone("hist_p");
      TH1D *hist_n  =(TH1D*)hist_vector_n[histid]->Clone("hist_n");
      TH1D *hist_data_dont_draw  =(TH1D*)hist_vector_data[histid]->Clone("hist_data_dont_draw");
      TH1D *hist_data  =(TH1D*)hist_vector_data[histid]->Clone("hist_data");

      if (!hist_data_dont_draw || !hist_p || !hist_n) {
	std::cerr << "Error: One of the histograms is null!" << std::endl;
	return;
      }

      FitHistogram fitHandler(hist_p,hist_n,xmin,xmax);

      fitHandler.setPolyOrder(PolyOrder);// setting the order for the polynomial background
      
      fitHandler.fitDataPoly(hist_data_dont_draw,initialParams);
      
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
     

      // //cout<<"shift_p in the loop: "<<shift_p<< " +/- "<<shift_p_err<<endl;
      // //cout<<"shift_n in the loop: "<<shift_n<< " +/- "<<shift_n_err<<endl;

      std::vector<double> poly_result;
      std::vector<double> poly_result_err;
      for (int i =0 ; i <=PolyOrder; i++)
	{
	  poly_result.push_back(fitHandler.poly_result[i]);
	  poly_result_err.push_back(fitHandler.poly_result_err[i]);
	}

      h_Rsf->Fill(Rsf);
      
      cout<<"histid: "<<histid<<endl;
      cout<<"scale_p = "<< scale_p <<endl;
      cout<<"Rsf = "<< Rsf <<" +/- "<< Rsf_err<<endl;
      cout<<"shift p: "<<shift_p<<" +/- "<<shift_p_err<<", shift_n: "<<shift_n<<"+/-"<<shift_n_err<<endl;
      cout<<"chi^2/ndf = "<<ChiSq<<"/"<<ndf<<" = "<<ChiSq/ndf<<endl;
      cout<<endl;

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
      // poly_result_vector_of_arrays.push_back(poly_result);
      poly_result_vector_of_vectors.push_back(poly_result);
      poly_result_err_vector_of_vectors.push_back(poly_result_err);
	  
      // delete Fit; //
    }// end loop over slices


 //     // Print the contents of the vector of vectors
 //  for (const auto& vec : poly_result_vector_of_vectors) {
 //    for (double val : vec) {
 //      // std::cout << val << " ";
 //    }
 //    std::cout << std::endl;
 //  }
  
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
  for (int i = 0 ; i < Rsf_vector.size() ; i++)
    {	    
      // utility function that scales and shifts the provided histogram and returns a new histogram 
      hist_temp_p = utilityHandler.ScaleAndShiftHistogram(hist_vector_p[i], scale_p_vector[i], shift_p_vector[i]);
      hist_temp_n = utilityHandler.ScaleAndShiftHistogram(hist_vector_n[i], scale_n_vector[i], shift_n_vector[i]);
	 
      // save histograms to the global vectors 
      hist_result_vector_p.push_back(hist_temp_p);
      hist_result_vector_n.push_back(hist_temp_n);	
    }// end loop over slices


  /// Make a canvas called c1 to overwrite the default canvas and avoid crashes
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  hist_vector_p[0]->Draw("colz");
  
  
  std::vector<TF1*> poly_fit_result;
  // loop over slices and make fits to plot the polynomial background result
  for (int i = 0 ; i < Rsf_vector.size() ; i++)
    {
      // get the fit results from the vector of vectors
      std::vector<double> poly_params_vector;
      poly_params_vector = poly_result_vector_of_vectors[i];

      // convert the vector into an array that TF1 can use
      int size = poly_params_vector.size();
      double poly_params_array[size]; 
      for (int i  = 0 ; i < size; i++){
	poly_params_array[i] = poly_params_vector[i];
      }

      TF1 *fit = new TF1(Form("poly_%i_",i),Form("pol%i",PolyOrder), xmin, xmax);
      fit->SetParameters(poly_params_array);
      fit->SetNpx(500);
      poly_fit_result.push_back(fit);
      
      if (utilityHandler.doesFunctionGoBelowZero(fit,xmin,xmax))
	{
	  // background function went below zero, which isn't physical
	  std::cout<<"On slice ID: "<<i<<std::endl;
	  std::cout<<"---------------------------"<<std::endl;
	}

    }// end loop over slices


  
  //// Plot Rsf
     //// Make arrays that TGraphErrors can use 
  double x[Rsf_vector.size()];
  double y[Rsf_vector.size()];
  double x_err[Rsf_vector.size()];
  double y_err[Rsf_vector.size()];
  for (int i = 0 ; i < Rsf_vector.size() ; i++)
    {
      //double xCenter =8 / static_cast<double>(binwidth_vector[i]);
      double xCenter =binwidth_vector[i];
      double xWidth = 0;
     
      x[i] = xCenter;
      y[i] = Rsf_vector[i];
      x_err[i] = xWidth;
      y_err[i] = Rsf_err_vector[i];

      //cout<< "i = "<<i<<", xCenter = "<<x[i]<<", x= "<<x[i]<<", y = "<<y[i]<<", x_err= "<<x_err[i]<<", y_err = "<<y_err[i]<<endl;
    }    
  TGraphErrors *Rsf_graph = new TGraphErrors(Rsf_vector.size(), x, y, x_err, y_err);
  utilityHandler.customizeGraphMore(Rsf_graph, 33, kBlue, 3,"",AxisTitle,"Rsf",1.4,1.4);

  // Fit a straight line to the graph
  TF1* fit_Rsf_graph = new TF1("fit_Rsf_graph", "pol1",Rsf_graph ->GetX()[0], Rsf_graph->GetX()[Rsf_graph->GetN()-1]);
  Rsf_graph->Fit(fit_Rsf_graph, " R");
  
  //// canvas
  TCanvas *graphcanvas = new TCanvas("graphcanvas","graphcanvas",1200,600);
  // utilityHandler.adjustCanvas(graphcanvas);
  graphcanvas->SetGrid();
  graphcanvas->Divide(2,1);
  graphcanvas->cd(1);
  utilityHandler.adjustCanvas(graphcanvas);
  Rsf_graph->Draw("AP");

  // Get fit parameters
  double constant = fit_Rsf_graph->GetParameter(0);       // The constant value
  double constantError = fit_Rsf_graph->GetParError(0);    // The error on the constant
   double slope = fit_Rsf_graph->GetParameter(1);       
  double slopeError = fit_Rsf_graph->GetParError(1);    
  double chi2 = fit_Rsf_graph->GetChisquare();             // The chi-squared value
  int ndf = fit_Rsf_graph->GetNDF();                       // The number of degrees of freedom
  double chi2_ndf = chi2 / ndf;                          // chi2/ndf

  // Draw the graph with the fit result again
  Rsf_graph->Draw("AP");

  // Use TLatex to add the fit result and chi²/ndf to the canvas
  TLatex latex;
  latex.SetNDC();  // Use normalized coordinates (0 to 1)
  latex.SetTextSize(0.04);
  latex.DrawLatex(0.2, 0.85, Form("Fit: y = %.5f x + %.5f", slope, constant));
  latex.DrawLatex(0.2, 0.80, Form("#chi^{2}/ndf = %.2f/%i = %.2f", chi2,ndf,chi2_ndf));
  latex.DrawLatex(0.2, 0.2, Form("Rsf Mean = %.5f #pm StDev = %.5f", Rsf_mean, Rsf_stdev));
  
  graphcanvas->Update();


//// Plot Chi2/ndf
  //// Make arrays that TGraphErrors can use
  
  double x_ch[Rsf_vector.size()];
  double y_ch[Rsf_vector.size()];
  double x_ch_err[Rsf_vector.size()];
  double y_ch_err[Rsf_vector.size()];
  for (int i = 0 ; i < Rsf_vector.size() ; i++)
    {
      double xCenter = binwidth_vector[i];
      double xWidth = 0;
      
      x_ch[i] = xCenter;
      y_ch[i] = ChiSq_vector[i] /ndf_vector[i];
      x_ch_err[i] = xWidth;
      y_ch_err[i] = 0;
    }    
  TGraphErrors *Chi2_ndf_graph = new TGraphErrors(Rsf_vector.size(), x_ch, y_ch, x_ch_err, y_ch_err);
  utilityHandler.customizeGraphMore(Chi2_ndf_graph, 33, kBlue, 3,"",AxisTitle,"Chi^{2}/ndf",1.4,1.4);
  graphcanvas->cd(2);
  Chi2_ndf_graph ->Draw("AP");
  graphcanvas->Update();
  
  graphcanvas->SaveAs(Form("%s/Rsf.pdf",outputfilelocation.Data() ) );

  
  TCanvas *Rsf_hist_canvas = new TCanvas("Rsf_hist_canvas","Rsf_hist_canvas",800,600);
  Rsf_hist_canvas->SetGrid();
  utilityHandler.adjustCanvas(Rsf_hist_canvas);
  h_Rsf->Draw("hist");

  TCanvas *pull_canvas = new TCanvas("pull_canvas","pull_canvas",800,600);
  pull_canvas->SetGrid();
  utilityHandler.adjustCanvas(pull_canvas);
  h_pull->Draw("hist");
  


  int nHist = Rsf_vector.size();
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
    hist_result_vector_p[i]->Draw("same hist");
    hist_result_vector_n[i] ->SetLineColor(kMagenta);
    hist_result_vector_n[i]->Draw("same hist");
    sum_histo ->SetLineColor(kRed);
    sum_histo->Draw("same");
     
  }// end loop over slices
  fits_canvas->Update();
  fits_canvas->SaveAs(Form("%s/fitted_slices.pdf",outputfilelocation.Data()));
  

  std::vector<TH1D*> hist_residual_vector;
 
  for (int i = 0 ; i < Rsf_vector.size() ; i++)
    {
      TH1D* hist_residual = (TH1D*)hist_vector_data[i]->Clone(Form("hist_residual_%d",i));  
      hist_residual  ->Add(overall_fit_as_histogram[i], -1);
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

  
  fout->Write(); 

  f1->Close();
  f2->Close();
  f3->Close();
  
}// End Main





