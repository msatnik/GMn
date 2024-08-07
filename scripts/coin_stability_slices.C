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

std::string output_dir = "testpdfs";


// these global histograms are going to temporarily hold the proton and neutron slices for when I call the fit function. 
TH1D *hist_p;
TH1D *hist_n;


 double xmin = -3; // -1.8
  double xmax = 2.5;//1;


//// struct to store everything for each slice
struct SliceResults{
  //// histograms
  //TH1D *hist_data; // dx from data
  // TH1D *hist_proton; // dx from proton mc
  //TH1D *hist_neutron; // dx from neutron mc

  //// statistics
  int nEntries_data = 0;
  int nEntries_proton = 0;
  int nEntries_neutron = 0;

  //// Fits. Maybe multuple fits? Maybe start with 4th order poly.
  TF1 *overall_fit;

};


// Functions 
void SliceAndProjectHistogram_xMaxFixed(TH2D* hist2D, const std::vector<double>& xSlices, std::vector<TH1D*>& histVector, double xMax, std::string xAxisName, std::string yAxisName, std::string type);
void SliceAndProjectHistogram_slices(TH2D* hist2D, const std::vector<double>& xSlices, std::vector<TH1D*>& histVector, std::string xAxisName, std::string yAxisName, std::string type);
Double_t mc_p_n_poly4_slice_fit(Double_t *x, Double_t *par);
TH1D* shiftHistogramX(TH1D* originalHist, double shiftValue);
Double_t poly4(Double_t *x, Double_t *par);
TH1D* GetResidualHistogram(TH1D* hist, TF1* fit);
void adjustCanvas(TCanvas* canvas, double leftMargin = 0.15, double rightMargin = 0.05, double bottomMargin = 0.15, double topMargin = 0.10);
void customizeGraph(TGraphErrors *graph, int markerStyle, int markerColor, double markerSize, 
                    const char* graphTitle ="", const char* xAxisLabel="", const char* yAxisLabel="",
		    double TitleOffsetX = 1.4, double TitleOffsetY = 2, 
		    double LabelOffsetX = 0.01, double LabelOffsetY = 0.01);
void printParsedTitle(TH2D* hist,  const char* outputname="");
TH1D* sumHistogramsWithPolynomial(TH1D* h1, TH1D* h2, TF1* poly);
void SliceAndProjectHistogram_AroundMean(TH2D* hist2D, const std::vector<double>& xSlices, std::vector<TH1D*>& histVector, double xMean, std::string xAxisName, std::string yAxisName, std::string type);


void coin_stability_slices(){ // main
  // bit of a test for now. Will need to make this more sophisticated in the future. 

  gStyle->SetNumberContours(255); 
  gStyle->SetOptStat(0110);
  


  //TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_test.root"); // data
  TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_tight_2Dhistos.root"); // data

  //TFile *f2 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_test_mc_p.root"); //proton
  TFile *f2 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_tight_2Dhistos.root"); //proton

  //TFile *f3 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_test_mc_n.root"); //neutron
  TFile *f3 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_tight_2Dhistos.root"); //neutron



 

  // set location and name for output file 
  // TString outputfilename = "../output/" + configfileinput + ".root";
  TString outputfilename = "../output/coin_stability_slices_test.root";

  // Declare outfile
  // TFile *fout = new TFile( Form("../output/sbs%d.root",kine), "RECREATE" );
  //TFile *fout = new TFile("../output/sbs4_30p.root", "RECREATE" );
  TFile *fout = new TFile(outputfilename,"RECREATE");
  cout<<"writing to file: "<< outputfilename <<endl;


  // Load Histograms
  TH2D *hist_data_orig = (TH2D*)f1->Get("hcal_dx__hcal_sh_atime_diff");
  TH1D *hist_proton_orig_1d  = (TH1D*)f2->Get("hcal_dx_1d_allcuts");
  TH1D *hist_neutron_orig_1d  =  (TH1D*)f3->Get("hcal_dx_1d_allcuts");

  // Make clones of the histograms
  TH2D *hist_2D_data = (TH2D*)hist_data_orig->Clone("hist_2D_data");
  hist_p = (TH1D*)hist_proton_orig_1d->Clone("hist_p"); // global histogram
  hist_n = (TH1D*)hist_neutron_orig_1d->Clone("hist_n"); // global histogram

 

 // Set up the slices for the Y projections. 
  std::vector<double> xSlices = {-8,-5,-3,-2.5,-2,-1.5,-1,-0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3,5,8};


  double mean = 0;
  int nSlices = xSlices.size()-1;

  //// Vectors to store the sliced histograms in. 
  std::vector<TH1D*> hist_vector_data; 
  std::vector<TH1D*> hist_vector_p; 
  std::vector<TH1D*> hist_vector_n; 

 // use the function to make the slices and projections
  SliceAndProjectHistogram_slices(hist_2D_data, xSlices, hist_vector_data, "coin time","dx","data");

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

  double initialParameters[9]={1,1,0,0,1,1,1,1,1};
  // initialParameters = {0};
      
  TH1D *hist_data;
  /// Fit each slice 
  for (int sliceid = 0 ; sliceid < nSlices ; sliceid++)
    {
      /// set up the global histograms that the fit function is going to use
      // this uses clones of the histograms
      hist_data  =(TH1D*)hist_vector_data[sliceid]->Clone("hist_data");

      //// this uses the histograms themselves: you can't uncouple the fit from them easily
      //  hist_p  =hist_vector_p[sliceid];
      //  hist_n  =hist_vector_n[sliceid];
      // hist_data  =hist_vector_data[sliceid];

	 
      TF1 *Fit = new TF1(Form("overall_fit_%i_",sliceid),mc_p_n_poly4_slice_fit, xmin, xmax, 9);
      Fit->SetParameters(initialParameters); 
      // set parameter limits. SetParLimits -> (par#, min, max)
      Fit->SetParLimits(0, 0, 2000); // scale_p greater than 0
      Fit->SetParLimits(1, 0,2000); // scale_n greater than 0
      Fit->SetParLimits(2, -0.02,0.02); // shift_p less than +- 10cm
      Fit->SetParLimits(3, -0.02,0.02); // shift_n less than +- 10cm
      Fit->SetNpx(500);
      hist_data->GetXaxis() ->SetRangeUser(xmin, xmax);
      hist_data->Fit(Fit,"Q");

      // retrieve fit results 
      double scale_p  = Fit ->GetParameter(0);
      double scale_p_err = Fit ->GetParError(0);
      double scale_n  = Fit ->GetParameter(1);
      double scale_n_err = Fit ->GetParError(1);

      double shift_p= Fit ->GetParameter(2); 
      double shift_n = Fit ->GetParameter(3);
      double shift_p_err= Fit ->GetParError(2); 
      double shift_n_err = Fit ->GetParError(3);
	  
      double ChiSq= Fit->GetChisquare();
      double ndf = Fit->GetNDF();

      std::vector<double> poly_result;
      std::vector<double> poly_result_err;
      for (int i =0 ; i < 5; i++)
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
      shift_p_vector.push_back(shift_p);
      shift_n_vector.push_back(shift_n);
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


  TH1D *hist_temp_p;
  TH1D *hist_temp_n;

  // loop over the slices to create shifted and scaled versions of the mc histograms to plot
  for (int sliceid = 0 ; sliceid < nSlices ; sliceid++)
    {	
      /// shift the histograms. Function returns a clone. 
      hist_temp_p = shiftHistogramX(hist_p, shift_p_vector[sliceid] );
      hist_temp_n = shiftHistogramX(hist_n, shift_n_vector[sliceid] );

      /// scale the histograms 
      hist_temp_p -> Scale(scale_p_vector[sliceid]);
      hist_temp_n -> Scale(scale_n_vector[sliceid]);
	 
      // save histograms to the global vectors 
      hist_result_vector_p.push_back(hist_temp_p);
      hist_result_vector_n.push_back(hist_temp_n);	
    }// end loop over slices


  std::vector<TF1*> poly_fit_result;
  // loop over slices and make fits to plot the polynomial background result
  for (int sliceid = 0 ; sliceid < nSlices ; sliceid++)
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

      TF1 *fit = new TF1(Form("poly4_%i_",sliceid),poly4, xmin, xmax, 5);
      fit->SetParameters(poly_params_array); 
      fit->SetNpx(500);
      poly_fit_result.push_back(fit);
    }// end loop over slices


     //// Plot Rsf
     //// Make arrays that TGraphErrors can use 
  double x[nSlices];
  double y[nSlices];
  double x_err[nSlices];
  double y_err[nSlices];
  for (int sliceid = 0 ; sliceid < nSlices ; sliceid++)
    {
      x[sliceid] = xSlices[sliceid];
      y[sliceid] = Rsf_vector[sliceid];
      x_err[sliceid] = 0;
      y_err[sliceid] = Rsf_err_vector[sliceid];
    }    
  TGraphErrors *Rsf_graph = new TGraphErrors(nSlices, x, y, x_err, y_err);
  customizeGraph(Rsf_graph, 33, kBlue, 3,"","coin time","Rsf");
  void customizeGraph(TGraphErrors *graph, int markerStyle, int markerColor, double markerSize, 
  		      const char* graphTitle, const char* xAxisLabel, const char* yAxisLabel);

  //// canvas
  TCanvas *graphcanvas = new TCanvas("graphcanvas","graphcanvas",800,600);  
  adjustCanvas(graphcanvas);
  Rsf_graph->Draw("AP");
  graphcanvas->Update();
  graphcanvas->SaveAs(Form("%s/Rsf.pdf",output_dir.c_str() ));


  int nHist = hist_vector_data.size();
  int nCols = 4;
  int nRows = (nHist + nCols - 1) / nCols;


  std::vector<TH1D*> overall_fit_as_histogram;

  TCanvas* fits_canvas = new TCanvas("fits_canvas", "fits_canvas", 1000, 600);
  fits_canvas->Divide(nCols, nRows); 
  for (int i = 0; i < nHist; ++i) {
    fits_canvas->cd(i + 1);
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

    hist_result_vector_p[i] ->SetLineColor(kGreen);
    hist_result_vector_p[i]->Draw("same");
    hist_result_vector_n[i] ->SetLineColor(kMagenta);
    hist_result_vector_n[i]->Draw("same");
    sum_histo ->SetLineColor(kRed);
    sum_histo->Draw("same");
     
  }// end loop over slices
  fits_canvas->Update();
  fits_canvas->SaveAs(Form("%s/fitted_slices.pdf",output_dir.c_str() ));
  

  std::vector<TH1D*> hist_residual_vector;
 
  for (int sliceid = 0 ; sliceid < nSlices ; sliceid++)
    {
      TH1D* hist_residual = (TH1D*)hist_vector_data[sliceid]->Clone("hist_residual");  
      hist_residual  ->Add(overall_fit_as_histogram[sliceid], -1);
      hist_residual->GetXaxis() ->SetRangeUser(xmin, xmax);
      hist_residual_vector.push_back(hist_residual);
    }// end loop over slices

   
  TCanvas* resid_canvas = new TCanvas("resid_canvas", "resid_canvas", 1000, 600);
  resid_canvas->Divide(nCols, nRows);
  for (int i = 0; i < nHist; ++i) {
    resid_canvas->cd(i + 1);
    hist_residual_vector[i]->Draw();
  }
  resid_canvas->Update();
  resid_canvas->SaveAs(Form("%s/residuals.pdf",output_dir.c_str() ));

     
  // Create a canvas to draw the histogram
  TCanvas* cutcanvas = new TCanvas("cutcanvas", "visulaized cuts", 800, 600);

  // Draw the TH2D histogram
  hist_2D_data->Draw("COLZ");

  // Draw vertical lines at each x value in xSlices
  for (size_t i = 0; i < nSlices+1; ++i) {
    double x = xSlices[i];
    TLine* line = new TLine(x, hist_2D_data->GetYaxis()->GetXmin(), x, hist_2D_data->GetYaxis()->GetXmax());
    line->SetLineColor(kRed);  // Set line color 
    line->SetLineWidth(2);     // Set line width
    line->Draw("SAME");
  }

  // Save the canvas to a file or display it
  cutcanvas->SaveAs(Form("%s/visulaized_cuts.pdf",output_dir.c_str() ));

  // Create a projection of the TH2D histogram onto the x-axis
  TH1D* hist1D_data = hist_2D_data->ProjectionX("dataProjX");
  //TH1D* hist1D_p = hist_2D_proton->ProjectionX("protonProjX");
  //TH1D* hist1D_n = hist_2D_neutron->ProjectionX("neutronProjX");


  // Create a canvas to draw the histogram
  TCanvas* cutcanvas1D = new TCanvas("cutcanvas1D", "X Projection with Vertical Lines", 800, 600);
 
  // Draw the TH1D histogram
  cutcanvas1D->cd(1);
  //hist1D_data->GetXaxis() ->SetRangeUser(0, 0.1);
  hist1D_data->SetLineWidth(2);
  hist1D_data->Draw();
  // Draw vertical lines at each x value in xSlices
  for (size_t i = 0; i < nSlices+1; ++i) {
    double x = xSlices[i];
    TLine* line = new TLine(x, hist1D_data->GetMinimum(), x, hist1D_data->GetMaximum());
    line->SetLineColor(kRed);  // Set line color (e.g., red)
    line->SetLineWidth(2);     // Set line width
    line->Draw("SAME");
  }

  

  cutcanvas1D->SaveAs(Form("%s/visulaized_cuts_1D.pdf",output_dir.c_str() ));
 

  // TCanvas* testcanvas1  = new TCanvas("testcanvas1", "testcanvas1", 800, 600);
  // TH1D* sumHistoTest = sumHistogramsWithPolynomial(hist_result_vector_p[0],hist_result_vector_n[0] , poly_fit_result[0]);
  // hist_vector_data[0]->Draw();
  // sumHistoTest ->Draw("same");

  
 
  //// extract the histogram title and print it 
  printParsedTitle(hist_2D_data,"data");

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


void SliceAndProjectHistogram_slices(TH2D* hist2D, const std::vector<double>& xSlices, std::vector<TH1D*>& histVector, std::string xAxisName, std::string yAxisName, std::string type) {
    // Clear the vector to ensure it's empty before filling
    histVector.clear();

    // loop over the slices 
    for (size_t i = 0; i < xSlices.size()-1 ; ++i) { 
      double xMin = xSlices[i];
      int binMin = hist2D ->GetXaxis()->FindBin(xMin);

      double xMax = xSlices[i+1];
      int binMax = hist2D ->GetXaxis()->FindBin(xMax);

      // Define the name and title for the TH1D histogram
      // Format the histogram name to display only two decimal places
      std::ostringstream stream;
      stream << std::fixed << std::setprecision(4) << "_"<<xMin << "_to_" << xMax;
      std::string histName = yAxisName+"__"+xAxisName + stream.str() + "_" +type;
      TH1D *projY = hist2D->ProjectionY(histName.c_str(), binMin, binMax);
      histVector.push_back(projY); // 

    }// end loop over slices
}// end SliceAndProjectHistogram_slices



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
                    const char* graphTitle ="", const char* xAxisLabel="", const char* yAxisLabel="",
		    double TitleOffsetX = 1.4, double TitleOffsetY = 2, 
		    double LabelOffsetX = 0.01, double LabelOffsetY = 0.01) {
  // Set marker style, color, and size
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerColor(markerColor);
  graph->SetMarkerSize(markerSize);

  // Set graph title and axis labels
  graph->SetTitle(graphTitle);
  graph->GetXaxis()->SetTitle(xAxisLabel);
  graph->GetYaxis()->SetTitle(yAxisLabel);

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
void printParsedTitle(TH2D* hist,  const char* outputname="") {
    // Get the histogram title
    std::string title = hist->GetTitle();

    cout<<title<<endl;
  
    // // Find the position of the first '{' character
  //   size_t pos = title.find('{');
    
  //   // Extract the y_axis:x_axis part
  // std::string axes = title.substr(0, pos);
  // //cout<<"axis = "<<axes<<endl;
    
  // // Extract the cuts part and remove '{' and '}'
  // std::string cuts = title.substr(pos + 1, title.size() - pos - 2);
  // //cout<<"cuts = "<<cuts<<endl;
   
  // //cout<<"broken up cuts"<<endl;

  // // Split the cuts into individual cut expressions
  // std::vector<std::string> cutList;
  // std::stringstream ss(cuts);
  // std::string cut;
  // while (std::getline(ss, cut, '&')) {
  //   // Remove leading and trailing whitespace
  //   cut.erase(0, cut.find_first_not_of(" \t"));
  //   cut.erase(cut.find_last_not_of(" \t") + 1);
        
  //   // Ensure the cut is not empty before processing
  //   if (!cut.empty() && cut.front() == '&') {
  //     cut.erase(cut.begin());
  //   }
        
  //   if (!cut.empty()) {
  //     cutList.push_back(cut);
  //     // cout<<cut<<endl;
  //   }
  // }

    
  // // Create a new canvas
  // TCanvas* cuts_canvas = new TCanvas("cuts_canvas", "Parsed Histogram Title", 1000, 600);
    
  // // Create a TLatex object to draw the text
  // TLatex latex;
  // latex.SetTextSize(0.03);  // Adjust text size
  // latex.SetTextAlign(13);   // Align text to top left
    
  // // Draw the axes part
  // latex.DrawLatex(0.1, 0.9, axes.c_str());
    
  // // Draw each cut expression on a new line
  // double yPos = 0.8;  // Start position for the first cut
  // for (const auto& cut : cutList) {
  //   latex.DrawLatex(0.1, yPos, cut.c_str());
  //   yPos -= 0.03;  // Move down for the next cut
  // }


    
  //   // Update the canvas
  //  cuts_canvas->Update();
    
  //   // Optionally, save the canvas as an image
  //   cuts_canvas->SaveAs(Form("%s/global_cuts_%s.pdf", output_dir.c_str(),outputname));
}
