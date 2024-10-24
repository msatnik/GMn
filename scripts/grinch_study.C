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


//// Vectors to store the sliced histograms in. Declaring them outside so all the functions can see them. 
std::vector<TH1D*> hist_vector_data; 
std::vector<TH1D*> hist_vector_p; 
std::vector<TH1D*> hist_vector_n; 

std::vector<TH1D*> hist_result_vector_p; 
std::vector<TH1D*> hist_result_vector_n; 

std::vector<TH1D*> hist_residual_vector;


std::vector<TH1D*> hist_vector_data_sum;


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
double correction_factor = 0.01;


//// ranges for slices. These are what we want to change.
////****** dy *************************************************************** 
std::string xMin_xMax_string_dy = "(-1,-0.5),(-0.5,0),(0,0.5),(0.5,1)";
////*****nsigdy***************************************************************
std::string xMin_xMax_string_nsigdy = "(-4.1,-2.1),(-3.1,-1.1),(-2,0.1),(-1,1),(0,2),(1.1,3),(2.1,4)";
///********vz**********************************************************************
std::string xMin_xMax_string_vz = "(-0.07,0.07),(-0.065,0.065),(-0.06,0.06),(-0.05,0.05),(-0.04,0.04),(-0.03,0.03)";
//// ******ps_e**********************************************************************
std::string xMin_xMax_string_ps =  "(0.1,0.2),(0.12,0.22),(0.14,0.24),(0.16,0.26),(0.18,0.28),(0.2,0.3),(0.22,0.32),(0.24,0.34),(0.26,0.36),(0.28,0.38)";
//std::string xMin_xMax_string_ps =  "(0.1,2),(0.12,2),(0.14,2),(0.16,2),(0.18,2),(0.2,2),(0.22,2),(0.24,2),(0.26,2),(0.28,2)";
//std::string xMin_xMax_string_ps ="(0.2,2)";
////******* hcal energy **************************************************************
//std::string xMin_xMax_string_hcal_e = "(0.015,0.025), (0.020,0.030), (0.025,0.035),(0.030,0.040)";
std::string xMin_xMax_string_hcal_e = "(0.015,1), (0.020,1), (0.025,1),(0.030,1)";
//// fiducial x
//// fiducial y
//// ******w2****************************************************************
//std::string xMin_xMax_string_w2 = "(0.4,0.6),(0.5,0.7),(0.6,0.8),(0.7,0.9),(0.8,1.0),(0.9,1.1),(1.0,1.2),(1.1,1.3) ";
//std::string xMin_xMax_string_w2 = "(0.4,1.2),(0.5,1.1),(0.6,1.0),(0.7,0.9)";
//std::string xMin_xMax_string_w2 = "(0.78,0.98),(0.68,1.08),(0.58,1.18)";
std::string xMin_xMax_string_w2 = "(0.4,0.5),(0.5,0.6),(0.6,0.7),(0.7,0.8),(0.8,0.9),(0.9,1.0),(1.0,1.1),(1.1,1.2),(1.2,1.3),(1.3,1.4),(1.4,1.5)";
// e_over_p
//// *** x_expected*******************************************
//std::string xMin_xMax_string_x_exp = "(-1.5,-1),(-1.25,-0.75),(-1,-0.5),(-0.75,-0.25),(-0.5,0),(-0.25,0.25),(0,0.5),(0.25,0.75),(0.5,1),(0.75,1.25)";
std::string xMin_xMax_string_x_exp = "(-0.8,0.0),(-0.9,0.1),(-1.0,0.2),(-1.1,0.3 )(-1.2,0.4) ,(-1.3,0.5),(-1.4,0.6)";
///****** nsigx_fid *******************************************
//std::string xMin_xMax_string_nsigx_fid = "(0,1),(0.5,1.5),(1,2),(1.5,2.5),(2,3),(2.5,3.5),(3.5,4.5),(4,5),(5,6),(5.5,6.5),(6,7),(6.5,7.5)";
//std::string xMin_xMax_string_nsigx_fid = "(0,7),(0.5,7),(1,7),(1.5,7),(2,7),(2.5,7),(3,7),(3.5,7),(4,7),(4.5,7),(5,7),(5.5,7),(6,7)";
std::string xMin_xMax_string_nsigx_fid = "(1.5,7),(2,7),(2.5,7),(3,7),(3.5,7)";
///******* y expected*************************************************************
//std::string xMin_xMax_string_y_exp = "(-1,-0.5),(-0.75,-0.25),(-0.5,0),(-0.25,0.25),(0,0.5),(0.25,0.75),(0.5,1)";
std::string xMin_xMax_string_y_exp = "(-0.35,0.4)";

//********* nsigy_fid ***************************************************************
std::string xMin_xMax_string_nsigy_fid = "(0,3),(0.5,3),(1,3),(1.5,3),(2,3)";


/// GRINCH study

std::string xMin_xMax_string_grinch_clus_size ="(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,10),(10,11),(11,12),(12,13),(13,14),(14,50)";
//std::string xMin_xMax_string_grinch_clus_adc = "(0,5),(5,10),(10,15),(15,20),(20,25),(25,30),(30,35),(35,40),(40,45),(45,50),(50,55),(55,60),(60,65),(65,70),(70,75),(75,80),(80,85),(85,500)";
std::string xMin_xMax_string_grinch_clus_adc = "(0,5),(5,10),(10,15),(15,20),(20,25),(25,30),(30,35),(35,40),(40,45),(45,50),(50,55),(55,60),(60,65),(65,70),(70,75),(75,80),(80,85),(85,90),(90,95),(95,100),(100,105),(105,110),(110,115),(115,120),(120,125),(125,130),(130,135),(135,140),(140,145),(145,150),(150,155),(155,160),(160,165),(165,170),(170,175),(175,180),(180,185),(185,190),(190,195),(195,200),(200,500)";


//std::string xMin_xMax_string_grinch_clus_size = "(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,10),(10,11),(11,12),(12,13),(13,14),(14,15),(15,16),(16,17),(17,18),(18,19),(19,20),(20,21),(21,22),(22,23),(23,24),(24,25),(25,26),(26,27),(27,28),(28,29),(29,30),(30,31),(31,32),(32,33),(33,34),(34,35),(35,36),(36,37),(37,38),(38,39),(39,40),(40,41),(41,42),(42,43),(43,44),(44,45),(45,46),(46,47),(47,48),(48,49),(49,50)";

//std::string xMin_xMax_string_grinch_clus_adc = "(0,5),(5,10),(10,15),(15,20),(20,25),(25,30),(30,35),(35,40),(40,45),(45,50),(50,55),(55,60),(60,65),(65,70),(70,75),(75,80),(80,85),(85,90),(90,95),(95,100),(100,105),(105,110),(110,115),(115,120),(120,125),(125,130),(130,135),(135,140),(140,145),(145,150),(150,155),(155,160),(160,165),(165,170),(170,175),(175,180),(180,185),(185,190),(190,195),(195,200),(200,205),(205,210),(210,215),(215,220),(220,225),(225,230),(230,235),(235,240),(240,245),(245,250),(250,255),(255,260),(260,265),(265,270),(270,275),(275,280),(280,285),(285,290),(290,295),(295,300),(300,305),(305,310),(310,315),(315,320),(320,325),(325,330),(330,335),(335,340),(340,345),(345,350),(350,355),(355,360),(360,365),(365,370),(370,375),(375,380),(380,385),(385,390),(390,395),(395,400),(400,405),(405,410),(410,415),(415,420),(420,425),(425,430),(430,435),(435,440),(440,445),(445,450),(450,455),(455,460),(460,465),(465,470),(470,475),(475,480),(480,485),(485,490),(490,495),(495,500)";



// SBS4 30p 
TString DataFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_2Dhistos.root";
TString ProtonFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_2Dhistos.root";
TString NeutronFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_2Dhistos.root";
double xmin_sbs4_30p = -2.1; // 
double xmax_sbs4_30p = 1.4;//
double correction_factor_sbs4_30p = 0.01;


// SBS4 50p 
TString DataFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_2Dhistos_sept19.root";
TString ProtonFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deep_2Dhistos.root";
TString NeutronFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deen_2Dhistos.root";
double xmin_sbs4_50p = -2.5; // 
double xmax_sbs4_50p = 1.4;//
double correction_factor_sbs4_50p= 0.01;

// SBS8 70p
TString DataFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_2Dhistos_ps_sept25.root";// sept25 has cuts on grinch clus >1 
TString ProtonFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deep_2Dhistos.root";
TString NeutronFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deen_2Dhistos.root";
double xmin_sbs8_70p = 0; // 
double xmax_sbs8_70p = 3;//
double correction_factor_sbs8_70p=  0.0055;

// SBS9 70p
TString DataFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_2Dhistos_ps_sept25.root"; // sept25 has cuts on grinch clus >1
TString ProtonFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_deep_2Dhistos_sept24_z.root";
TString NeutronFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_deen_2Dhistos_sept24_z.root";
double xmin_sbs9_70p = 0; // need to adjust range
double xmax_sbs9_70p = 3;//
double correction_factor_sbs9_70p= 0.05;




// Functions 


void grinch_study(TString KineString="sbs9_70p", TString VarString = "grinch_clus_adc"){ // main
  // bit of a test for now. Will need to make this more sophisticated in the future. 

  gStyle->SetNumberContours(255); 
  gStyle->SetOptStat(0110);

  //std::vector<int> color_vector = {kRed,kBlue,kGreen,kMagenta,kCyan,kYellow};
  std::vector<int> color_vector = {kRed,kRed};

  Utility utilityHandler; // class that gives us access to various functions to use between programs.

  
  // Taking user input on the kinematic and what variable we want to study. Going to try keeping it all here instead of config files.

  //double correction = 0.1;//  sbs9 with no electron cuts 
  //double correction = 0.019;//  sbs9 with electron cuts no track match

  //double correction = 0.075;//  sbs9 with no electron cuts track match
  //double correction = 0.0187;//  sbs9 with electron cuts track match

 
  //double correction = 0;//  sbs8 with electron cuts 
  //double correction = 0.0105;//  sbs8 with no electron cuts

  // double correction = 0.01005;//  sbs8 with no electron cuts and track match


  // SBS9 with more preshower
  // double correction = 0.025;// electron cuts NO track match
  //double correction = 0.02;// NO electron cuts NO track match

  //double correction = 0.0;// zero!

  //double correction = 0.0235;// electron cuts track match

  double correction = 0.02445;// electron cuts track match  _trackmatch
  //double correction = 0;

  //double correction = 0.1103;// no electron cuts track match _nocuts_trackmatch
  
  if (VarString == "grinch_clus_size" ){
    HistogramName = "ps_e__grinch_clus_size_trackmatch";
    AxisTitle ="grinch_clus_size";
    xMin_xMax_string = xMin_xMax_string_grinch_clus_size;
  }else if (VarString == "grinch_clus_adc" ){
    HistogramName = "ps_e__grinch_clus_adc_trackmatch";
    AxisTitle ="grinch_clus_adc";
    xMin_xMax_string = xMin_xMax_string_grinch_clus_adc;
  }  else {
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
  }else if (KineString == "sbs9_70p"){
    DataFileString = DataFileString_sbs9_70p;
    ProtonFileString =ProtonFileString_sbs9_70p;
    NeutronFileString = NeutronFileString_sbs9_70p;
    xmin = xmin_sbs9_70p;
    xmax = xmax_sbs9_70p;
  }else {
    std::cout<<"Error with kinematic setting"<<std::endl;
    return;
  }


  cout<<"What we're running: "<<KineString<< ", "<<VarString<<endl;
  cout<<DataFileString<<endl;
  //cout<<ProtonFileString<<endl;
  // cout<<NeutronFileString<<endl;
  cout<<xmin<< ", "<<xmax<<endl;
  cout<<HistogramName<<endl;
  cout<<AxisTitle<<endl;
  cout<<xMin_xMax_string<<endl;
 
      
      
  
  TFile *f1 = TFile::Open(DataFileString); // data
  TFile *f2 = TFile::Open(ProtonFileString); //proton
  TFile *f3 = TFile::Open(NeutronFileString); //neutron



  //// set location and name for root output file 
  TString outputfilelocation="../output/grinch_output/"+KineString+"/"+VarString;
  TString outputfilename = outputfilelocation +"/"+ KineString+ "_" +VarString+".root";
  
  // TFile *fout = new TFile("../output/stability/testoutput.root","RECREATE");

  TFile *fout = new TFile(outputfilename,"RECREATE");
  cout<<"writing to file: "<< outputfilename <<endl;


  // Load Histograms
  TH2D *hist_data_orig = (TH2D*)f1->Get(HistogramName);
  TH2D *hist_proton_orig = (TH2D*)f2->Get("hcal_dx__ps_e");
  TH2D *hist_neutron_orig = (TH2D*)f3->Get("hcal_dx__ps_e");
 

  // Make clones of the histograms
  TH2D *hist_2D_data = (TH2D*)hist_data_orig->Clone("hist_2D_data");
  TH2D *hist_2D_proton = (TH2D*)hist_proton_orig->Clone("hist_2D_proton");
  TH2D *hist_2D_neutron = (TH2D*)hist_neutron_orig->Clone("hist_2D_neutron");

  TH1D* hist1D_data = hist_2D_data->ProjectionX("dataProjX");
  TH1D* hist1D_proton = hist_2D_proton->ProjectionX("protonProjX");
  TH1D* hist1D_neutron = hist_2D_neutron->ProjectionX("neutronProjX");




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
 


  //// Vectors to store the sliced histograms in. 
  std::vector<TH1D*> hist_vector_data; 
  std::vector<TH1D*> hist_vector_p; 
  std::vector<TH1D*> hist_vector_n;

  
  utilityHandler.SliceAndProjectHistogram_xMinxMax_inclusiveMin(hist_2D_data, xMinSlices,xMaxSlices, hist_vector_data, AxisTitle,"dx","data");

 

  std::vector<double> nPions_vector;
  std::vector<double> nElectrons_vector;
  std::vector<double> electrons_to_pions_ratio;
  std::vector<double> pions_to_electrons_ratio;
  std::vector<double> electron_detection_eff;
  std::vector<double> electron_detection_eff_corrected;
  std::vector<double> pion_rejection_eff;
  
  std::vector<double> pion_rejection_eff_test;
  std::vector<double> pion_rejection_eff_corrected;

  std::vector<double> nPions_vector_corrected;
  std::vector<double> nElectrons_vector_corrected;

  double ps_cut = 0.2;

  // Project the 2D histogram onto the Y-axis over the entire X-axis range
  TH1D* full_Y_projection = hist_2D_data->ProjectionY("full_Y_projection", 0, -1);
  // Estimate that pions are anything with ps_e =< ps_cut
  /// Estimate that electrons are anything ps_e >ps_cut
  // Get the number of bins in the histogram
  int nBins = full_Y_projection->GetNbinsX();
  // Get the bin corresponding to the lower edge of ps_cut
  int binLowEdge = full_Y_projection->FindBin(ps_cut); // This bin contains the value ps_cut
  // Get the bin corresponding to the upper edge of the histogram
  int lastBin = full_Y_projection->GetNbinsX();
  // Integral from [0, ps_cut] (inclusive)
  double total_pions = full_Y_projection->Integral(1, binLowEdge); // from the first bin to the bin containing ps_cut
  // Integral from (ps_cut, max] (exclusive of the ps_cut bin)
  double total_electrons = full_Y_projection->Integral(binLowEdge + 1, lastBin); // from the next bin to the last bin
  TH1D* full_X_projection = hist_2D_data->ProjectionX("full_X_projection", 0, -1);


  // Project the 2D histogram onto the X-axis over the entire Y-axis range
  TH1D* full_X_projection_proton = hist_2D_proton->ProjectionX("full_X_projection_proton", 0, -1);
  TH1D* full_X_projection_neutron = hist_2D_neutron->ProjectionX("full_X_projection_neutron", 0, -1);
  TH1D* full_X_projection_sim = (TH1D*)full_X_projection_proton->Clone("full_X_projection_sim");
  full_X_projection_proton->Add(full_X_projection_neutron,1);
  // Estimate that pions are anything with ps_e =< ps_cut
  /// Estimate that electrons are anything ps_e >ps_cut
  // Get the number of bins in the histogram
  int nBins_proton = full_X_projection_proton->GetNbinsX();
  // Get the bin corresponding to the lower edge of ps_cut
  int binLowEdge_proton = full_X_projection_proton->FindBin(ps_cut); // This bin contains the value ps_cut
  // Get the bin corresponding to the upper edge of the histogram
  int lastBin_proton = full_X_projection_proton->GetNbinsX();
  // Integral from [0, ps_cut] (inclusive)
  double total_pions_proton = full_X_projection_proton->Integral(1, binLowEdge); // from the first bin to the bin containing ps_cut
  // Integral from (ps_cut, max] (exclusive of the ps_cut bin)
  double total_electrons_proton = full_X_projection_proton->Integral(binLowEdge + 1, lastBin); // from the next bin to the last bin
  TH1D* full_Y_projection_proton = hist_2D_proton->ProjectionX("full_Y_projection_proton", 0, -1);

  double pion_to_electron_ratio_from_proton = total_pions_proton/total_electrons_proton;
  cout<<endl<<"Testing finding the ''pion to electron ratio'' from proton simulation"<<endl;
  cout<<"ratio = "<<total_pions_proton<<" / "<<total_electrons_proton<<" = "<<pion_to_electron_ratio_from_proton <<endl<<endl;

  
  TH1D *hist_data;

  //double correction = 0.0099611;

  // Elastic cuts 
  //double correction = 0.05;
  //double correction = 0.0062;

  //double correction = 0.010;
  //double correction = 0.05;


  
  //double correction = 0.0476; // sbs9 with electron cuts
  // double correction = 0.079; // sbs9 with no electron cuts

   double total_electrons_corrected = total_electrons + correction*total_electrons;
   
  double total_pions_corrected = total_pions - correction*total_electrons;

  cout<<endl;
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {
      // this uses clones of the histograms
      hist_data  =(TH1D*)hist_vector_data[sliceid]->Clone("hist_data");
      // Estimating that pions are anything ps_e < ps_cut, and electrons are anything over
      
      double nPions = hist_data->Integral(1,binLowEdge); // from the first bin to the bin containing ps_cut
      double nElectrons= hist_data->Integral(binLowEdge+1, lastBin);// from the next bin to the last bin
      double nPions_corrected = nPions - correction*nElectrons;
      double nElectrons_corrected = nElectrons + correction*nElectrons;
      // cout<<"slice ID: "<<sliceid<<endl;
      // cout<<"nPions in slice: "<<nPions<< endl;
      // cout<<"nElectrons in slice: "<<nElectrons<<endl;
      // cout<<"nPions - correction*nElectrons: "<<nPions_corrected<<endl;
      // cout<<"nElectrons + correction * nElectrons: "<<nElectrons_corrected<<endl<<endl;
      nPions_vector.push_back(nPions);
      nPions_vector_corrected.push_back(nPions_corrected);
      nElectrons_vector.push_back(nElectrons);
      nElectrons_vector_corrected.push_back(nElectrons_corrected);
    }// end loop over histogram slices

  std::vector<double> detected_electrons_sum;
  std::vector<double> detected_electrons_sum_corrected;
  std::vector<double> detected_pions_sum;
  std::vector<double> detected_pions_sum_corrected;
 
  std::vector<double> rejected_electrons_sum;
  std::vector<double> rejected_pions_sum;
  std::vector<double> rejected_pions_sum_corrected;
  
  TH1D* hist_vector_data_sum_temp;

  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {
      detected_electrons_sum.push_back(0);
      detected_electrons_sum_corrected.push_back(0);
      detected_pions_sum.push_back(0);
      detected_pions_sum_corrected.push_back(0);
      rejected_electrons_sum.push_back(0);
      rejected_pions_sum.push_back(0);
      rejected_pions_sum_corrected.push_back(0);
      pion_rejection_eff.push_back(0);
      pion_rejection_eff_test.push_back(0);
      pion_rejection_eff_corrected.push_back(0);
      electron_detection_eff.push_back(0);
      electron_detection_eff_corrected.push_back(0);
      if(nPions_vector[sliceid]!=0&&nElectrons_vector[sliceid]!=0){
	electrons_to_pions_ratio.push_back(nElectrons_vector[sliceid]/nPions_vector[sliceid]);
	pions_to_electrons_ratio.push_back(nPions_vector[sliceid]/nElectrons_vector[sliceid]);
      }
      //hist_vector_data_sum_temp->Reset();
      for (int inner_loop_cnt = sliceid; inner_loop_cnt < xMinSlices.size() ; inner_loop_cnt++)
	{// sum up the number of "electrons" and "pions" in the slices that would be accepted if we made the cut 
	  detected_electrons_sum[sliceid] += nElectrons_vector[inner_loop_cnt];
	  detected_electrons_sum_corrected[sliceid] += nElectrons_vector_corrected[inner_loop_cnt];
	  detected_pions_sum[sliceid] += nPions_vector[inner_loop_cnt];
	  detected_pions_sum_corrected[sliceid] += nPions_vector_corrected[inner_loop_cnt];
	  // hist_vector_data_sum_temp += hist_vector_data[inner_loop_cnt]
	}//end inner loop
      for (int inner_loop_cnt = sliceid-1; inner_loop_cnt >= 0 ; inner_loop_cnt--) // trying counting backwards
	{// sum up the number of "electrons" and "pions" in the slices that would be REJECTED
	  //cout<<"inner loop backwards count"<<inner_loop_cnt<<endl;
	  rejected_electrons_sum[sliceid] += nElectrons_vector[inner_loop_cnt];
	  rejected_pions_sum[sliceid] += nPions_vector[inner_loop_cnt];
	  rejected_pions_sum_corrected[sliceid] += nPions_vector_corrected[inner_loop_cnt];
	  // hist_vector_data_sum_temp += hist_vector_data[inner_loop_cnt]
	}//end inner loop backwards
      //cout<<endl;
      pion_rejection_eff_test[sliceid] = rejected_pions_sum[sliceid]/total_pions;
      // pion_rejection_eff_corrected[sliceid] = rejected_pions_sum_corrected[sliceid]/total_pions_corrected;
      pion_rejection_eff[sliceid] =(total_pions-detected_pions_sum[sliceid])/total_pions;
      pion_rejection_eff_corrected[sliceid] =(total_pions_corrected-detected_pions_sum_corrected[sliceid])/total_pions_corrected;
      electron_detection_eff[sliceid] = detected_electrons_sum[sliceid]/total_electrons;
      electron_detection_eff_corrected[sliceid] = detected_electrons_sum_corrected[sliceid]/total_electrons_corrected;
      //hist_vector_data_sum.pushback(hist_vector_data_sum_temp);
    }// end loop over histogram slices 


  cout<<endl;
  cout<<"Total Electrons in Preshower: "<<total_electrons<<endl;
  cout<<"Total Pions in Preshower: "<<total_pions<<endl;
  cout<<"Total Pions - correction * Total Electrons: "<<total_pions_corrected<<endl;
  cout<<"Total Electrons + correction * Total Electrons: "<<total_electrons_corrected<<endl;
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {
      // cout<<"Slice ID: "<<sliceid<<endl;
      // cout<<"GRINCH slice range: ("<<xMinSlices[sliceid]<<","<<xMaxSlices[sliceid]<<")"<<endl;
      // cout<<"Electrons in Slice: "<<nElectrons_vector[sliceid]<<endl;
      // cout<<"Pions in Slice: "<<nPions_vector[sliceid]<<endl;
      // cout<<"Pions - correction * Electrons: "<<nPions_vector_corrected[sliceid]<<endl;
      //  cout<<"Electrons + correction * Electrons: "<<nElectrons_vector_corrected[sliceid]<<endl;
      // cout<<"Electron to Pion ratio in slice: "<<electrons_to_pions_ratio[sliceid]<<endl;
      // cout<<"Pion to Electron ratio in slice: "<<pions_to_electrons_ratio[sliceid]<<endl;
      // cout<<"------------------------------------------------------------------"<<endl;
      // cout<<"Electron Detection Eff using "<<xMinSlices[sliceid]<<" as cut: "<<electron_detection_eff[sliceid]<<endl;
      //  cout<<"Electron Detection Eff Corrected using "<<xMinSlices[sliceid]<<" as cut: "<<electron_detection_eff_corrected[sliceid]<<endl;
      // cout<<"Pion Rejection Eff using "<<xMinSlices[sliceid]<<" as cut: "<<pion_rejection_eff[sliceid]<<endl;
      //  cout<<"Pion Rejection Eff TEST using "<<xMinSlices[sliceid]<<" as cut: "<<pion_rejection_eff_test[sliceid]<<endl;
      //   cout<<"Pion Rejection Eff Corrected using "<<xMinSlices[sliceid]<<" as cut: "<<pion_rejection_eff_corrected[sliceid]<<endl;
      // cout<<"-----------------------------------------------------------------"<<endl;
    }// end loop over histogram slices
  cout<<endl;

  int max = xMinSlices.size()-1;

  cout<<"UNCORRECTED"<<endl;
  cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
  cout<<"Slice ID 0"<<endl;
  cout<<"Electron Detection Eff  using "<<xMinSlices[0]<<" as cut: "<<electron_detection_eff[0]<<endl;
  cout<<"Pion Rejection Eff using "<<xMinSlices[0]<<" as cut: "<<pion_rejection_eff[0]<<endl;
  cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
  cout<<"Slice ID "<<max<<endl;
  cout<<"Electron Detection Eff using "<<xMinSlices[max]<<" as cut: "<<electron_detection_eff[max]<<endl;
  cout<<"Pion Rejection Eff using "<<xMinSlices[max]<<" as cut: "<<pion_rejection_eff[max]<<endl;
  cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;

  cout<<endl<<"CORRECTED"<<endl;
  cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
  cout<<"Slice ID 0"<<endl;
  cout<<"Electron Detection Eff Corrected using "<<xMinSlices[0]<<" as cut: "<<electron_detection_eff_corrected[0]<<endl;
  cout<<"Pion Rejection Eff Corrected using "<<xMinSlices[0]<<" as cut: "<<pion_rejection_eff_corrected[0]<<endl;
  cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
  cout<<"Slice ID "<<max<<endl;
  cout<<"Electron Detection Eff Corrected using "<<xMinSlices[max]<<" as cut: "<<electron_detection_eff_corrected[max]<<endl;
  cout<<"Pion Rejection Eff Corrected using "<<xMinSlices[max]<<" as cut: "<<pion_rejection_eff_corrected[max]<<endl;
  cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;

  

   
   

  TCanvas* ps_gr_2d_canvas = new TCanvas("ps_gr_2d_canvas", "ps_gr_2d_canvas", 1000, 600);
  hist_2D_data->Draw();


  TCanvas* ps_gr_canvas = new TCanvas("ps_gr_canvas", "ps_gr_canvas", 1000, 600);
  ps_gr_canvas->Divide(2,1);
  ps_gr_canvas->cd(1);
  full_Y_projection->Draw();
  ps_gr_canvas->cd(2);
  full_X_projection->Draw();

   TCanvas* dx_ps_canvas = new TCanvas("dx_ps_canvas", "dx_ps_canvas", 1000, 600);
  dx_ps_canvas->Divide(2,1);
  dx_ps_canvas->cd(1);
  full_Y_projection_proton->Draw();
  dx_ps_canvas->cd(2);
  full_X_projection_proton->Draw();
  
  
  
  int nHist = hist_vector_data.size();
  int nCols = 4;
  int nRows = (nHist + nCols - 1) / nCols;

 
 
  TCanvas* slice_canvas = new TCanvas("slice_canvas", "slice_canvas", 1000, 600);
  if (hist_vector_data.size()!=1)
    {
      slice_canvas->Divide(nCols, nRows);
    }
  for (int i = 0; i < nHist; ++i) {
    slice_canvas->cd(i + 1);
    hist_vector_data[i]->GetXaxis() ->SetRangeUser(xmin, xmax);
    hist_vector_data[i]->Draw();
  }
  
  
  /// graph the electron detection eff
  double x[xMinSlices.size()];
  double y[xMinSlices.size()];
  double x_err[xMinSlices.size()];
  double y_err[xMinSlices.size()];
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {
      double xCenter = xMinSlices[sliceid];
     
      x[sliceid] = xCenter;
      y[sliceid] = electron_detection_eff[sliceid];
      x_err[sliceid] = 0;
      y_err[sliceid] = 0;

    }    
  TGraphErrors *elec_detec_eff_graph = new TGraphErrors(xMinSlices.size(), x, y, x_err, y_err);
  utilityHandler.customizeGraphMore(elec_detec_eff_graph, 33, kBlue, 3,"",AxisTitle,"electron detection efficiency");


  

  //// canvas
  TCanvas *graphcanvas = new TCanvas("graphcanvas","graphcanvas",800,600);
  graphcanvas->SetGrid();
  utilityHandler.adjustCanvas(graphcanvas);
  elec_detec_eff_graph->Draw("AP");
  graphcanvas->Update();
  graphcanvas->SaveAs(Form("%s/elec_detec_eff.pdf",outputfilelocation.Data() ) );


  /// graph the pion rejection eff  
  double x_pion_rej_eff[xMinSlices.size()];
  double y_pion_rej_eff[xMinSlices.size()];
  double x_pion_rej_eff_err[xMinSlices.size()];
  double y_pion_rej_eff_err[xMinSlices.size()];
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {
      double xCenter = xMinSlices[sliceid];
     
      x_pion_rej_eff[sliceid] = xCenter;
      y_pion_rej_eff[sliceid] = pion_rejection_eff[sliceid];
      x_pion_rej_eff_err[sliceid] = 0;
      y_pion_rej_eff_err[sliceid] = 0;

    }    
  TGraphErrors *pion_rej_eff_graph = new TGraphErrors(xMinSlices.size(), x_pion_rej_eff, y_pion_rej_eff, x_pion_rej_eff_err, y_pion_rej_eff_err);
  utilityHandler.customizeGraphMore(pion_rej_eff_graph, 33, kRed, 3,"",AxisTitle,"pion rejection eff");


  /// graph the pion rejection eff  
  double x_pion_rej_eff_corrected[xMinSlices.size()];
  double y_pion_rej_eff_corrected[xMinSlices.size()];
  double x_pion_rej_eff_err_corrected[xMinSlices.size()];
  double y_pion_rej_eff_err_corrected[xMinSlices.size()];
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {
      double xCenter = xMinSlices[sliceid];
     
      x_pion_rej_eff_corrected[sliceid] = xCenter;
      y_pion_rej_eff_corrected[sliceid] = pion_rejection_eff_corrected[sliceid];
      x_pion_rej_eff_err_corrected[sliceid] = 0;
      y_pion_rej_eff_err_corrected[sliceid] = 0;
    }
  TGraphErrors *pion_rej_eff_graph_corrected = new TGraphErrors(xMinSlices.size(), x_pion_rej_eff_corrected, y_pion_rej_eff_corrected, x_pion_rej_eff_err_corrected, y_pion_rej_eff_err_corrected);
  utilityHandler.customizeGraph(pion_rej_eff_graph_corrected, 33, kMagenta, 3);


  /// graph the elec detec eff  
  double x_elec_detec_eff_corrected[xMinSlices.size()];
  double y_elec_detec_eff_corrected[xMinSlices.size()];
  double x_elec_detec_eff_err_corrected[xMinSlices.size()];
  double y_elec_detec_eff_err_corrected[xMinSlices.size()];
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {
      double xCenter = xMinSlices[sliceid];
     
      x_elec_detec_eff_corrected[sliceid] = xCenter;
      y_elec_detec_eff_corrected[sliceid] = electron_detection_eff_corrected[sliceid];
      x_elec_detec_eff_err_corrected[sliceid] = 0;
      y_elec_detec_eff_err_corrected[sliceid] = 0;
    }
  TGraphErrors *elec_detec_eff_graph_corrected = new TGraphErrors(xMinSlices.size(), x_elec_detec_eff_corrected, y_elec_detec_eff_corrected, x_elec_detec_eff_err_corrected, y_elec_detec_eff_err_corrected);
  utilityHandler.customizeGraphMore(elec_detec_eff_graph_corrected, 33, kGreen, 3,"",AxisTitle,"");


  
  //// canvas
  TCanvas *graphcanvas_uncorrected = new TCanvas("graphcanvas_uncorrected","graphcanvas_uncorrected",800,600);
  graphcanvas_uncorrected->SetGrid();
  utilityHandler.adjustCanvas(graphcanvas_uncorrected);
  pion_rej_eff_graph->GetYaxis()->SetRangeUser(0,1.1);
  pion_rej_eff_graph->Draw("AP");
  elec_detec_eff_graph->Draw("P SAME");
  graphcanvas_uncorrected->Update();
  graphcanvas_uncorrected->SaveAs(Form("%s/pion_rej_eff.pdf",outputfilelocation.Data() ) );



  TCanvas *graphcanvas_corrected = new TCanvas("graphcanvas_corrected","graphcanvas_corrected",800,600);
  graphcanvas_corrected->SetGrid();
  utilityHandler.adjustCanvas(graphcanvas_corrected);
  pion_rej_eff_graph_corrected->GetYaxis()->SetRangeUser(0,1.1);
  elec_detec_eff_graph_corrected->Draw("AP"); 
  pion_rej_eff_graph_corrected->Draw("P SAME");
  graphcanvas_corrected->Update();
  //graphcanvas_corrected->SaveAs(Form("%s/pion_rej_eff.pdf",outputfilelocation.Data() ) );


  TCanvas *graphcanvas_both = new TCanvas("graphcanvas_both","graphcanvas_both",800,600);
  graphcanvas_corrected->SetGrid();
  utilityHandler.adjustCanvas(graphcanvas_both);
  pion_rej_eff_graph_corrected->GetYaxis()->SetRangeUser(0,1.1);
  elec_detec_eff_graph_corrected->GetYaxis()->SetRangeUser(0,1.1);
  elec_detec_eff_graph_corrected->Draw("AP"); 
  pion_rej_eff_graph_corrected->Draw("P SAME");
  pion_rej_eff_graph->Draw("P SAME");
  elec_detec_eff_graph->Draw("P SAME");
  graphcanvas_both->Update();
  //graphcanvas_corrected->SaveAs(Form("%s/pion_rej_eff.pdf",outputfilelocation.Data() ) );
  
 
  //// Make arrays that TGraphErrors can use 
  double x_PE[xMinSlices.size()];
  double y_PE[xMinSlices.size()];
  double x_PE_err[xMinSlices.size()];
  double y_PE_err[xMinSlices.size()];
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {
      double xCenter = xMinSlices[sliceid];
     
      x_PE[sliceid] = xCenter;
      y_PE[sliceid] = pions_to_electrons_ratio[sliceid];
      x_PE_err[sliceid] = 0;
      y_PE_err[sliceid] = 0;

    }    
  TGraphErrors *PE_graph = new TGraphErrors(xMinSlices.size(), x_PE, y_PE, x_PE_err, y_PE_err);
  utilityHandler.customizeGraphMore(PE_graph, 33, kBlue, 3,"",AxisTitle,"Pion to electron ratio per slice");

  //// canvas
  TCanvas *graphcanvas_PE = new TCanvas("graphcanvas_PE","graphcanvas_PE",800,600);
  graphcanvas_PE->SetGrid();
  utilityHandler.adjustCanvas(graphcanvas_PE);
  PE_graph->Draw("AP");
  graphcanvas_PE->Update();
  graphcanvas_PE->SaveAs(Form("%s/PE.pdf",outputfilelocation.Data() ) );
 
  //// extract the histogram title and print it
  std::string title = hist_data_orig->GetTitle();
  utilityHandler.printParsedTitle(title,outputfilelocation.Data(),"data");

  
  //  double initialParameters[7]={1,1,-0.02,-0.09,1,1,-1};
  //  // initialParameters = {0};
      
  //  TH1D *hist_data;
  //  /// Fit each slice 
  //  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
  //    {
  //      /// set up the global histograms that the fit function is going to use
  //      // this uses clones of the histograms
  //      hist_p  =(TH1D*)hist_vector_p[sliceid]->Clone("hist_p");
  //      hist_n  =(TH1D*)hist_vector_n[sliceid]->Clone("hist_n");
  //      hist_data  =(TH1D*)hist_vector_data[sliceid]->Clone("hist_data");

  //      //// this uses the histograms themselves: you can't uncouple the fit from them easily 
  //      //  hist_p  =hist_vector_p[sliceid];
  //      //  hist_n  =hist_vector_n[sliceid];
  //      // hist_data  =hist_vector_data[sliceid];

	 
  //      TF1 *Fit = new TF1(Form("overall_fit_%i_",sliceid),mc_p_n_poly2_slice_fit, xmin, xmax, 7);
  //      Fit->SetParameters(initialParameters);
  //       // set parameter limits. SetParLimits -> (par#, min, max)
  //      Fit->SetParLimits(0, 0, 2000); // scale_p greater than 0
  //      Fit->SetParLimits(1, 0,2000); // scale_n greater than 0
  //      Fit->SetParLimits(2, -0.10,0.10); // shift_p less than +- 10cm
  //      Fit->SetParLimits(3, -0.10,0.10); // shift_n less than +- 10cm
  //      Fit->SetParLimits(6,-10000000,-0.0000000001); // x^2 term negative to force downward concavity 
  //      Fit->SetNpx(500);
  //      hist_data->GetXaxis() ->SetRangeUser(xmin, xmax);
  //      hist_data->Fit(Fit,"R Q");

  //      // retrieve fit results 
  //      double scale_p  = Fit ->GetParameter(0);
  //      double scale_p_err = Fit ->GetParError(0);
  //      double scale_n  = Fit ->GetParameter(1);
  //      double scale_n_err = Fit ->GetParError(1);

  //      double shift_p= Fit ->GetParameter(2);
  //      double shift_p_err= Fit ->GetParError(2); 
  //      double shift_n = Fit ->GetParameter(3);
  //      double shift_n_err = Fit ->GetParError(3);
	  
  //      double ChiSq= Fit->GetChisquare();
  //      double ndf = Fit->GetNDF();

  //      //cout<<"shift_p in the loop: "<<shift_p<< " +/- "<<shift_p_err<<endl;
  //      //cout<<"shift_n in the loop: "<<shift_n<< " +/- "<<shift_n_err<<endl;

  //      std::vector<double> poly_result;
  //      std::vector<double> poly_result_err;
  //      for (int i =0 ; i < 3; i++)
  // 	{
  // 	  poly_result.push_back( Fit->GetParameter(4+i) );
  // 	  poly_result_err.push_back(Fit->GetParError(4+i) );
  // 	}

  //      //compute results
  //      double Rsf = scale_n/scale_p;
  //      double Rsf_err = Rsf * sqrt( pow( (scale_n_err / scale_n), 2) + pow( (scale_p_err / scale_p),2) ); //just adding the uncert from the fit parameters in quadrature for now. 

  //      cout<<"sliceid: "<<sliceid<< " ratio: scale_n / scale_p = "<< scale_n <<" / " <<scale_p <<" = "<< Rsf <<" +/- "<< Rsf_err<<endl;
  //      cout<<"sliceid: "<<sliceid<<" shift p: "<<shift_p<<" +/- "<<shift_p_err<<" shift_n: "<<shift_n<<"+/-"<<shift_n_err<<endl;
  //      cout<<endl;

  //      //// save results into the vectors 
  //      //fit_vector.push_back(Fit);
  //      scale_p_vector.push_back(scale_p);
  //      scale_n_vector.push_back(scale_n);
  //      scale_p_err_vector.push_back(scale_p_err);
  //      scale_n_err_vector.push_back(scale_n_err);
  //      shift_p_vector.push_back(shift_p);
  //      shift_n_vector.push_back(shift_n);
  //      shift_p_err_vector.push_back(shift_p_err);
  //      shift_n_err_vector.push_back(shift_n_err);
  //      ChiSq_vector.push_back(ChiSq);
  //      ndf_vector.push_back(ndf);
  //      Rsf_vector.push_back(Rsf);
  //      Rsf_err_vector.push_back(Rsf_err);
  //      // poly_result_vector_of_arrays.push_back(poly_result);
  //      poly_result_vector_of_vectors.push_back(poly_result);
  //      poly_result_err_vector_of_vectors.push_back(poly_result_err);
	  
  //      delete Fit; //
  //    }// end loop over slices


  //     // // Print the contents of the vector of vectors
  //     // for (const auto& vec : poly_result_vector_of_vectors) {
  //     //   for (double val : vec) {
  //     // 	 std::cout << val << " ";
  //     //   }
  //     //   std::cout << std::endl;
  //     // }



  //  double Rsf_mean = CalculateMean(Rsf_vector);
  //  double Rsf_stdev = CalculateStDev(Rsf_vector);
  //  double Rsf_mean_w = CalculateWeightedMean(Rsf_vector,Rsf_err_vector);
  //  double Rsf_stdev_w = CalculateWeightedStDev(Rsf_vector,Rsf_err_vector);
  //  cout<<"Rsf mean = " << Rsf_mean<<endl;
  //  cout<<"Rsf StDev = " << Rsf_stdev<<endl;
  //  cout<<"weighted Rsf mean = "<<Rsf_mean_w<<endl;
  //  cout<<"weighted Rsf StDev = "<<Rsf_stdev_w<<endl;
  //  cout<<endl;


 
  
  
  //  TH1D *hist_temp_p;
  //  TH1D *hist_temp_n;

  //  // loop over the slices to create shifted and scaled versions of the mc histograms to plot
  //  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
  //    {	
  //      /// shift the histograms. Function returns a clone. 
  //      hist_temp_p = shiftHistogramX(hist_vector_p[sliceid], shift_p_vector[sliceid] );
  //      hist_temp_n = shiftHistogramX(hist_vector_n[sliceid], shift_n_vector[sliceid] );

  //      /// scale the histograms 
  //      hist_temp_p -> Scale(scale_p_vector[sliceid]);
  //      hist_temp_n -> Scale(scale_n_vector[sliceid]);
	 
  //      // save histograms to the global vectors 
  //      hist_result_vector_p.push_back(hist_temp_p);
  //      hist_result_vector_n.push_back(hist_temp_n);	
  //    }// end loop over slices


  //  std::vector<TF1*> poly_fit_result;
  //  // loop over slices and make fits to plot the polynomial background result
  //  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
  //    {
  //      // get the fit results from the vector of vectors
  //      std::vector<double> poly_params_vector;
  //      poly_params_vector = poly_result_vector_of_vectors[sliceid];

  //      // convert the vector into an array that TF1 can use
  //      int size = poly_params_vector.size();
  //      double poly_params_array[size]; 
  //      for (int i  = 0 ; i < size; i++){
  // 	poly_params_array[i] = poly_params_vector[i];
  //      }

  //      TF1 *fit = new TF1(Form("poly2_%i_",sliceid),poly2, xmin, xmax, 3);
  //      fit->SetParameters(poly_params_array); 
  //      fit->SetNpx(500);
  //      poly_fit_result.push_back(fit);
  //    }// end loop over slices


  //     //// Plot Rsf
  //     //// Make arrays that TGraphErrors can use 
  //  double x[xMinSlices.size()];
  //  double y[xMinSlices.size()];
  //  double x_err[xMinSlices.size()];
  //  double y_err[xMinSlices.size()];
  //  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
  //    {
  //      double xCenter = (xMaxSlices[sliceid] + xMinSlices[sliceid] )/2;
  //      double xWidth = xMaxSlices[sliceid] - xMinSlices[sliceid];
     
  //      x[sliceid] = xCenter;
  //      y[sliceid] = Rsf_vector[sliceid];
  //      x_err[sliceid] = xWidth/2;
  //      y_err[sliceid] = Rsf_err_vector[sliceid];

  //      //cout<< "sliceid = "<<sliceid<<", xCenter = "<<x[sliceid]<<", x= "<<x[sliceid]<<", y = "<<y[sliceid]<<", x_err= "<<x_err[sliceid]<<", y_err = "<<y_err[sliceid]<<endl;
  //    }    
  //  TGraphErrors *Rsf_graph = new TGraphErrors(xMinSlices.size(), x, y, x_err, y_err);
  //  utilityHandler.customizeGraphMore(Rsf_graph, 33, kBlue, 3,"",AxisTitle,"Rsf");

  
  
  //  //// canvas
  //  TCanvas *graphcanvas = new TCanvas("graphcanvas","graphcanvas",800,600);
  //  graphcanvas->SetGrid();
  //  utilityHandler.adjustCanvas(graphcanvas);
  //  Rsf_graph->Draw("AP");
  //  graphcanvas->Update();
  //  graphcanvas->SaveAs(Form("%s/Rsf_xCenter.pdf",outputfilelocation.Data() ) );


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


  //   double x_max[xMinSlices.size()];
  //  double y_max[xMinSlices.size()];
  //  double x_max_err[xMinSlices.size()];
  //  double y_max_err[xMinSlices.size()];
  //  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
  //    {
  //      x_max[sliceid] = xMaxSlices[sliceid];
  //      y_max[sliceid] = Rsf_vector[sliceid];
  //      x_max_err[sliceid] = 0;
  //      y_max_err[sliceid] = Rsf_err_vector[sliceid];
  //    }    
  //  TGraphErrors *Rsf_xMax_graph = new TGraphErrors(xMinSlices.size(), x_max, y_max, x_max_err, y_max_err);
  //  utilityHandler.customizeGraphMore(Rsf_xMax_graph, 33, kBlue, 3,"",AxisTitle,"Rsf");
  
  // //// canvas
  //  TCanvas *Rsf_xMaxCanvas = new TCanvas("Rsf_xMaxCanvas","Rsf_xMaxCanvas",800,600);
  //  Rsf_xMaxCanvas->SetGrid();
  //  utilityHandler.adjustCanvas(Rsf_xMaxCanvas);
  //  Rsf_xMax_graph->Draw("AP");
  //  Rsf_xMaxCanvas->Update();
  //  Rsf_xMaxCanvas->SaveAs(Form("%s/Rsf_xMax.pdf",outputfilelocation.Data()));
  

  //   double x_m[xMinSlices.size()];
  //  double y_m[xMinSlices.size()];
  //  double x_m_err[xMinSlices.size()];
  //  double y_m_err[xMinSlices.size()];
  //  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
  //    {
  //      x_m[sliceid] = xMinSlices[sliceid];
  //      y_m[sliceid] = Rsf_vector[sliceid];
  //      x_m_err[sliceid] = 0;
  //      y_m_err[sliceid] = Rsf_err_vector[sliceid];
  //    }    
  //  TGraphErrors *Rsf_xMin_graph = new TGraphErrors(xMinSlices.size(), x_m, y_m, x_m_err, y_m_err);
  //  utilityHandler.customizeGraphMore(Rsf_xMin_graph, 33, kBlue, 3,"",AxisTitle,"Rsf");
  
  // //// canvas
  //  TCanvas *Rsf_xMinCanvas = new TCanvas("Rsf_xMinCanvas","Rsf_xMinCanvas",800,600);
  //  Rsf_xMinCanvas->SetGrid();
  //  utilityHandler.adjustCanvas(Rsf_xMinCanvas);
  //  Rsf_xMin_graph->Draw("AP");
  //  Rsf_xMinCanvas->Update();
  //  Rsf_xMinCanvas->SaveAs(Form("%s/Rsf_xMin.pdf",outputfilelocation.Data()));

  

  //  int nHist = hist_vector_data.size();
  //  int nCols = 4;
  //  int nRows = (nHist + nCols - 1) / nCols;


  //  std::vector<TH1D*> overall_fit_as_histogram;

  //  TCanvas* fits_canvas = new TCanvas("fits_canvas", "fits_canvas", 1000, 600);
  //  if (hist_vector_data.size()!=1)
  //    {
  //      fits_canvas->Divide(nCols, nRows);
  //    }
  //  for (int i = 0; i < nHist; ++i) {
  //    fits_canvas->cd(i + 1);
  //    hist_vector_data[i]->GetXaxis() ->SetRangeUser(xmin, xmax);
  //    hist_vector_data[i]->Draw();
      
  //    //// make a histogram that is the sum of the poly fit and the scaled mc histograms. 
  //    TH1D* sum_histo = sumHistogramsWithPolynomial(hist_result_vector_p[i],hist_result_vector_n[i] , poly_fit_result[i]);
  //    overall_fit_as_histogram.push_back(sum_histo);

  //    // // Create and customize the legend
  //    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  //    // legend->AddEntry(hist_vector_data[i], hist_vector_data[i]->GetName(), "l");
  //    legend->AddEntry("", Form("R= %.4f +/- %.4f ", Rsf_vector[i],Rsf_err_vector[i]), "");
  //    legend->AddEntry("", Form("#chi^{2}/ndf = %.2f / %.0f  ", ChiSq_vector[i] ,ndf_vector[i]), "");
  //    legend->Draw();

  //    poly_fit_result[i]->SetLineColor(kCyan);
  //    poly_fit_result[i]->Draw("same");
  //    hist_result_vector_p[i] ->SetLineColor(kGreen);
  //    hist_result_vector_p[i]->Draw("same");
  //    hist_result_vector_n[i] ->SetLineColor(kMagenta);
  //    hist_result_vector_n[i]->Draw("same");
  //    sum_histo ->SetLineColor(kRed);
  //    sum_histo->Draw("same");
     
  //  }// end loop over slices
  //  fits_canvas->Update();
  //   fits_canvas->SaveAs(Form("%s/fitted_slices.pdf",outputfilelocation.Data()));
  

  //  std::vector<TH1D*> hist_residual_vector;
 
  //  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
  //    {
  //      TH1D* hist_residual = (TH1D*)hist_vector_data[sliceid]->Clone("hist_residual");  
  //      hist_residual  ->Add(overall_fit_as_histogram[sliceid], -1);
  //      hist_residual->GetXaxis() ->SetRangeUser(xmin, xmax);
  //      hist_residual_vector.push_back(hist_residual);
  //    }// end loop over slices

   
  //  TCanvas* resid_canvas = new TCanvas("resid_canvas", "resid_canvas", 1000, 600);
  //  if (hist_vector_data.size()!=1)
  //    {
  //      resid_canvas->Divide(nCols, nRows);
  //    }
  //  for (int i = 0; i < nHist; ++i) {
  //    resid_canvas->cd(i + 1);
  //    hist_residual_vector[i]->Draw("E sames");
  //  }
  //  resid_canvas->Update();
  //  resid_canvas->SaveAs(Form("%s/residuals.pdf",outputfilelocation.Data() ));

     
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


  // Create a canvas to draw the histogram
  TCanvas* datacanvas1d = new TCanvas("datacanvas1d", "datacanvas1d", 800, 600);
  // Draw vertical lines at each x value in xSlices
  hist1D_data->SetLineWidth(2);
  hist1D_data->Draw();
  for (size_t i = 0; i < xMinSlices.size(); ++i) {
    int color = color_vector[i%(color_vector.size()-1)];
    double xLeft = xMinSlices[i];
    double xRight = xMaxSlices[i];
    TLine* lineLeft = new TLine(xLeft, hist1D_data->GetMinimum(), xLeft, hist1D_data->GetMaximum());
    TLine* lineRight = new TLine(xRight, hist1D_data->GetMinimum(), xRight, hist1D_data->GetMaximum());
    lineLeft->SetLineColor(color);  // Set line color (e.g., red)
    lineLeft->SetLineWidth(1);     // Set line width
    lineLeft->Draw("SAME");
    lineRight->SetLineColor(color);  // Set line color (e.g., red)
    lineRight->SetLineWidth(2);     // Set line width
    lineRight->Draw("SAME");
  }
  hist1D_data->Draw("SAME");
  
  datacanvas1d ->SaveAs(Form("%s/visulaized_cuts_1D_data.pdf",outputfilelocation.Data() ));

  
  //  // Create a projection of the TH2D histogram onto the x-axis
  //  TH1D* hist1D_data = hist_2D_data->ProjectionX("dataProjX");
  //  TH1D* hist1D_p = hist_2D_proton->ProjectionX("protonProjX");
  //  TH1D* hist1D_n = hist_2D_neutron->ProjectionX("neutronProjX");

  //  cout<<"data:"<<endl;
  //  cout<< "Mean = "<<hist1D_data->GetMean()<<", StdDev = "<<hist1D_data->GetStdDev()<<endl;
  //  cout<<"proton:"<<endl;
  //  cout<< "Mean = "<<hist1D_p->GetMean()<<", StdDev = "<<hist1D_p->GetStdDev()<<endl;
  //  cout<<"neutron:"<<endl;
  //  cout<< "Mean = "<<hist1D_n->GetMean()<<", StdDev = "<<hist1D_n->GetStdDev()<<endl;

  //  // Create a canvas to draw the histogram
  //  TCanvas* cutcanvas1D = new TCanvas("cutcanvas1D", "X Projection with Vertical Lines", 800, 600);
  //  cutcanvas1D->Divide(1,3);

  //  // Draw the TH1D histogram
  //  cutcanvas1D->cd(1);
  //  //hist1D_data->GetXaxis() ->SetRangeUser(0, 0.1);
  //  hist1D_data->Draw();
  //  // Draw vertical lines at each x value in xSlices
  //  for (size_t i = 0; i < xMinSlices.size(); ++i) {
  //    int color = color_vector[i%(color_vector.size()-1)];
  //    double xLeft = xMinSlices[i];
  //    double xRight = xMaxSlices[i];
  //    TLine* lineLeft = new TLine(xLeft, hist1D_data->GetMinimum(), xLeft, hist1D_data->GetMaximum());
  //    TLine* lineRight = new TLine(xRight, hist1D_data->GetMinimum(), xRight, hist1D_data->GetMaximum());
  //    lineLeft->SetLineColor(color);  // Set line color (e.g., red)
  //    lineLeft->SetLineWidth(2);     // Set line width
  //    lineLeft->Draw("SAME");
  //    lineRight->SetLineColor(color);  // Set line color (e.g., red)
  //    lineRight->SetLineWidth(2);     // Set line width
  //    lineRight->Draw("SAME");
  //  }

  //  // Draw the TH1D histogram
  //  cutcanvas1D->cd(2);
  //  //hist1D_p->GetXaxis() ->SetRangeUser(0, 0.1);
  //  hist1D_p->Draw("hist");
  //  // Draw vertical lines at each x value in xSlices
  //  for (size_t i = 0; i < xMinSlices.size(); ++i) {
  // int color = color_vector[i%(color_vector.size()-1)];
  //    double xLeft = xMinSlices[i];
  //    double xRight = xMaxSlices[i];
  //    TLine* lineLeft = new TLine(xLeft, hist1D_p->GetMinimum(), xLeft, hist1D_p->GetMaximum());
  //    TLine* lineRight = new TLine(xRight, hist1D_p->GetMinimum(), xRight, hist1D_p->GetMaximum());
  //    lineLeft->SetLineColor(color);  // Set line color (e.g., red)
  //    lineLeft->SetLineWidth(2);     // Set line width
  //    lineLeft->Draw("SAME");
  //    lineRight->SetLineColor(color);  // Set line color (e.g., red)
  //    lineRight->SetLineWidth(2);     // Set line width
  //    lineRight->Draw("SAME"); 
  //  }
  
  //  // Draw the TH1D histogram
  //  cutcanvas1D->cd(3);
  //  //hist1D_n->GetXaxis() ->SetRangeUser(0, 0.1);
  //  hist1D_n->Draw("hist");
  //  // Draw vertical lines at each x value in xSlices
  //  for (size_t i = 0; i < xMinSlices.size(); ++i) { 
  //     int color = color_vector[i%(color_vector.size()-1)];
  //    double xLeft = xMinSlices[i];
  //    double xRight = xMaxSlices[i];
  //    TLine* lineLeft = new TLine(xLeft, hist1D_n->GetMinimum(), xLeft, hist1D_n->GetMaximum());
  //    TLine* lineRight = new TLine(xRight, hist1D_n->GetMinimum(), xRight, hist1D_n->GetMaximum());
  //    lineLeft->SetLineColor(color);  // Set line color (e.g., red)
  //    lineLeft->SetLineWidth(2);     // Set line width
  //    lineLeft->Draw("SAME");
  //    lineRight->SetLineColor(color);  // Set line color (e.g., red)
  //    lineRight->SetLineWidth(2);     // Set line width
  //    lineRight->Draw("SAME"); 
  //  }

  //  cutcanvas1D->SaveAs(Form("%s/visulaized_cuts_1D.pdf",outputfilelocation.Data() ));
  

  //  // // TCanvas* testcanvas1  = new TCanvas("testcanvas1", "testcanvas1", 800, 600);
  //  // // TH1D* sumHistoTest = sumHistogramsWithPolynomial(hist_result_vector_p[0],hist_result_vector_n[0] , poly_fit_result[0]);
  //  // // hist_vector_data[0]->Draw();
  //  // // sumHistoTest ->Draw("same");

  
 
  //  // //// extract the histogram title and print it
  //  std::string title = hist_2D_data->GetTitle();
  //  cout<<"data title"<<endl;
  //  printParsedTitle(title,outputfilelocation);
  // std::string proton_title = hist_2D_proton->GetTitle();
  // cout<<"proton title: "<<endl;
  // cout<< proton_title<<endl;
  // std::string neutron_title =hist_2D_neutron->GetTitle();
  // cout<<"neutron title: "<<endl;
  // cout<< neutron_title<<endl;


  fout->Write(); 

  f1->Close();
  // f2->Close();
  //f3->Close();

}// End Main


