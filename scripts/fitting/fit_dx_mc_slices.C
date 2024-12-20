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
Double_t wide_x_range = 0.2;
Double_t nsigma = 0.6;
const int numbins = 200;

std::string AxisTitle ="hcal dx";

//// default values to be overwritten 
TString DataFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_tight_2Dhistos.root";
TString ProtonFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_tight_2Dhistos.root";
TString NeutronFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_tight_2Dhistos.root";
double xmin = -2.1; // -2.1 for SBS4 30p, maybe -2.5 for 50p
double xmax = 1.0;//1.4 for SBS4 40p
int PolyOrder =2;
TString config_title = "Testing out";



// SBS4 30p 
TString DataFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_2Dhistos_sept26.root";
//TString ProtonFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_2Dhistos_sept26.root";
TString ProtonFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_2Dhistos_100_199.root";
//TString NeutronFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_2Dhistos_sept26.root";
TString NeutronFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_2Dhistos_100_199.root";
int PolyOrder_sbs4_30p = 2;
double xmin_sbs4_30p = -2.15; 
double xmax_sbs4_30p = 1.4;
TString config_title_sbs4_30p = "SBS4 30%";

// SBS4 50p 
TString DataFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_2Dhistos_sept26.root";
TString ProtonFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deep_2Dhistos_sept26.root";
TString NeutronFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deen_2Dhistos_sept26.root";
int PolyOrder_sbs4_50p = 2;
double xmin_sbs4_50p = -2.63; 
double xmax_sbs4_50p = 1.4;
TString config_title_sbs4_50p = "SBS4 50%";

// SBS8 70p
TString DataFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_2Dhistos_sept26.root";
TString ProtonFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deep_2Dhistos_sept26.root";
TString NeutronFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deen_2Dhistos_sept26.root";
int PolyOrder_sbs8_70p = 4;
double xmin_sbs8_70p = -1.9; 
double xmax_sbs8_70p = 1.0;
TString config_title_sbs8_70p = "SBS8 70%";

// SBS8 50p
TString DataFileString_sbs8_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_50p_cuts_2Dhistos_sept26.root";
TString ProtonFileString_sbs8_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_50p_cuts_deep_2Dhistos_sept26.root";
TString NeutronFileString_sbs8_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_50p_cuts_deen_2Dhistos_sept26.root";
int PolyOrder_sbs8_50p = 4;
double xmin_sbs8_50p = -1.6; 
double xmax_sbs8_50p = 1.0;
TString config_title_sbs8_50p = "SBS8 50%";


// SBS8 100p
TString DataFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_2Dhistos_sept24.root";
TString ProtonFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_deep_2Dhistos_sept24.root";
TString NeutronFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_deen_2Dhistos_sept24.root";
int PolyOrder_sbs8_100p = 4;
double xmin_sbs8_100p = -2.2; 
double xmax_sbs8_100p = 1.0;
TString config_title_sbs8_100p = "SBS8 100%";


// SBS9 70p
TString DataFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_2Dhistos_sept26_1.root";
TString ProtonFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_deep_2Dhistos_sept26_z_1.root";
TString NeutronFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_deen_2Dhistos_sept26_z_1.root";
int PolyOrder_sbs9_70p = 4;
double xmin_sbs9_70p = -1.9; 
double xmax_sbs9_70p = 1.0;
TString config_title_sbs9_70p = "SBS9 70%";


// Function Declarations



// MAIN
void fit_dx_mc_slices(TString KineString="sbs4_30p",int PolyOrder_input = 2){//main
  //gStyle->SetOptFit(0);
  // gStyle->SetCanvasPreferGL(1);
  // gStyle -> SetOptStat(0);
  //gStyle ->SetEndErrorSize(0);

 

  //// Histogram to store the Rsf distriubtion
  TH1D *h_Rsf = new TH1D("h_Rsf","h_Rsf",300, 0.8,1.2);

  /// Histogram for the pull distribution
  TH1D *h_pull = new TH1D("h_pull","h_pull",10,-1,1);

  Utility utilityHandler; // class that gives us access to various functions to use between programs.

  FileNames fileNamesHandler;

  if (KineString == "sbs4_30p"){
    DataFileString = DataFileString_sbs4_30p;
    //ProtonFileString =ProtonFileString_sbs4_30p;
    //NeutronFileString = NeutronFileString_sbs4_30p;
    xmin = fileNamesHandler.xmin_sbs4_30p;
    xmax =fileNamesHandler.xmax_sbs4_30p;
    config_title = config_title_sbs4_30p;
  }
  // else if (KineString == "sbs4_50p"){
  //   DataFileString = DataFileString_sbs4_50p;
  //   ProtonFileString =ProtonFileString_sbs4_50p;
  //   NeutronFileString = NeutronFileString_sbs4_50p;
  //   xmin = xmin_sbs4_50p;
  //   xmax =xmax_sbs4_50p;
  //   config_title = config_title_sbs4_50p;
  //}
  else if (KineString == "sbs8_70p"){
    DataFileString = DataFileString_sbs8_70p;
    // ProtonFileString =ProtonFileString_sbs8_70p;
    //NeutronFileString = NeutronFileString_sbs8_70p;
    xmin = fileNamesHandler.xmin_sbs8_70p;
    xmax =fileNamesHandler.xmax_sbs8_70p;
    config_title = config_title_sbs8_70p;
  }
//else if (KineString == "sbs8_50p"){
  //   DataFileString = DataFileString_sbs8_50p;
  //   ProtonFileString =ProtonFileString_sbs8_50p;
  //   NeutronFileString = NeutronFileString_sbs8_50p;
  //   ProtonFileString =ProtonFileString_sbs8_50p;
  //   NeutronFileString = NeutronFileString_sbs8_50p;
  //   xmin = xmin_sbs8_50p;
  //   xmax =xmax_sbs8_50p;
  //   config_title = config_title_sbs8_50p;
  // }else if (KineString == "sbs8_100p"){
  //   DataFileString = DataFileString_sbs8_100p;
  //   ProtonFileString =ProtonFileString_sbs8_100p;
  //   NeutronFileString = NeutronFileString_sbs8_100p;
  //   ProtonFileString =ProtonFileString_sbs8_100p;
  //   NeutronFileString = NeutronFileString_sbs8_100p;
  //   xmin = xmin_sbs8_100p;
  //   xmax =xmax_sbs8_100p;
  //   config_title = config_title_sbs8_100p;
  // }
  else if (KineString == "sbs9_70p"){
    DataFileString = DataFileString_sbs9_70p;
    //ProtonFileString =ProtonFileString_sbs9_70p;
    //NeutronFileString = NeutronFileString_sbs9_70p;
    xmin = fileNamesHandler.xmin_sbs9_70p;
    xmax =fileNamesHandler.xmax_sbs9_70p;
    config_title = config_title_sbs9_70p;
  }
  else {
    std::cout<<"Error with kinematic setting"<<std::endl;
    return;
  }

  if((PolyOrder_input < 0) && PolyOrder_input > 6){
    cout<<"error with poly order input. Setting to 2."<<endl;
    PolyOrder = 2;
  }else PolyOrder=PolyOrder_input; 

  cout<<endl;
  cout<<"What we're running: "<<KineString<<endl;
  cout<<DataFileString<<endl;
  //cout<<ProtonFileString<<endl;
  //cout<<NeutronFileString<<endl;
  cout<<"poly order: "<<PolyOrder <<endl;
  cout<<"dx range: ("<<xmin<<", "<<xmax<<")"<<endl;
  cout<<endl;

  int nGroups;
  TString base_file_string_p;
  TString base_file_string_n;
  std::vector<TString> n_filename;
  std::vector<TString> p_filename;
  
    if (KineString == "sbs4_30p"){
      // Set up file paths
      nGroups = 5;
      base_file_string_p = "/lustre24/expphy/volatile/halla/sbs/msatnik/output/sbs4_30p_cuts_deep_2Dhistos_";
      base_file_string_n = "/lustre24/expphy/volatile/halla/sbs/msatnik/output/sbs4_30p_cuts_deen_2Dhistos_";
      n_filename = {
	base_file_string_n + "0_99.root",
	base_file_string_n + "100_199.root",
	base_file_string_n + "200_299.root",
	base_file_string_n + "300_387.root",
	base_file_string_n + "389_499.root"
      };
      p_filename = {
	base_file_string_p + "0_99.root",
	base_file_string_p + "100_199.root",
	base_file_string_p + "200_299.root",
	base_file_string_p + "300_387.root",
	base_file_string_p + "389_499.root"
      };
    }else if(KineString == "sbs8_70p"){
      // Set up file paths
      nGroups = 3;
      base_file_string_p = "/lustre24/expphy/volatile/halla/sbs/msatnik/output/sbs8_70p_cuts_deep_subset_2Dhistos_";
      base_file_string_n = "/lustre24/expphy/volatile/halla/sbs/msatnik/output/sbs8_70p_cuts_deen_subset_2Dhistos_";
      n_filename = {
	base_file_string_n + "66_132.root",
	base_file_string_n + "133_199.root",
	base_file_string_n + "0_65.root"
      };
      p_filename = {
	base_file_string_p + "66_132.root",
	base_file_string_p + "133_199.root",
	base_file_string_p + "0_65.root"
      };
    }else if(KineString == "sbs9_70p"){
      // Set up file paths
      nGroups = 3;
      base_file_string_p = "/lustre24/expphy/volatile/halla/sbs/msatnik/output/sbs9_70p_cuts_deep_2Dhistos_";
      base_file_string_n = "/lustre24/expphy/volatile/halla/sbs/msatnik/output/sbs9_70p_cuts_deen_2Dhistos_";
      n_filename = {
	base_file_string_n + "0_76.root",
	base_file_string_n + "77_153.root",
	base_file_string_n + "154_229.root"
      };
      p_filename = {
	base_file_string_p + "0_76.root",
	base_file_string_p + "77_153.root",
	base_file_string_p + "154_229.root"
      };
    }

    
    
  for (int i = 0; i<p_filename.size();i++)
    {
      cout<<"Proton files: "<<endl;
      cout<<p_filename[i]<<endl;
    }
  cout<<endl;
  for (int i = 0; i<p_filename.size();i++)
    {
      cout<<"Neutron files: "<<endl;
      cout<<n_filename[i]<<endl;
    }


  // testing stdev
  std::vector<double> test1;
  std::vector<double> test2;

  test1 = { 5, 20 , 40};
  test2 = {40, 20, 5};

  double test1_stdev = utilityHandler.CalculateStDev(test1);
  double test2_stdev = utilityHandler.CalculateStDev(test2);
  cout<<endl;
  cout<<"testing stdev because I feel like something is wrong"<<endl;
  for (int i  = 0; i< test1.size();i++){
    cout<<test1[i]<<", ";
  }
  cout<<endl;
  cout << "StDev = "<<test1_stdev<<endl;
  
    for (int i  = 0; i< test2.size();i++){
    cout<<test2[i]<<", ";
  }
  cout<<endl;
  cout << "StDev = "<<test2_stdev<<endl;
  
  
  
  
  // Vectors to store the histograms
  std::vector<TH1D*> hist_vector_p;
  std::vector<TH1D*> hist_vector_n;

  TString hist_name = "hcal_dx_1d_allcuts";

   
  // TString DataFileString  = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_2Dhistos_sept26.root";
  

  // // Load data file and histograms 
  // TFile *data_file  = TFile::Open(DataFileString,"READ"); // data
   
  // TH1D *hist_data_orig = (TH1D*)data_file->Get(hist_name);
  // hist_data_orig->SetDirectory(0);  // Detach histogram from file so it persists after file is closed
  // TH1D *hist_data = (TH1D*)hist_data_orig->Clone("hist_data");
  // TH1D *hist_data_dont_draw = (TH1D*)hist_data_orig->Clone("hist_data_dont_draw");

  TH1D *hist_data_orig, *hist_data, *hist_data_dont_draw;
  
  // Load data file and histograms 
  TFile *data_file  = TFile::Open(DataFileString,"READ"); // data
  if (data_file && !data_file->IsZombie()) {
    hist_data_orig = (TH1D*)data_file->Get(hist_name);
    if (hist_data_orig) {
      hist_data = (TH1D*)hist_data_orig->Clone("hist_data");
      hist_data_dont_draw = (TH1D*)hist_data_orig->Clone("hist_data_dont_draw");
      hist_data_orig->SetDirectory(0);  // Detach histogram from file so it persists after file is closed
    } else {
      std::cerr << "Warning: Histogram " << hist_name << " not found in file " << DataFileString << std::endl;
    }
    // data_file->Close();
  } else {
    std::cerr << "Error: Unable to open file " << DataFileString << std::endl;
  }
     
    
  // Loop over p files
  for (int i = 0; i < p_filename.size(); i++) {
    TFile *p_file = TFile::Open(p_filename[i], "READ");
    if (p_file && !p_file->IsZombie()) {
      TH1D *hist_p_orig = (TH1D*)p_file->Get(hist_name);
      if (hist_p_orig) {
	TH1D *p_hist = (TH1D*)hist_p_orig->Clone(Form("p_hist_%d", i));
	hist_vector_p.push_back(p_hist);
	hist_p_orig->SetDirectory(0);  // Detach histogram from file so it persists after file is closed
      } else {
	std::cerr << "Warning: Histogram " << hist_name << " not found in file " << p_filename[i] << std::endl;
      }
      //p_file->Close();
    } else {
      std::cerr << "Error: Unable to open file " << p_filename[i] << std::endl;
    }
  }

  // Loop over n files
  for (int i = 0; i < n_filename.size(); i++) {
    TFile *n_file = TFile::Open(n_filename[i], "READ");
    if (n_file && !n_file->IsZombie()) {
      TH1D *hist_n_orig = (TH1D*)n_file->Get(hist_name);
      if (hist_n_orig) {
	TH1D *n_hist = (TH1D*)hist_n_orig->Clone(Form("n_hist_%d", i));
	hist_vector_n.push_back(n_hist);
	hist_n_orig->SetDirectory(0);  // Detach histogram from file so it persists after file is closed
      } else {
	std::cerr << "Warning: Histogram " << hist_name << " not found in file " << n_filename[i] << std::endl;
      }
      // n_file->Close();
    } else {
      std::cerr << "Error: Unable to open file " << n_filename[i] << std::endl;
    }
  }

  // Do something with histograms if needed
  std::cout << "Loaded " << hist_vector_p.size() << " p histograms and " << hist_vector_n.size() << " n histograms." << std::endl;

  if (hist_vector_p.size()!=hist_vector_n.size())
    {
      cout<<"Error. Different number of proton and neutron histograms loaded"<<endl;
      return; 
    }
  
 
 
  // set location and name for output file 
  TString outputfilelocation="../../output/mc_slice_study/"+KineString+"/poly"+PolyOrder;
  TString outputfilename = outputfilelocation +"/mc_slice_study_"+ KineString+"_poly"+PolyOrder+".root";
  
  // Declare outfile
  TFile *fout = new TFile(outputfilename,"RECREATE");
  cout<<"writing to file: "<< outputfilename <<endl;



  std::vector<TF1*> fit_vector;

  // std::vector<TH1D*> hist_result_p_vector;
  // std::vector<TH1D*> hist_result_n_vector;

  std::vector<TH1D*> hist_result_vector_p; 
  std::vector<TH1D*> hist_result_vector_n; 

  std::vector<TH1D*> hist_residual_vector; 

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
  hist_result_vector_p.clear();
  hist_result_vector_n.clear();
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

  int nEntries_proton_sum =0;
  int nEntries_neutron_sum = 0;
  

  double initialParameters[7]={1,1,0,0,1,1,-1};
  for (int sliceid = 0 ; sliceid <  hist_vector_n.size() ; sliceid++)
    {
      TH1D *hist_p  =(TH1D*)hist_vector_p[sliceid]->Clone("hist_p");
      TH1D *hist_n  =(TH1D*)hist_vector_n[sliceid]->Clone("hist_n");
      if (!hist_data_dont_draw || !hist_p || !hist_n) {
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
      
      int nEntries_data = hist_data->GetEntries();
      int nEntries_neutron = hist_n->GetEntries();
      int nEntries_proton = hist_p ->GetEntries();

      nEntries_proton_sum += nEntries_proton;
      nEntries_neutron_sum +=nEntries_neutron;
      
      cout<<endl;
      cout<<"sliceid: "<<sliceid<<endl;
      cout<<"nEntries in data histogram: "<<nEntries_data<<endl;
      cout<<"nEntries in proton histogram (before scaling): "<<nEntries_proton<<endl;
      cout<<"nEntries in neutron histogram (before scaling): "<<nEntries_neutron<<endl;;
      cout<<"shift_p in the loop: "<<shift_p<< " +/- "<<shift_p_err<<endl;
      cout<<"shift_n in the loop: "<<shift_n<< " +/- "<<shift_n_err<<endl;
    
       
      std::vector<double> poly_result;
      std::vector<double> poly_result_err;
      for (int i =0 ; i <=PolyOrder; i++)
	{
	  poly_result.push_back(fitHandler.poly_result[i]);
	  poly_result_err.push_back(fitHandler.poly_result_err[i]);
	}

      h_Rsf->Fill(Rsf);

      cout<<"scale_p = "<< scale_p <<endl;
      cout<<"Rsf = "<< Rsf <<" +/- "<< Rsf_err<<endl;
      cout<<"shift p: "<<shift_p<<" +/- "<<shift_p_err<<", shift_n: "<<shift_n<<"+/-"<<shift_n_err<<endl;
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
       
    }// end for loop over slices

  cout<<"sum of nEntries for proton before scaling: "<<nEntries_proton_sum<<endl;
  cout<<"sum of nEntreis for neutron before scaling: "<<nEntries_neutron_sum<<endl;  
  
  cout<<"---------------------------------------------------------------------------------------------"<<endl;
  double Rsf_mean = utilityHandler.CalculateMean(Rsf_vector);
  double Rsf_stdev = utilityHandler.CalculateStDev(Rsf_vector);
  double Rsf_mean_w = utilityHandler.CalculateWeightedMean(Rsf_vector,Rsf_err_vector);
  double Rsf_stdev_w = utilityHandler.CalculateWeightedStDev(Rsf_vector,Rsf_err_vector);
  cout<<"Rsf mean = " << Rsf_mean<<endl;
  cout<<"Rsf StDev = " << Rsf_stdev<<endl;
  //cout<<"weighted Rsf mean = "<<Rsf_mean_w<<endl;
  //cout<<"weighted Rsf StDev = "<<Rsf_stdev_w<<endl;
  cout<<"StDev/sqrt("<<Rsf_vector.size()<<") = "<<Rsf_stdev / sqrt(Rsf_vector.size())<<endl;
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
  for (int sliceid = 0 ; sliceid < hist_vector_p.size() ; sliceid++)
    {	
    
      // utility function that scales and shifts the provided histogram and returns a new histogram 
      hist_temp_p = utilityHandler.ScaleAndShiftHistogram(hist_vector_p[sliceid], scale_p_vector[sliceid], shift_p_vector[sliceid]);
      hist_temp_n = utilityHandler.ScaleAndShiftHistogram(hist_vector_n[sliceid], scale_n_vector[sliceid], shift_n_vector[sliceid]);
	 
      // save histograms to the global vectors 
      hist_result_vector_p.push_back(hist_temp_p);
      hist_result_vector_n.push_back(hist_temp_n);	
    }// end loop over slices

  std::vector<TF1*> poly_fit_result;
  // loop over slices and make fits to plot the polynomial background result
  for (int sliceid = 0 ; sliceid < hist_vector_p.size() ; sliceid++)
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
  double x[hist_vector_p.size()];
  double y[hist_vector_p.size()];
  double x_err[hist_vector_p.size()];
  double y_err[hist_vector_p.size()];
  for (int sliceid = 0 ; sliceid < hist_vector_p.size() ; sliceid++)
    {
      double xCenter = sliceid+1;
      double xWidth = 0;
     
      x[sliceid] = xCenter;
      y[sliceid] = Rsf_vector[sliceid];
      x_err[sliceid] =0;
      y_err[sliceid] = Rsf_err_vector[sliceid];

      //cout<< "sliceid = "<<sliceid<<", xCenter = "<<x[sliceid]<<", x= "<<x[sliceid]<<", y = "<<y[sliceid]<<", x_err= "<<x_err[sliceid]<<", y_err = "<<y_err[sliceid]<<endl;
    }    
  TGraphErrors *Rsf_graph = new TGraphErrors(hist_vector_p.size(), x, y, x_err, y_err);
  utilityHandler.customizeGraphMore(Rsf_graph, 33, kBlue, 3,"","Group #","Rsf");

  // Fit a straight line to the graph
  TF1* fit_Rsf_graph = new TF1("fit_Rsf_graph", "[0]",Rsf_graph ->GetX()[0], Rsf_graph->GetX()[Rsf_graph->GetN()-1]);
  Rsf_graph->Fit(fit_Rsf_graph, "Q R");
  
  //// canvas
  TCanvas *graphcanvas = new TCanvas("graphcanvas","graphcanvas",800,600);
  graphcanvas->SetGrid();
  utilityHandler.adjustCanvas(graphcanvas);
  //Rsf_graph->GetXaxis()->SetRangeUser(-1,5);
  Rsf_graph->Draw("AP");

  // Get fit parameters
  double constant = fit_Rsf_graph->GetParameter(0);       // The constant value
  double constantError = fit_Rsf_graph->GetParError(0);    // The error on the constant
  double chi2 = fit_Rsf_graph->GetChisquare();             // The chi-squared value
  int ndf = fit_Rsf_graph->GetNDF();                       // The number of degrees of freedom
  double chi2_ndf = chi2 / ndf;                          // chi2/ndf

  // Draw the graph with the fit result again
  Rsf_graph->Draw("AP");
  graphcanvas->Update();

  // Use TLatex to add the fit result and chi²/ndf to the canvas
  TLatex latex;
  latex.SetNDC();  // Use normalized coordinates (0 to 1)
  latex.SetTextSize(0.04);
  //latex.DrawLatex(0.2, 0.85, Form("FIT: R_{sf} = %.5f #pm %.5f, #chi^{2}/ndf = %.2f/%i = %.2f", constant, constantError, chi2,ndf, chi2_ndf));
  //latex.DrawLatex(0.2,0.80, Form("(Fit Uncert)/#sqrt{%zu} = %.5f", hist_vector_p.size(),constantError / sqrt(5) ) );

   
  latex.DrawLatex(0.2, 0.25, Form("R_{sf} Mean = %.6f #pm StDev = %.6f", Rsf_mean, Rsf_stdev));
  latex.DrawLatex(0.2,0.20, Form("StDev/#sqrt{%zu} = %.6f",hist_vector_p.size(), Rsf_stdev / sqrt(hist_vector_p.size() ) ));

  
   
  graphcanvas->Update();
  graphcanvas->SaveAs(Form("%s/Rsf_xCenter.pdf",outputfilelocation.Data() ) );


//// Plot Chi2/ndf
  //// Make arrays that TGraphErrors can use
  
  double x_ch[hist_vector_p.size()];
  double y_ch[hist_vector_p.size()];
  double x_ch_err[hist_vector_p.size()];
  double y_ch_err[hist_vector_p.size()];
  for (int sliceid = 0 ; sliceid < hist_vector_p.size() ; sliceid++)
    {
      double xCenter = sliceid +1;
      double xWidth = 0;
      
      x_ch[sliceid] = xCenter;
      y_ch[sliceid] = ChiSq_vector[sliceid] /ndf_vector[sliceid];
      x_ch_err[sliceid] = xWidth;
      y_ch_err[sliceid] = 0;
    }    
  TGraphErrors *Chi2_ndf_graph = new TGraphErrors(hist_vector_p.size(), x_ch, y_ch, x_ch_err, y_ch_err);
  utilityHandler.customizeGraphMore(Chi2_ndf_graph, 33, kBlue, 3,"","Group #","Chi^{2}/ndf");


  //// canvas
  TCanvas *chi2_ndf_canvas = new TCanvas("chi2_ndf_canvas","chi2_ndf_canvas",800,600);  
  utilityHandler.adjustCanvas(chi2_ndf_canvas);
  chi2_ndf_canvas ->SetGrid();
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
  double x_nEntries_p[hist_vector_p.size()];
  double y_nEntries_p[hist_vector_p.size()];
  double x_nEntries_p_err[hist_vector_p.size()];
  double y_nEntries_p_err[hist_vector_p.size()];
  for (int sliceid = 0 ; sliceid < hist_vector_p.size() ; sliceid++)
    {
      double xCenter_nEntries_p = sliceid +1;
      double xWidth_nEntries_p =0;
     
      x_nEntries_p[sliceid] = xCenter_nEntries_p;
      y_nEntries_p[sliceid] = hist_vector_p[sliceid]->GetEntries();
      x_nEntries_p_err[sliceid] = xWidth_nEntries_p/2;
      y_nEntries_p_err[sliceid] = 0;
    }    
  TGraphErrors *nEntries_p_graph = new TGraphErrors(hist_vector_p.size(), x_nEntries_p, y_nEntries_p, x_nEntries_p_err, y_nEntries_p_err);
  utilityHandler.customizeGraphMore(nEntries_p_graph, 33, kGreen, 3,"",AxisTitle,"nEntries_p");

  //// Plot nEntries
  //// Make arrays that TGraphErrors can use 
  double x_nEntries_n[hist_vector_p.size()];
  double y_nEntries_n[hist_vector_p.size()];
  double x_nEntries_n_err[hist_vector_p.size()];
  double y_nEntries_n_err[hist_vector_p.size()];
  for (int sliceid = 0 ; sliceid < hist_vector_p.size() ; sliceid++)
    {
      double xCenter_nEntries_n = sliceid+1;
      double xWidth_nEntries_n =0;
     
      x_nEntries_n[sliceid] = xCenter_nEntries_n;
      y_nEntries_n[sliceid] = hist_vector_n[sliceid]->GetEntries();
      x_nEntries_n_err[sliceid] = xWidth_nEntries_n/2;
      y_nEntries_n_err[sliceid] = 0;
    }    
  TGraphErrors *nEntries_n_graph = new TGraphErrors(hist_vector_p.size(), x_nEntries_n, y_nEntries_n, x_nEntries_n_err, y_nEntries_n_err);
  utilityHandler.customizeGraphMore(nEntries_n_graph, 33, kMagenta, 3,"",AxisTitle,"nEntries_n");
  
  //// canvas
  TCanvas *nEntriesCanvas = new TCanvas("nEntriesCanvas","nEntriesCanvas",800,600);  nEntriesCanvas->SetGrid();
  utilityHandler.adjustCanvas(nEntriesCanvas);
  nEntries_p_graph ->GetYaxis()->SetRangeUser(35000,46000);
  nEntries_p_graph->Draw("AP");
  nEntries_n_graph ->GetYaxis()->SetRangeUser(35000,46000);
  nEntries_n_graph->Draw("P SAME");
  
  nEntriesCanvas->Update();
  nEntriesCanvas->SaveAs(Form("%s/nEntries.pdf",outputfilelocation.Data()));

 int nHist = hist_vector_p.size();
  int nCols = 4;
  int nRows = (nHist + nCols - 1) / nCols;


  std::vector<TH1D*> overall_fit_as_histogram;

  TCanvas* fits_canvas = new TCanvas("fits_canvas", "fits_canvas", 1000, 600);
  if (hist_vector_p.size()!=1)
    {
      fits_canvas->Divide(nCols, nRows);
    }
   hist_data->GetXaxis() ->SetRangeUser(xmin, xmax);
  for (int i = 0; i < nHist; ++i) {
    fits_canvas->cd(i + 1);
    fits_canvas->SetGrid();
    hist_data->Draw();
    fits_canvas->Update();

    int nEntries = hist_vector_p[i]->GetEntries() + hist_vector_n[i]->GetEntries();
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
 
 
  for (int sliceid = 0 ; sliceid < hist_vector_p.size() ; sliceid++)
    {
      TH1D* hist_residual = (TH1D*)hist_data->Clone(Form("hist_residual_%d",sliceid));  
      hist_residual  ->Add(overall_fit_as_histogram[sliceid], -1);
      hist_residual->GetXaxis() ->SetRangeUser(xmin, xmax);
      hist_residual_vector.push_back(hist_residual);
    }// end loop over slices

   
  TCanvas* resid_canvas = new TCanvas("resid_canvas", "resid_canvas", 1000, 600);
  if (hist_vector_p.size()!=1)
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

     
// //// extract the histogram title and print it
  std::string title = hist_data->GetTitle();
  cout<<"data title"<<endl;
  utilityHandler.printParsedTitle(title,outputfilelocation,"data");
  std::string proton_title = hist_vector_p[0]->GetTitle();
  cout<<"proton title: "<<endl;
  cout<< proton_title<<endl;
  std::string neutron_title =hist_vector_n[0]->GetTitle();
  cout<<"neutron title: "<<endl;
  cout<< neutron_title<<endl;

  
  /// Make a canvas called c1 to overwrite the default canvas and avoid crashes
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  hist_data ->Draw();
  c1->Update();
     

 
  fout->Write(); 
  
}// end main



