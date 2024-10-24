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
TString HistogramName = "hcal_dx__hcal_nsigdy";
std::string AxisTitle ="hcal_nsigdy";
TString note = "";
std::string xMin_xMax_string = "(-4.1,-2.1),(-3.1,-1.1),(-2,0.1),(-1,1),(0,2),(1.1,3),(2.1,4)";
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


//// ranges for slices. These are what we want to change.
////****** dy *************************************************************** 
//std::string xMin_xMax_string_dy = "(-1,-0.5),(-0.5,0),(0,0.5),(0.5,1)";
//std::string xMin_xMax_string_dy = "(-0.9,-0.7),(-0.7,-0.5)(-0.5,-0.3),(-0.3,-0.1)(-0.1,0.1),(0.1,0.3),(0.3,0.5),(0.5,0.7),(0.7,0.9)";
std::string xMin_xMax_string_dy = "(-0.5,-0.3),(-0.3,-0.1)(-0.1,0.1),(0.1,0.3),(0.3,0.5)";
//std::string xMin_xMax_string_dy = "(-0.5,-0.4),(-0.4,-0.3)(-0.3,-0.2),(-0.2,-0.1),(-0.1,0),(0,0.1),(0.1,0.2),(0.2,0.3),(0.3,0.4)(0.4,0.5)";
//std::string xMin_xMax_string_dy = "(-1,-0.5),(-0.75,-0.25),(-0.5,0),(-0.25,0.25),(0,0.5),(0.25,0.75),(0.5,1)";
//std::string xMin_xMax_string_dy = "(-1,1),(-0.9,0.9),(-0.8,0.8),(-0.7,0.7),(-0.6,0.6),(-0.5,0.5),(-0.4,0.4),(-0.3,0.3),(-0.2,0.2),(-0.1,0.1),(-0.05,0.05)";
////*****nsigdy***************************************************************
//std::string xMin_xMax_string_nsigdy = "(-4.1,-2.1),(-3.1,-1.1),(-2,0.1),(-1,1),(0,2),(1.1,3),(2.1,4)";
std::string xMin_xMax_string_nsigdy = "(-4,-3),(-3,-2),(-1,0),(0,1),(1,2),(2,3),(3,4)";
///********vz**********************************************************************
//std::string xMin_xMax_string_vz = "(-0.08,0.08),(-0.07,0.07),(-0.065,0.065),(-0.06,0.06),(-0.05,0.05),(-0.04,0.04),(-0.03,0.03)";
//std::string xMin_xMax_string_vz = "(-0.08,-0.07),(-0.07,-0.06),(-0.06,-0.05),(-0.05,-0.04),(-0.04,-0.03),(-0.03,-0.02),(-0.02,-0.01),(-0.01,0),(0,0.01),(0.01,0.02),(0.02,0.03),(0.03,0.04),(0.04,0.05),(0.05,0.06),(0.06,0.07),(0.07,0.08)";
//std::string xMin_xMax_string_vz = "(-0.08,-0.06),(-0.06,-0.04),(-0.04,-0.02),(-0.02,0.00),(0.00,0.02),(0.02,0.04),(0.04,0.06),(0.06,0.08)";
std::string xMin_xMax_string_vz = "(-0.07,-0.05),(-0.05,-0.03),(-0.03,-0.01),(-0.01,0.01),(0.01,0.03),(0.03,0.05),(0.05,0.07)";
//std::string xMin_xMax_string_vz = "(-0.08,-0.06),(-0.07,-0.05),(-0.06,-0.04),(-0.05,-0.03),(-0.04,-0.02),(-0.03,-0.01),(-0.02,0),(-0.01,0.01),(0,0.02),(0.01,0.03),(0.02,0.04),(0.03,0.05),(0.04,0.06),(0.05,0.07),(0.06,0.08)";
//std::string xMin_xMax_string_vz = "(-0.065,0.065)";
//// ******ps_e**********************************************************************
//std::string xMin_xMax_string_ps =  "(0.1,0.2),(0.12,0.22),(0.14,0.24),(0.16,0.26),(0.18,0.28),(0.2,0.3),(0.22,0.32),(0.24,0.34),(0.26,0.36),(0.28,0.38)";
//std::string xMin_xMax_string_ps =  "(0.1,2),(0.12,2),(0.14,2),(0.16,2),(0.18,2),(0.2,2),(0.22,2),(0.24,2),(0.26,2),(0.28,2)";
//std::string xMin_xMax_string_ps ="(0.2,2)";
std::string xMin_xMax_string_ps =  "(0.10,0.15),(0.15,0.20),(0.20,0.25),(0.25,0.3),(0.30,0.35),(0.40,0.45)";
////******* hcal energy **************************************************************
//std::string xMin_xMax_string_hcal_e = "(0.015,0.025), (0.020,0.030), (0.025,0.035),(0.030,0.040)";
//std::string xMin_xMax_string_hcal_e = "(0.015,1), (0.020,1), (0.025,1),(0.030,1)";
std::string xMin_xMax_string_hcal_e = "(0.010,0.015), (0.015,0.020), (0.020,0.025), (0.025,0.030), (0.030,0.035), (0.035,0.040), (0.040,0.045), (0.045,0.050),(0.050,0.055),(0.055,0.060)";

//// ******w2****************************************************************
//std::string xMin_xMax_string_w2 = "(0.4,0.6),(0.5,0.7),(0.6,0.8),(0.7,0.9),(0.8,1.0),(0.9,1.1),(1.0,1.2),(1.1,1.3) ";
//std::string xMin_xMax_string_w2 = "(0.4,1.2),(0.5,1.1),(0.6,1.0),(0.7,0.9)";
//std::string xMin_xMax_string_w2 = "(0.78,0.98),(0.68,1.08),(0.58,1.18)";
std::string xMin_xMax_string_w2 = "(0.4,0.5),(0.5,0.6),(0.6,0.7),(0.7,0.8),(0.8,0.9),(0.9,1.0),(1.0,1.1),(1.1,1.2),(1.2,1.3),(1.3,1.4),(1.4,1.5),(1.5,1.6)";
//std::string xMin_xMax_string_w2 = "(0.6,0.7),(0.7,0.8),(0.8,0.9),(0.9,1.0),(1.0,1.1),(1.1,1.2)";

//// *** x_expected*******************************************
//std::string xMin_xMax_string_x_exp = "(-1.5,-1),(-1.25,-0.75),(-1,-0.5),(-0.75,-0.25),(-0.5,0),(-0.25,0.25),(0,0.5),(0.25,0.75),(0.5,1),(0.75,1.25)";
//std::string xMin_xMax_string_x_exp = "(-0.8,0.0),(-0.9,0.1),(-1.0,0.2),(-1.1,0.3 )(-1.2,0.4) ,(-1.3,0.5),(-1.4,0.6)";
std::string xMin_xMax_string_x_exp = "(-1.5,-1.3),(-1.3,-1.1),(-1.1,-0.9),(-0.9,-0.7 ),(-0.7,-0.5) ,(-0.5,-0.3),(-0.3,-0.1),(-0.1,0.1),(0.1,0.3),(0.3,0.5),(0.5,0.7)";
///****** nsigx_fid *******************************************
//std::string xMin_xMax_string_nsigx_fid = "(0,1),(0.5,1.5),(1,2),(1.5,2.5),(2,3),(2.5,3.5),(3.5,4.5),(4,5),(5,6),(5.5,6.5),(6,7),(6.5,7.5)";
//std::string xMin_xMax_string_nsigx_fid = "(0,7),(0.5,7),(1,7),(1.5,7),(2,7),(2.5,7),(3,7),(3.5,7),(4,7),(4.5,7),(5,7),(5.5,7),(6,7)";
std::string xMin_xMax_string_nsigx_fid = "(0,0.5),(0.5,1),(1,1.5),(1.5,2),(2,2.5),(2.5,3),(3,3.5),(3.5,4),(4,4.5),(4.5,5),(5,5.5),(5.5,6),(6,6.5),(6.5,7)";
///******* y expected*************************************************************
std::string xMin_xMax_string_y_exp ="(-0.8,-0.6),(-0.6,-0.4),(-0.4,-0.2),(-0.2,0),(0,0.2),(0.2,0.4),(0.4,0.6),(0.6,0.8)";
//std::string xMin_xMax_string_y_exp = "(-0.35,0.4)";
//std::string xMin_xMax_string_y_exp = "(-0.6,0.65),(-0.55,0.6),(-0.5,0.55),(-0.45,0.5),(-0.4,0.45),(-0.35,0.4),(-0.3,0.35),(-0.25,0.3),(-0.2,0.25),(-0.15,0.2),(-0.1,0.15)";



//********* nsigy_fid ***************************************************************
std::string xMin_xMax_string_nsigy_fid = "(0,0.5),(0.5,1),(1.5,2),(2,2.5),(2.5,3)";

//********* W2 Range Study ***************************************************************
double left_w2=0.66;
double right_w2 = 1.10;
double plus_or_minus_w2 = 0.1;
double stepsize_w2 = 0.01;
int nDecimals_w2 = 2;
//********* vz Range Study ***************************************************************
double left_vz=-0.06;
double right_vz = 0.06;
double plus_or_minus_vz = 0.01;
double stepsize_vz = 0.001;
int nDecimals_vz = 3;
//********* ps Range Study ***************************************************************
double left_ps=0.2;
double right_ps = 3;
double plus_or_minus_ps = 0.1;
double stepsize_ps = 0.01;
int nDecimals_ps = 2;
//********* hcal_e Range Study ***************************************************************
double left_hcal_e=0.025;
double right_hcal_e = 1;
double plus_or_minus_hcal_e = 0.01;
double stepsize_hcal_e = 0.001;
int nDecimals_hcal_e = 3;
//********* hcal_dy Range Study ***************************************************************
double left_dy=-0.32;
double right_dy = 0.32;
double plus_or_minus_dy = 0.2;
double stepsize_dy = 0.04;
int nDecimals_dy = 2;
//*********nsigx_fid  Range Study ***************************************************************
double left_nsigx_fid=2;
double right_nsigx_fid= 7;
double plus_or_minus_nsigx_fid = 1;
double stepsize_nsigx_fid = 0.1;
int nDecimals_nsigx_fid = 2;
//*********x_exp  Range Study ***************************************************************
double left_x_exp=-1.3;
double right_x_exp= 0.5;
double plus_or_minus_x_exp = 0.1;
double stepsize_x_exp = 0.01;
int nDecimals_x_exp = 2;
//*********nsigy_fid  Range Study ***************************************************************
double left_nsigy_fid=1.5;
double right_nsigy_fid= 3;
double plus_or_minus_nsigy_fid = 1;
double stepsize_nsigy_fid = 0.1;
int nDecimals_nsigy_fid = 2;

//*********y_exp  Range Study ***************************************************************
double left_y_exp=-0.6;
double right_y_exp= 0.6;
double plus_or_minus_y_exp = 0.1;
double stepsize_y_exp = 0.05;
int nDecimals_y_exp = 2;



// SBS4 30p 
// TString DataFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_2Dhistos_sept26.root";
// TString ProtonFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_2Dhistos_sept26.root";
// TString NeutronFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_2Dhistos_sept26.root";
// double xmin_sbs4_30p = -2.15; // 
// double xmax_sbs4_30p = 1.4;//
// int PolyOrder_sbs4_30p = 2;


// SBS4 50p 
// TString DataFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_2Dhistos.root";
// TString ProtonFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deep_2Dhistos.root";
// TString NeutronFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deen_2Dhistos.root";
// double xmin_sbs4_50p = -2.5; // 
// double xmax_sbs4_50p = 1.4;//
// int PolyOrder_sbs4_50p = 2;

// SBS8 70p
// TString DataFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_2Dhistos_sept11.root";
// TString ProtonFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deep_2Dhistos.root";
// TString NeutronFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deen_2Dhistos.root";
// double xmin_sbs8_70p = -2.5; // need to adjust range
// double xmax_sbs8_70p = 1.4;//
// int PolyOrder_sbs8_70p = 4;

// SBS9 70p
// TString DataFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_2Dhistos_sept13.root";
// TString ProtonFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_deep_2Dhistos.root";
// TString NeutronFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_deen_2Dhistos.root";
// double xmin_sbs9_70p = -2.5; // need to adjust range
// double xmax_sbs9_70p = 1.4;//
// int PolyOrder_sbs9_70p = 4;




// Functions 


void stability_study_R(TString KineString="sbs4_30p", TString VarString = "dy", int PolyOrder_input=2){ // main
  // bit of a test for now. Will need to make this more sophisticated in the future. 

  gStyle->SetNumberContours(255); 
  // gStyle->SetOptStat(0110);
  gStyle->SetOptStat(0);
  
  FileNames fileNamesHandler;
  

  //std::vector<int> color_vector = {kRed,kBlue,kGreen,kMagenta,kCyan,kYellow};
  std::vector<int> color_vector = {kRed,kRed};
  

  //// Histogram to store the Rsf distriubtion
  TH1D *h_Rsf = new TH1D("h_Rsf","h_Rsf",300, 0.8,1.2);

  /// Histogram for the pull distribution
  TH1D *h_pull = new TH1D("h_pull","h_pull",10,-1,1);

  Utility utilityHandler; // class that gives us access to various functions to use between programs.

  
  // Taking user input on the kinematic and what variable we want to study. Going to try keeping it all here instead of config files. 

  if (VarString == "nsigdy"){
    HistogramName = "hcal_dx__hcal_nsigdy";
    AxisTitle ="hcal_nsigdy";
    xMin_xMax_string = xMin_xMax_string_nsigdy;
  }else if (VarString == "dy"){
    HistogramName = "hcal_dx__hcal_dy";
    AxisTitle ="hcal_dy";
    xMin_xMax_string = xMin_xMax_string_dy;
  }else if (VarString == "dy_range_left"){
    HistogramName = "hcal_dx__hcal_dy";
    AxisTitle ="hcal_dy";
    xMin_xMax_string  = utilityHandler.incrementRangeStringStepsize(left_dy, right_dy,plus_or_minus_dy, 0,stepsize_dy, nDecimals_dy);
  }else if (VarString == "dy_range_right"){
    HistogramName = "hcal_dx__hcal_dy";
    AxisTitle ="hcal_dy";
    xMin_xMax_string  = utilityHandler.incrementRangeStringStepsize(left_dy, right_dy,0,plus_or_minus_dy,stepsize_dy, nDecimals_dy);
  }else if (VarString == "vz"){
    HistogramName = "hcal_dx__tr_vz";
    AxisTitle ="track vz";
    xMin_xMax_string = xMin_xMax_string_vz;
  }else if (VarString == "vz_range_left"){
    HistogramName = "hcal_dx__tr_vz";
    AxisTitle ="track vz";
    xMin_xMax_string  = utilityHandler.incrementRangeStringStepsize(left_vz, right_vz,plus_or_minus_vz, 0,stepsize_vz, nDecimals_vz);
  }else if (VarString == "vz_range_right"){
    HistogramName = "hcal_dx__tr_vz";
    AxisTitle ="track vz";
    xMin_xMax_string  = utilityHandler.incrementRangeStringStepsize(left_vz, right_vz,0, plus_or_minus_vz,stepsize_vz, nDecimals_vz);
  } else if (VarString == "ps_range"){
    HistogramName = "hcal_dx__ps_e";
    AxisTitle ="preshower energy";
    xMin_xMax_string  = utilityHandler.incrementRangeStringStepsize(left_ps, right_ps,plus_or_minus_ps, 0,stepsize_ps, nDecimals_ps);
  }else if (VarString == "ps"){
    HistogramName = "hcal_dx__ps_e";
    AxisTitle ="preshower energy";
    xMin_xMax_string = xMin_xMax_string_ps;
  } else if (VarString == "hcal_e"){
    HistogramName = "hcal_dx__hcal_e";
    AxisTitle ="hcal energy";
    xMin_xMax_string = xMin_xMax_string_hcal_e;
  }else if (VarString == "hcal_e_range"){
    HistogramName = "hcal_dx__hcal_e";
    AxisTitle ="hcal energy";
    xMin_xMax_string  = utilityHandler.incrementRangeStringStepsize(left_hcal_e, right_hcal_e,plus_or_minus_hcal_e, 0,stepsize_hcal_e, nDecimals_hcal_e);
  } else if (VarString == "w2"){
    HistogramName = "hcal_dx__W2";
    AxisTitle ="W^{2}";
    xMin_xMax_string = xMin_xMax_string_w2;
  }else if (VarString == "w2_range_left"){
    HistogramName = "hcal_dx__W2";
    AxisTitle ="W^{2}";
    xMin_xMax_string  = utilityHandler.incrementRangeStringStepsize(left_w2, right_w2,plus_or_minus_w2, 0,stepsize_w2, nDecimals_w2);
  }else if (VarString == "w2_range_right"){
    HistogramName = "hcal_dx__W2";
    AxisTitle ="W^{2}";
    xMin_xMax_string = utilityHandler.incrementRangeStringStepsize(left_w2,right_w2,0, plus_or_minus_w2,stepsize_w2, nDecimals_w2);
  }// else if (VarString == "coin"){
  //   HistogramName = "hcal_dx__hcal_sh_atime_diff";
  //   AxisTitle ="hcal - sh time ";
  //   xMin_xMax_string = xMin_xMax_string_coin;
  // }
  // else if (VarString == "coin_range_left"){
  //   HistogramName = "hcal_dx__hcal_sh_atime_diff";
  //   AxisTitle ="hcal - sh time }";
  //   xMin_xMax_string  = utilityHandler.incrementRangeStringStepsize(left_w2, right_w2,plus_or_minus_w2, 0,stepsize_w2, nDecimals_w2);
  // }
  // else if (VarString == "coin_range_right"){
  //   HistogramName = "hcal_dx__hcal_sh_atime_diff";
  //   AxisTitle ="hcal - sh time";
  //   xMin_xMax_string = utilityHandler.incrementRangeStringStepsize(left_w2,right_w2,0, plus_or_minus_w2,stepsize_w2, nDecimals_w2);
  // }
  else if (VarString == "x_exp"){
    HistogramName = "hcal_dx__hcal_x_exp";
    AxisTitle ="x expected";
    xMin_xMax_string = xMin_xMax_string_x_exp;
  }else if (VarString == "x_exp_range_left"){
    HistogramName = "hcal_dx__hcal_x_exp";
    AxisTitle ="x expected";
    xMin_xMax_string  = utilityHandler.incrementRangeStringStepsize(left_x_exp, right_x_exp,plus_or_minus_x_exp, 0,stepsize_x_exp, nDecimals_x_exp);
  }else if (VarString == "x_exp_range_right"){
    HistogramName = "hcal_dx__hcal_x_exp";
    AxisTitle ="x expected";
    xMin_xMax_string  = utilityHandler.incrementRangeStringStepsize(left_x_exp, right_x_exp,0,plus_or_minus_x_exp,stepsize_x_exp, nDecimals_x_exp);
  }else if (VarString == "nsigx_fid"){
    HistogramName = "hcal_dx__nsigx_fid";
    AxisTitle ="nsigx_fid";
    xMin_xMax_string = xMin_xMax_string_nsigx_fid;
  }else if (VarString == "nsigx_fid_range"){
    HistogramName = "hcal_dx__nsigx_fid";
    AxisTitle ="nsigx_fid";
    xMin_xMax_string  = utilityHandler.incrementRangeStringStepsize(left_nsigx_fid, right_nsigx_fid,plus_or_minus_nsigx_fid, 0,stepsize_nsigx_fid, nDecimals_nsigx_fid);
  } else if (VarString == "y_exp"){
    HistogramName = "hcal_dx__hcal_y_exp";
    AxisTitle ="y expected";
    xMin_xMax_string = xMin_xMax_string_y_exp;
  }else if (VarString == "y_exp_range_left"){
    HistogramName = "hcal_dx__hcal_y_exp";
    AxisTitle ="y expected";
    xMin_xMax_string  = utilityHandler.incrementRangeStringStepsize(left_y_exp, right_y_exp,plus_or_minus_y_exp, 0,stepsize_y_exp, nDecimals_y_exp);
  }else if (VarString == "y_exp_range_right"){
    HistogramName = "hcal_dx__hcal_y_exp";
    AxisTitle ="y expected";
    xMin_xMax_string  = utilityHandler.incrementRangeStringStepsize(left_y_exp, right_y_exp,0,plus_or_minus_y_exp,stepsize_y_exp, nDecimals_y_exp);
  }  else if (VarString == "nsigy_fid"){
    HistogramName = "hcal_dx__nsigy_fid";
    AxisTitle ="nsigy_fid";
    xMin_xMax_string = xMin_xMax_string_nsigy_fid;
  }else if (VarString == "nsigy_fid"){
    HistogramName = "hcal_dx__nsigy_fid";
    AxisTitle ="nsigy_fid";
    xMin_xMax_string = xMin_xMax_string_nsigy_fid;
  }else if (VarString == "nsigy_fid_range"){
    HistogramName = "hcal_dx__nsigy_fid";
    AxisTitle ="nsigy_fid";
    xMin_xMax_string  = utilityHandler.incrementRangeStringStepsize(left_nsigy_fid, right_nsigy_fid,plus_or_minus_nsigy_fid, 0,stepsize_nsigy_fid, nDecimals_nsigy_fid);
  } else {
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

  

  // // testing incrementRange here quick for vz:
  // double min=-0.055;
  // double max = 0.055;
  // double plus_or_minus_xmin = 0.01;
  // double plus_or_minus_xmax = 0;
  // double stepsize = 0.001;
  // int nDecimals = 3;
 

  //  xMin_xMax_string = utilityHandler.incrementRangeStringStepsize(min, max,plus_or_minus_xmin, plus_or_minus_xmax,stepsize, nDecimals);
      
  // cout<<"Testing writing as a string"<<endl;
  // cout<<xMin_xMax_string<<endl;
  // cout<<endl;



  
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
 

    
  std::vector<FitHistogram> fitHandler_vector; // vector of objects fo the FitHistogram class that will help us do data-mc comparion. . 
    

   

  cout<<endl<<"Slicing Data"<<endl;
  utilityHandler.SliceAndProjectHistogram_xMinxMax(hist_2D_data, xMinSlices,xMaxSlices, hist_vector_data, AxisTitle,"dx","data");
  cout<<endl<<"Slicing Proton"<<endl;
  utilityHandler.SliceAndProjectHistogram_xMinxMax(hist_2D_proton, xMinSlices, xMaxSlices, hist_vector_p, AxisTitle,"dx","p");
  cout<<endl<<"Slicing Neutron"<<endl;
  utilityHandler.SliceAndProjectHistogram_xMinxMax(hist_2D_neutron, xMinSlices, xMaxSlices, hist_vector_n, AxisTitle,"dx","n");
  cout<<endl;
  


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
  // initialParameters = {0};
      
  /// Fit each slice 
  for (int sliceid = 0 ; sliceid < xMinSlices.size() ; sliceid++)
    {
      // this uses clones of the histograms
      TH1D *hist_p  =(TH1D*)hist_vector_p[sliceid]->Clone("hist_p");
      TH1D *hist_n  =(TH1D*)hist_vector_n[sliceid]->Clone("hist_n");
      TH1D *hist_data_dont_draw  =(TH1D*)hist_vector_data[sliceid]->Clone("hist_data_dont_draw");
      TH1D *hist_data  =(TH1D*)hist_vector_data[sliceid]->Clone("hist_data");

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
      // poly_result_vector_of_arrays.push_back(poly_result);
      poly_result_vector_of_vectors.push_back(poly_result);
      poly_result_err_vector_of_vectors.push_back(poly_result_err);
	  
      // delete Fit; //
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
      hist_temp_p = utilityHandler.ScaleAndShiftHistogram(hist_vector_p[sliceid], scale_p_vector[sliceid], shift_p_vector[sliceid]);
      hist_temp_n = utilityHandler.ScaleAndShiftHistogram(hist_vector_n[sliceid], scale_n_vector[sliceid], shift_n_vector[sliceid]);
	 
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
  // utilityHandler.customizeGraphMore(Rsf_graph, 33, kBlue, 3,"",AxisTitle,"Rsf");
  utilityHandler.customizeGraphMore(Rsf_graph, 33, kBlue, 3,"",AxisTitle,"Rsf",1.4,1.4);

  // Fit a straight line to the graph
  TF1* fit_Rsf_graph = new TF1("fit_Rsf_graph", "[0]",Rsf_graph ->GetX()[0], Rsf_graph->GetX()[Rsf_graph->GetN()-1]);
  fit_Rsf_graph->SetLineColor(kRed);
  Rsf_graph->Fit(fit_Rsf_graph, "Q RN"); // N supresses the drawing of it automatically

    // Get fit parameters
  double constant = fit_Rsf_graph->GetParameter(0);       // The constant value
  double constantError = fit_Rsf_graph->GetParError(0);    // The error on the constant
  double chi2 = fit_Rsf_graph->GetChisquare();             // The chi-squared value
  int ndf = fit_Rsf_graph->GetNDF();                       // The number of degrees of freedom
  double chi2_ndf = chi2 / ndf;                          // chi2/ndf

  
  // Fit a pol1 to the graph 
  TF1* fit_pol1_Rsf_graph = new TF1("fit_pol1_Rsf_graph", "pol1",Rsf_graph ->GetX()[0], Rsf_graph->GetX()[Rsf_graph->GetN()-1]);
  fit_pol1_Rsf_graph->SetLineColor(kViolet);
  fit_pol1_Rsf_graph->SetLineWidth(2);
  fit_pol1_Rsf_graph->SetLineStyle(9); 
  Rsf_graph->Fit(fit_pol1_Rsf_graph, "Q RN"); // N supresses the drawing of it automatically

  
  //// canvas
  TCanvas *graphcanvas = new TCanvas("graphcanvas","graphcanvas",1200,600);
  graphcanvas->SetGrid();
  graphcanvas->Divide(2,1);
  graphcanvas->cd(1);
  utilityHandler.adjustCanvas(graphcanvas);
  Rsf_graph->Draw("AP");
  fit_Rsf_graph->Draw("same");
  fit_pol1_Rsf_graph->Draw("same");
  graphcanvas->Update();
  

  // Use TLatex to add the fit result and chi²/ndf to the canvas
  TLatex latex;
  latex.SetNDC();  // Use normalized coordinates (0 to 1)
  latex.SetTextSize(0.04);
  latex.DrawLatex(0.15, 0.85, Form("pol0: y = %.5f #pm %.5f, #chi^{2}/ndf = %.2f", constant, constantError,chi2_ndf));
  // latex.DrawLatex(0.15, 0.80, Form("#chi^{2}/ndf = %.2f/%i = %.2f", chi2,ndf,chi2_ndf));
  latex.DrawLatex(0.15, 0.2, Form("Rsf Mean = %.5f #pm StDev = %.5f", Rsf_mean, Rsf_stdev));
  
  graphcanvas->Update();
  //graphcanvas->SaveAs(Form("%s/Rsf_xCenter.pdf",outputfilelocation.Data() ) );

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
  utilityHandler.customizeGraphMore(Chi2_ndf_graph, 33, kBlue, 3,"","Bin Width","Chi^{2}/ndf",1.4,1.4);

  graphcanvas->cd(2);
  Chi2_ndf_graph ->Draw("AP");
  graphcanvas->Update();
  graphcanvas->SaveAs(Form("%s/Rsf.pdf",outputfilelocation.Data() ) );

  /// /// canvas
  // TCanvas *chi2_ndf_canvas = new TCanvas("chi2_ndf_canvas","chi2_ndf_canvas",800,600);  
  // utilityHandler.adjustCanvas(chi2_ndf_canvas);
  // Chi2_ndf_graph ->Draw("AP");
  // chi2_ndf_canvas->Update();
  // chi2_ndf_canvas->SaveAs(Form("%s/Chi2_ndf.pdf",outputfilelocation.Data()));
  
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
  utilityHandler.customizeGraphMore(Rsf_xMax_graph, 33, kBlue, 3,"",AxisTitle,"Rsf",1.4,1.4);
  
  // Fit a straight line to the graph
  TF1* fit_Rsf_xMax_graph = new TF1("fit_Rsf_xMax_graph", "[0]",Rsf_xMax_graph ->GetX()[0], Rsf_xMax_graph->GetX()[Rsf_xMax_graph->GetN()-1]);
  Rsf_xMax_graph->Fit(fit_Rsf_xMax_graph, "Q R");
  
  //// canvas
  TCanvas *Rsf_xMaxCanvas = new TCanvas("Rsf_xMaxCanvas","Rsf_xMaxCanvas",1200,600);
  Rsf_xMaxCanvas->SetGrid();
  Rsf_xMaxCanvas ->Divide(2,1);
  Rsf_xMaxCanvas ->cd(1);
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
  latex_xMax.DrawLatex(0.15, 0.85, Form("pol0: y = %.5f #pm %.5f, #chi^{2}/ndf = %.2f", constant_xMax, constantError_xMax,chi2_ndf_xMax ));
  // latex_xMax.DrawLatex(0.2, 0.80, Form("#chi^{2}/ndf = %.2f/%i = %.2f", chi2_xMax,ndf_xMax,chi2_ndf_xMax));
  latex_xMax.DrawLatex(0.15, 0.2, Form("Rsf Mean = %.5f #pm StDev = %.5f", Rsf_mean, Rsf_stdev));

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
  latex_xMin.DrawLatex(0.15, 0.85, Form("pol0: y = %.5f #pm %.5f, #chi^{2}/ndf = %.2f", constant_xMin, constantError_xMin,chi2_ndf_xMin));
  // latex_xMin.DrawLatex(0.2, 0.80, Form("#chi^{2}/ndf = %.2f/%i = %.2f", chi2_xMin,ndf_xMin,chi2_ndf_xMin));
  latex_xMin.DrawLatex(0.15, 0.2, Form("Rsf Mean = %.5f #pm StDev = %.5f", Rsf_mean, Rsf_stdev));
  

  int nHist = hist_vector_data.size();
  int nCols = 3;
  int nRows = (nHist + nCols - 1) / nCols;


  std::vector<TH1D*> overall_fit_as_histogram;

  TCanvas* fits_canvas = new TCanvas("fits_canvas", "fits_canvas", 1100, 700);
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
    hist_result_vector_p[i] ->SetLineColor(kGreen+2);
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
  

  //   // // TCanvas* testcanvas1  = new TCanvas("testcanvas1", "testcanvas1", 800, 600);
  //   // // TH1D* sumHistoTest = utilityHandler.sumHistogramsWithPolynomial(hist_result_vector_p[0],hist_result_vector_n[0] , poly_fit_result[0]);
  //   // // hist_vector_data[0]->Draw();
  //   // // sumHistoTest ->Draw("same");

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
  std::string proton_title = hist_2D_proton->GetTitle();
  cout<<"proton title: "<<endl;
  cout<< proton_title<<endl;
  std::string neutron_title =hist_2D_neutron->GetTitle();
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


