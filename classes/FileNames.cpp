#include "FileNames.h"
#include <iostream>

// Constructor
FileNames::FileNames() {
  
  // SBS4 0p LD2
  //DataFileString_sbs4_0p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/";
  // ProtonFileString_sbs4_0p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_2Dhistos_sept26.root";
  //NeutronFileString_sbs4_0p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_2Dhistos_sept26.root";
  //InelasticFileString_sbs4_0p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_inel_2Dhistos.root";

  // SBS4 0p LH2
  DataFileString_sbs4_0p_LH2 = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_0p_cuts_LH2_2Dhistos_oct14.root";// oct 14 is using average mass 

  // SBS4 0p  Fit Range (probably could be adjusted
  xmin_sbs4_0p = -2.15;
  xmax_sbs4_0p = 1.4;
  
  /// SBS4 30p LD2
  DataFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_2Dhistos_oct10.root";
  ProtonFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_2Dhistos_oct10.root";
  NeutronFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_2Dhistos_oct10.root";
  InelasticFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_inel_2Dhistos_oct14.root";

  /// SBS4 30p LD2 BinStudy
  DataFileString_sbs4_30p_BinStudy = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_BinStudy_oct4.root";
  ProtonFileString_sbs4_30p_BinStudy = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_BinStudy_oct4.root";
  NeutronFileString_sbs4_30p_BinStudy = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_BinStudy_oct4.root";

  // SBS4 30p LH2
  DataFileString_sbs4_30p_LH2="/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_LH2_2Dhistos_sept26.root";

  // SBS4 30p Fit Range 
  xmin_sbs4_30p = -2.15;
  xmax_sbs4_30p = 1.4;

  //SBS4 30p initial fit params 
  initialParams_sbs4_30p ={1,1,-0.05,-0.042,1,1,-1};

  // SBS4 50p LD2
  DataFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_2Dhistos_sept26.root";
  ProtonFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deep_2Dhistos_sept26.root";
  NeutronFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deen_2Dhistos_sept26.root";
  InelasticFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_inel_2Dhistos.root";

  /// SBS4 50p LD2 BinStudy
  DataFileString_sbs4_50p_BinStudy = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_BinStudy_oct8.root";
  ProtonFileString_sbs4_50p_BinStudy = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deep_BinStudy_oct8.root";
  NeutronFileString_sbs4_50p_BinStudy = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deen_BinStudy_oct8.root";

  // SBS4 50p LH2
  DataFileString_sbs4_50p_LH2= "/w/halla-scshelf2102/sbs/msatnik/GMn/outputsbs4_50p_cuts_LH2_2Dhistos_sept26.root";

  // SBS4 50p Fit Range 
  xmin_sbs4_50p = -2.63;
  xmax_sbs4_50p = 1.4;

  // SBS4 50p initial fit params
  initialParams_sbs4_50p ={1,1,-0.05,-0.03,1,1,-1};

  // SBS8 70p LD2
  DataFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_2Dhistos_sept26.root";
  ProtonFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deep_2Dhistos_sept26.root";
  NeutronFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deen_2Dhistos_sept26.root";
  InelasticFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_inel_2Dhistos.root";

  // SBS8 70p LD2 BinStudy
  DataFileString_sbs8_70p_BinStudy = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_BinStudy_oct8.root" ; 
  ProtonFileString_sbs8_70p_BinStudy = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deep_BinStudy_oct8.root" ; 
  NeutronFileString_sbs8_70p_BinStudy=  "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deen_BinStudy_oct8.root";

  // SBS8 70p LH2
  DataFileString_sbs8_70p_LH2 = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_LH2_2Dhistos_sept26.root";

  // SBS8 70p Fit Range 
  xmin_sbs8_70p = -1.9;
  xmax_sbs8_70p = 1.0;

  // SBS8 70p initial fit params 
  initialParams_sbs8_70p= {1,1,-0.04,0.01,1,1,-1};
  
  // SBS8 50p LD2
  DataFileString_sbs8_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_50p_cuts_2Dhistos_sept26.root";
  ProtonFileString_sbs8_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_50p_cuts_deep_2Dhistos_sept26.root";
  NeutronFileString_sbs8_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_50p_cuts_deen_2Dhistos_sept26.root";
  InelasticFileString_sbs8_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_inel_2Dhistos.root"; // NEED TO RUN

  // SBS8 50p LD2 BinStudy
  DataFileString_sbs8_50p_BinStudy = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_50p_cuts_BinStudy_oct8.root" ; 
  ProtonFileString_sbs8_50p_BinStudy = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_50p_cuts_deep_BinStudy_oct8.root" ; 
  NeutronFileString_sbs8_50p_BinStudy=  "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_50p_cuts_deen_BinStudy_oct8.root";

  // SBS8 50p LH2
  DataFileString_sbs8_50p_LH2 = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_50p_cuts_LH2_2Dhistos_sept26.root";

  // SBS8 50p Fit Range 
  xmin_sbs8_50p = -1.6;
  xmax_sbs8_50p = 1.0;

  // SBS8 50p initial fit params 
  initialParams_sbs8_50p= {1,1,-0.03,0.005,1,1,-1};

  // SBS8 100p LD2
  // DataFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_2Dhistos_sept24.root";
  // ProtonFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_deep_2Dhistos_sept24.root";
  // NeutronFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_deen_2Dhistos_sept24.root";
  // InelasticFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_inel_2Dhistos.root"; // NEED TO RUN
  
  DataFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_2Dhistos_oct8.root";
  ProtonFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_deep_2Dhistos_oct8.root";
  NeutronFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_deen_2Dhistos_oct8.root";
  InelasticFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_inel_2Dhistos.root"; // NEED TO RUN

  // SBS8 100p LD2 BinStudy
  DataFileString_sbs8_100p_BinStudy = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_BinStudy_oct8.root" ; 
  ProtonFileString_sbs8_100p_BinStudy = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_deep_BinStudy_oct8.root" ; 
  NeutronFileString_sbs8_100p_BinStudy=  "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_deen_BinStudy_oct8.root";

  // SBS8 100p LH2
  DataFileString_sbs8_100p_LH2= "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_LH2_2Dhistos_sept26.root";

  // SBS8 100p Fit Range 
  xmin_sbs8_100p = -2.2;
  xmax_sbs8_100p = 1.0;

  // SBS8 100p initial fit params 
  initialParams_sbs8_100p ={1,1,-0.04,0.01,1,1,-1};

  // SBS9 70p LD2
  DataFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_2Dhistos_sept26_1.root";
  ProtonFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_deep_2Dhistos_sept26_z_1.root";
  NeutronFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_deen_2Dhistos_sept26_z_1.root";
  InelasticFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_inel_2Dhistos.root";

  // SBS9 70p LD2 BinStudy
  DataFileString_sbs9_70p_BinStudy = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_BinStudy_oct8.root" ; 
  ProtonFileString_sbs9_70p_BinStudy = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_deep_BinStudy_oct8.root" ; 
  NeutronFileString_sbs9_70p_BinStudy=  "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_deen_BinStudy_oct8.root";

  // SBS 9 70p LH2
  DataFileString_sbs9_70p_LH2="/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_LH2_2Dhistos.root"; /// Not sure why I don't have a date on this one.. may want to rerun.

  // SBS9 70p Fit Range 
  xmin_sbs9_70p = -1.9;
  xmax_sbs9_70p = 1.0;

  // SBS9 70p initial fit params 
  initialParams_sbs9_70p= {1,1,-0.05,-0.025,1,1,-1};
  
}// end constructor 

// Destructor
FileNames::~FileNames() {
  // Destructor body (can be empty)
}
