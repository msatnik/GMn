#ifndef FILENAMES_H
#define FILENAMES_H

class FileNames {
private:

public:

  // SBS4 0p LD2
  //TString DataFileString_sbs4_0p;
  // TString ProtonFileString_sbs4_0p;
  //TString NeutronFileString_sbs4_0p;
  //TString InelasticFileString_sbs4_0p;
  
  // SBS4 0p LH2
  TString DataFileString_sbs4_0p_LH2;

  // SBS4 0p  Fit Range (probably could be adjusted
  double xmin_sbs4_0p;
  double xmax_sbs4_0p;
  
  /// SBS4 30p LD2
  TString DataFileString_sbs4_30p; 
  TString ProtonFileString_sbs4_30p; 
  TString NeutronFileString_sbs4_30p;
  TString InelasticFileString_sbs4_30p;

  // SBS4 30p LD2 BinStudy
  TString DataFileString_sbs4_30p_BinStudy; 
  TString ProtonFileString_sbs4_30p_BinStudy; 
  TString NeutronFileString_sbs4_30p_BinStudy;
  
  // SBS4 30p LH2
  TString DataFileString_sbs4_30p_LH2;

  // SBS4 30p Fit Range 
  double xmin_sbs4_30p;
  double xmax_sbs4_30p;

  //SBS4 30p initial fit params 
  std::vector<double> initialParams_sbs4_30p;

  // SBS4 50p LD2
  TString  DataFileString_sbs4_50p; 
  TString  ProtonFileString_sbs4_50p;
  TString NeutronFileString_sbs4_50p; 
  TString InelasticFileString_sbs4_50p;

  // SBS4 50p LD2 BinStudy
  TString DataFileString_sbs4_50p_BinStudy; 
  TString ProtonFileString_sbs4_50p_BinStudy; 
  TString NeutronFileString_sbs4_50p_BinStudy;

  // SBS4 50p LH2
  TString DataFileString_sbs4_50p_LH2;

  // SBS4 50p Fit Range 
  double xmin_sbs4_50p;
  double xmax_sbs4_50p;

  //SBS4 50p initial fit params 
  std::vector<double> initialParams_sbs4_50p;
  
  // SBS8 70p LD2
  TString DataFileString_sbs8_70p;
  TString ProtonFileString_sbs8_70p; 
  TString NeutronFileString_sbs8_70p; 
  TString InelasticFileString_sbs8_70p;

  // SBS8 70p LD2 BinStudy
  TString DataFileString_sbs8_70p_BinStudy; 
  TString ProtonFileString_sbs8_70p_BinStudy; 
  TString NeutronFileString_sbs8_70p_BinStudy;

  // SBS8 70p LH2
  TString DataFileString_sbs8_70p_LH2;

  // SBS8 70p Fit Range 
  double xmin_sbs8_70p;
  double xmax_sbs8_70p;
 
  //SBS8 70p initial fit params 
  std::vector<double> initialParams_sbs8_70p;
  
  // SBS8 50p LD2
  TString DataFileString_sbs8_50p;
  TString ProtonFileString_sbs8_50p;  
  TString NeutronFileString_sbs8_50p; 
  TString InelasticFileString_sbs8_50p;

  // SBS8 50p LD2 BinStudy
  TString DataFileString_sbs8_50p_BinStudy; 
  TString ProtonFileString_sbs8_50p_BinStudy; 
  TString NeutronFileString_sbs8_50p_BinStudy;

  // SBS8 50p LH2
  TString DataFileString_sbs8_50p_LH2;

  // SBS8 50p Fit Range 
  double xmin_sbs8_50p;
  double xmax_sbs8_50p;

  //SBS8 50p initial fit params 
  std::vector<double> initialParams_sbs8_50p;
  
  // SBS8 100p LD2
  TString DataFileString_sbs8_100p; 
  TString ProtonFileString_sbs8_100p;  
  TString NeutronFileString_sbs8_100p; 
  TString InelasticFileString_sbs8_100p;


  // SBS8 100p LD2 BinStudy
  TString DataFileString_sbs8_100p_BinStudy; 
  TString ProtonFileString_sbs8_100p_BinStudy; 
  TString NeutronFileString_sbs8_100p_BinStudy;

  // SBS8 100p LH2
  TString DataFileString_sbs8_100p_LH2;

  // SBS8 100p Fit Range 
  double xmin_sbs8_100p;
  double xmax_sbs8_100p;

  //SBS8 50p initial fit params 
  std::vector<double> initialParams_sbs8_100p;
    
  // SBS9 70p LD2
  TString DataFileString_sbs9_70p;
  TString ProtonFileString_sbs9_70p; 
  TString NeutronFileString_sbs9_70p;   
  TString InelasticFileString_sbs9_70p;

  // SBS9 70p LD2 BinStudy
  TString DataFileString_sbs9_70p_BinStudy; 
  TString ProtonFileString_sbs9_70p_BinStudy; 
  TString NeutronFileString_sbs9_70p_BinStudy;

  // SBS9 70p LH2
  TString DataFileString_sbs9_70p_LH2;

  // SBS9 70p Fit Range 
  double xmin_sbs9_70p;
  double xmax_sbs9_70p;

  //SBS9 70p initial fit params 
  std::vector<double> initialParams_sbs9_70p;
  
  FileNames();  // Constructor
  ~FileNames();  // Destructor
  
};
#endif
