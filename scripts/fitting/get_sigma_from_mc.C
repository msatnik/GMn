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


// global params

// functions 


//
/// Porogram to just read in a histogram and print out it's mean, mean error, stdev, and stdev error
//
void get_sigma_from_mc(){//main
  gStyle->SetOptFit(11111);
  gStyle->SetCanvasPreferGL(1);
  gStyle -> SetOptStat(0);
  gStyle ->SetEndErrorSize(0);



  TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deep_2Dhistos.root"); // 
  

  TH1D *hist1 = (TH1D*)f1->Get("hcal_dx_1d_allcuts");
  hist1->GetXaxis()->SetRangeUser(-1.8,0.2);
  

  cout<< "Mean = "<<hist1->GetMean()<<" +/- "<<hist1->GetMeanError()<<", StdDev = "<<hist1->GetStdDev()<<" +/- "<<hist1->GetStdDevError()<<endl;

  hist1->Draw();
  

}//end main
