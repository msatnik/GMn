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
using namespace std;

void MakeBinStudyHistos_mc(TString configfileinput="sbs4_30p_test_mc"){ // main

  gStyle->SetNumberContours(255); 

  cout<<"Loading branches. Hold tight! "<<endl;
  

  TChain *C = new TChain("P");

  //string configfilename = Form("../config/sbs%d.cfg",kine);
  TString configfilename = "../config/" + configfileinput + ".cfg";
  // string configfilename = "/w/halla-scshelf2102/sbs/msatnik/GMn/config/sbs4_30p.cfg";
  cout<<"reading from config file: "<<configfilename<<endl;


  // set location and name for output file 
  TString outputfilename = "../output/" + configfileinput + "_BinStudy.root";
  //TString outputfilename = "../output/test.root";

  // Declare outfile
  // TFile *fout = new TFile( Form("../output/sbs%d.root",kine), "RECREATE" );
  //TFile *fout = new TFile("../output/sbs4_30p.root", "RECREATE" );
  TFile *fout = new TFile(outputfilename,"RECREATE");
  cout<<"writing to file: "<< outputfilename <<endl;

  // C->Add("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_may20.root");
  // C->Add("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_may18.root");
  //C->Add("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_simc_deen.root");
  // C->Add("/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_simc_deep.root");


  double E_e=0;

  // Setting up cuts 
  std::string EnergyCutString = "bb_ps_e>0.2&&(bb_ps_e+bb_sh_e)>1.7"; // PS and Shower Energy cut
  std::string TrackQualityCutString = "bb_tr_n==1&&bb_gem_track_nhits>=3"; //  coudl also put chi2 cut
  std::string TargetVertexCutString = "abs(bb_tr_vz)<0.07";
  // eventually need optics cuts
  std::string W2CutString = "W2<2"; 
  std::string FidXCutString = "nsigx_fid>0.5";
  std::string FidYCutString = "nsigy_fid>0.5";
  std::string dyCutString = "abs(hcal_dy)<3*dysig";
  std::string e_over_p_CutString = "abs(e_over_p - 1)<0.5";
  std::string HCal_Shower_atime_CutString= "abs(hcal_sh_atime_diff)<10";
  std::string HCal_Energy_CutString = "hcal_e>0.01";
  std::string Optics_CutString = "abs(bb_tr_r_x-bb_tr_r_th*0.9)<0.3";
  // data vs sim
  std::string is_data_CutString = "is_hydrogen==1||is_deuterium==1";
  std::string is_simulation_CutString = "is_proton==1||is_neutron==1";
  std::string is_proton_CutString= "is_proton==1";
  std::string is_neutron_CutString= "is_neutron==1";


  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      if(!currentline) cout << "WARNING: No file exists at " << currentline << "." << endl;
      C->Add(currentline);
      cout << "Loaded file at: " << currentline << endl;
    }    
  }
  TCut globalcut = "";
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
    }    
    cout<< "Global Cut: "<<globalcut<<endl;
  }

  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
    TObjArray *tokens = currentline.Tokenize(" ");
    Int_t ntokens = tokens->GetEntries();
    if( ntokens>1 ){
      TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
      if( skey == "E_e" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	E_e = sval.Atof();
	cout << "Loading beam energy: " << E_e << endl;
      }
      if( skey == "EnergyCutString" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        EnergyCutString = sval.Data();
	cout << "Loading EnergyCutString: " << EnergyCutString << endl;
      }
      if( skey == "TrackQualityCutString" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        TrackQualityCutString = sval.Data();
	cout << "Loading TrackQualityCutString: " << TrackQualityCutString << endl;
      }
      if( skey == "TargetVertexCutString" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        TargetVertexCutString = sval.Data();
	cout << "Loading TargetVertexCutString: " << TargetVertexCutString << endl;
      }
      if( skey == "W2CutString" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	W2CutString = sval.Data();
	cout << "Loading W2CutString: " <<W2CutString << endl;
      }
      if( skey == "FidXCutString" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	FidXCutString = sval.Data();
	cout << "Loading FidXCutString: " <<FidXCutString << endl;
      }
      if( skey == "FidYCutString" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	FidYCutString = sval.Data();
	cout << "Loading FidYCutString: " <<FidYCutString << endl;
      }
      if( skey == "dyCutString" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dyCutString = sval.Data();
	cout << "Loading dyCutString: " <<dyCutString << endl;
      }
      if( skey == "e_over_p_CutString" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        e_over_p_CutString = sval.Data();
	cout << "Loading e_over_p_CutString: " <<e_over_p_CutString << endl;
      }
      if( skey == "HCal_Energy_CutString" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        HCal_Energy_CutString= sval.Data();
	cout << "Loading HCal_Energy_CutString: " << HCal_Energy_CutString<< endl;
      }
      if( skey == "HCal_Shower_atime_CutString" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        HCal_Shower_atime_CutString= sval.Data();
	cout << "Loading HCal_Shower_atime_CutString: " << HCal_Shower_atime_CutString<< endl;
      }
      if( skey == "Optics_CutString" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	Optics_CutString= sval.Data();
	cout << "Loading Optics_CutString: " << Optics_CutString<< endl;
      }


    }
    delete tokens;
  }

  configfile.close();


 

  double  bb_tr_n, bb_tr_vz, bb_tr_p, bb_tr_th, bb_tr_ph; // track variables
  double bb_tr_r_x, bb_tr_r_y, bb_tr_r_th, bb_tr_r_ph;// optics variables 
  double bb_ps_e, bb_ps_rowblk, bb_ps_colblk, bb_sh_e, bb_sh_rowblk, bb_sh_colblk; //shower variables
  double bb_gem_track_nhits, bb_gem_track_ngoodhits, bb_gem_track_chi2ndf ;// track quality
  double hcal_x, hcal_y,hcal_e; // hcal 
  double  hcal_dx, hcal_dy, hcal_x_exp, hcal_y_exp;// hcal projections 
  double Q2, W2, nu, tau, epsilon; // physics quantities 
  int is_deuterium, is_hydrogen; // truth variables on what target was analyzed. 
  double dxsig_p, dxsig_n, dysig, dx_pn; // variables used in fid cuts
  double nsigx_fid, nsigy_fid; // fiducial cut size 
  double e_over_p;
  double mag_field; 

  double hcal_clus_atime, bb_sh_atimeblk, hcal_sh_atime_diff; // time for coinicidene cuts 

  // MC variables 
  int is_proton, is_neutron; // if it was proton simulaiton or neutron simulation
  double corrected_weight;

  

  // set branches 
  C->SetBranchAddress("bb_tr_n",&bb_tr_n);
  C->SetBranchAddress("bb_tr_vz",&bb_tr_vz);
  C->SetBranchAddress("bb_tr_p",&bb_tr_p); 
  C->SetBranchAddress("bb_tr_th",&bb_tr_th);
  C->SetBranchAddress("bb_tr_ph",&bb_tr_ph);
  C->SetBranchAddress("bb_tr_r_x",&bb_tr_r_x);
  C->SetBranchAddress("bb_tr_r_y",&bb_tr_r_y);
  C->SetBranchAddress("bb_tr_r_th",&bb_tr_r_th);
  C->SetBranchAddress("bb_tr_r_ph",&bb_tr_r_ph);

  C->SetBranchAddress("bb_ps_e",&bb_ps_e);
  C->SetBranchAddress("bb_ps_rowblk",&bb_ps_rowblk);
  C->SetBranchAddress("bb_ps_colblk",&bb_ps_colblk);
  C->SetBranchAddress("bb_sh_e",&bb_sh_e);
  C->SetBranchAddress("bb_sh_rowblk",&bb_sh_rowblk);
  C->SetBranchAddress("bb_sh_colblk",&bb_sh_colblk);

  C->SetBranchAddress("bb_gem_track_nhits",&bb_gem_track_nhits);
  C->SetBranchAddress("bb_gem_track_ngoodhits",&bb_gem_track_ngoodhits);
  C->SetBranchAddress("bb_gem_track_chi2ndf",&bb_gem_track_chi2ndf);

  C->SetBranchAddress("hcal_x",&hcal_x);
  C->SetBranchAddress("hcal_y",&hcal_y);
  C->SetBranchAddress("hcal_dx",&hcal_dx);
  C->SetBranchAddress("hcal_dy",&hcal_dy);
  C->SetBranchAddress("hcal_x_exp",&hcal_x_exp);
  C->SetBranchAddress("hcal_y_exp",&hcal_y_exp);
  C->SetBranchAddress("hcal_e",&hcal_e);

  C->SetBranchAddress("Q2",&Q2);
  C->SetBranchAddress("W2",&W2);
  C->SetBranchAddress("nu",&nu);
  C->SetBranchAddress("tau",&tau);
  C->SetBranchAddress("epsilon",&epsilon);

  C->SetBranchAddress("is_deuterium",&is_deuterium); // data only
  C->SetBranchAddress("is_hydrogen",&is_hydrogen);

  C->SetBranchAddress("dxsig_p",&dxsig_p);
  C->SetBranchAddress("dxsig_n",&dxsig_n);
  C->SetBranchAddress("dysig",&dysig);
  C->SetBranchAddress("dx_pn",&dx_pn);

  C->SetBranchAddress("nsigx_fid",&nsigx_fid);
  C->SetBranchAddress("nsigy_fid",&nsigy_fid);

  C->SetBranchAddress("e_over_p",&e_over_p);
  C->SetBranchAddress("mag_field",&mag_field);

  C->SetBranchAddress("hcal_clus_atime",&hcal_clus_atime);
  C->SetBranchAddress("bb_sh_atimeblk",&bb_sh_atimeblk);
  C->SetBranchAddress("hcal_sh_atime_diff",&hcal_sh_atime_diff);

  C->SetBranchAddress("is_proton",&is_proton);
  C->SetBranchAddress("is_neutron",&is_neutron);
  C->SetBranchAddress("corrected_weight",&corrected_weight);




  // make TH1D with all the cuts 
  // since this is mc data, we don't want to cut on adc time. But we need something with all the other cuts for us to compare and fit to when looking at adc time on data. 
  std::string allcuts_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" +W2CutString+ "&&"+ FidXCutString + "&&"+ FidYCutString +"&&"+e_over_p_CutString +  "&&"+ HCal_Energy_CutString+ "&&"+dyCutString;


  // Make with various numbers of bins. 
  //// Draw the 1D histogram
  C->Draw("hcal_dx>>hcal_dx_1d_allcuts_50(50, -4, 4)", Form("corrected_weight*(%s)", allcuts_study_string.c_str()), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_dx_1d_allcuts_50= (TH1D*)gDirectory->Get("hcal_dx_1d_allcuts_50");
  if (hcal_dx_1d_allcuts_50) {
    hcal_dx_1d_allcuts_50->SetXTitle("hcal_dx");
  }

  C->Draw("hcal_dx>>hcal_dx_1d_allcuts_100(100, -4, 4)", Form("corrected_weight*(%s)", allcuts_study_string.c_str()), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_dx_1d_allcuts_100= (TH1D*)gDirectory->Get("hcal_dx_1d_allcuts_100");
  if (hcal_dx_1d_allcuts_100) {
    hcal_dx_1d_allcuts_100->SetXTitle("hcal_dx");
  }
  

  C->Draw("hcal_dx>>hcal_dx_1d_allcuts_200(200, -4, 4)", Form("corrected_weight*(%s)", allcuts_study_string.c_str()), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_dx_1d_allcuts_200= (TH1D*)gDirectory->Get("hcal_dx_1d_allcuts_200");
  if (hcal_dx_1d_allcuts_200) {
    hcal_dx_1d_allcuts_200->SetXTitle("hcal_dx");
  }

   C->Draw("hcal_dx>>hcal_dx_1d_allcuts_300(300, -4, 4)", Form("corrected_weight*(%s)", allcuts_study_string.c_str()), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_dx_1d_allcuts_300= (TH1D*)gDirectory->Get("hcal_dx_1d_allcuts_300");
  if (hcal_dx_1d_allcuts_300) {
    hcal_dx_1d_allcuts_300->SetXTitle("hcal_dx");
  }

  C->Draw("hcal_dx>>hcal_dx_1d_allcuts_400(400, -4, 4)", Form("corrected_weight*(%s)", allcuts_study_string.c_str()), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_dx_1d_allcuts_400= (TH1D*)gDirectory->Get("hcal_dx_1d_allcuts_400");
  if (hcal_dx_1d_allcuts_400) {
    hcal_dx_1d_allcuts_400->SetXTitle("hcal_dx");
  }

  C->Draw("hcal_dx>>hcal_dx_1d_allcuts_500(500, -4, 4)", Form("corrected_weight*(%s)", allcuts_study_string.c_str()), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_dx_1d_allcuts_500= (TH1D*)gDirectory->Get("hcal_dx_1d_allcuts_500");
  if (hcal_dx_1d_allcuts_500) {
    hcal_dx_1d_allcuts_500->SetXTitle("hcal_dx");
  }

   C->Draw("hcal_dx>>hcal_dx_1d_allcuts_600(600, -4, 4)", Form("corrected_weight*(%s)", allcuts_study_string.c_str()), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_dx_1d_allcuts_600= (TH1D*)gDirectory->Get("hcal_dx_1d_allcuts_600");
  if (hcal_dx_1d_allcuts_600) {
    hcal_dx_1d_allcuts_600->SetXTitle("hcal_dx");
  }

  C->Draw("hcal_dx>>hcal_dx_1d_allcuts_700(700, -4, 4)", Form("corrected_weight*(%s)", allcuts_study_string.c_str()), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_dx_1d_allcuts_700= (TH1D*)gDirectory->Get("hcal_dx_1d_allcuts_700");
  if (hcal_dx_1d_allcuts_700) {
    hcal_dx_1d_allcuts_700->SetXTitle("hcal_dx");
  }

  C->Draw("hcal_dx>>hcal_dx_1d_allcuts_800(800, -4, 4)", Form("corrected_weight*(%s)", allcuts_study_string.c_str()), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_dx_1d_allcuts_800= (TH1D*)gDirectory->Get("hcal_dx_1d_allcuts_800");
  if (hcal_dx_1d_allcuts_800) {
    hcal_dx_1d_allcuts_800->SetXTitle("hcal_dx");
  }

  C->Draw("hcal_dx>>hcal_dx_1d_allcuts_900(900, -4, 4)", Form("corrected_weight*(%s)", allcuts_study_string.c_str()), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_dx_1d_allcuts_900= (TH1D*)gDirectory->Get("hcal_dx_1d_allcuts_900");
  if (hcal_dx_1d_allcuts_900) {
    hcal_dx_1d_allcuts_900->SetXTitle("hcal_dx");
  }

  C->Draw("hcal_dx>>hcal_dx_1d_allcuts_1000(1000, -4, 4)", Form("corrected_weight*(%s)", allcuts_study_string.c_str()), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_dx_1d_allcuts_1000= (TH1D*)gDirectory->Get("hcal_dx_1d_allcuts_1000");
  if (hcal_dx_1d_allcuts_1000) {
    hcal_dx_1d_allcuts_1000->SetXTitle("hcal_dx");
  }


  

 
  fout ->Write();
}// end main
