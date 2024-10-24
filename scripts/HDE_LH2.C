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
#include <string>
using namespace std;



#include "/work/halla//sbs/msatnik/GMn/classes/Utility.h"
#include "/work/halla//sbs/msatnik/GMn/classes/Utility.cpp"


// Positions (mc)
static const Double_t hcalposXi_mc = -2.655;    //m, distance from beam center to top of HCal w/75cm offset
static const Double_t hcalposXf_mc = 1.155;     //m, distance from beam center to bottom of HCal w/75cm offset
static const Double_t hcalposYi_mc = -0.92964;  //m, distance from beam center to opposite-beam side of HCal
static const Double_t hcalposYf_mc = 0.92964;   //m, distance from beam center to beam side of HCal

const double hcal_x_exp_min = -1.3;
const double hcal_x_exp_max = 0.5;
const double hcal_y_exp_min = -0.35;
const double hcal_y_exp_max = 0.4;

// functions
void adjustCanvas(TCanvas* canvas,
                  double leftMargin = 0.15, double rightMargin = 0.05, 
                  double bottomMargin = 0.15, double topMargin = 0.10);
void printParsedTitle(const std::string& title);



void HDE_LH2(TString configfileinput="sbs4_30p_cuts_LH2"){ // main
  //set draw params
  gStyle->SetNumberContours(255); 
  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(1111);
  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0110);
  gStyle->SetStatTextColor(kBlack);
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.06);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);

  cout<<"Loading Rootfiles. Hold tight!"<<endl;

  Utility utilityHandler; // class that gives us access to various functions to use between programs.

  TChain *C = new TChain("P");

  //string configfilename = Form("../config/sbs%d.cfg",kine);
  TString configfilename = "../config/" + configfileinput + ".cfg";
  // string configfilename = "/w/halla-scshelf2102/sbs/msatnik/GMn/config/sbs4_30p.cfg";
  cout<<"reading from config file: "<<configfilename<<endl;


  // set location and name for output file 
  TString outputfilename = "../output/" + configfileinput + "_HDE.root";
  //TString outputfilename = "../output/test.root";

  // Declare outfile
  // TFile *fout = new TFile( Form("../output/sbs%d.root",kine), "RECREATE" );
  //TFile *fout = new TFile("../output/sbs4_30p.root", "RECREATE" );
  TFile *fout = new TFile(outputfilename,"RECREATE");
  cout<<"writing to file: "<< outputfilename <<endl;


  double E_e=0;

  // Setting up cuts. These will be overwritten from the config files. 
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
  std::string GRINCH_CutString= "bb_grinch_clus_size>=3";
  std::string Optics_CutString = "abs(bb_tr_r_x-bb_tr_r_th*0.9)<0.3";
  std::string ProtonSpot_CutString ="pow((hcal_dx+0.71)/0.21,2)+pow(hcal_dy/0.3,2)<=1";
  std::string NeutronSpot_CutString="pow(hcal_dx/0.20,2)+pow(hcal_dy/0.3,2)<=1";
  double ProtonSpot_offset = 0; 
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
    //cout<< "Global Cut: "<<globalcut<<endl;
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
      if( skey == "GRINCH_CutString" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        GRINCH_CutString= sval.Data();
	cout << "Loading GRINCH_CutString: " << GRINCH_CutString<< endl;
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
      if( skey == "ProtonSpot_CutString" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	ProtonSpot_CutString= sval.Data();
	cout << "Loading ProtonSpot_Cutstring: " << ProtonSpot_CutString<< endl;
      }
      if( skey == "NeutronSpot_CutString" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	NeutronSpot_CutString= sval.Data();
	cout << "Loading NeutronSpot_Cutstring: " << NeutronSpot_CutString<< endl;
      }
      if( skey == "ProtonSpot_offset" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	ProtonSpot_offset = sval.Atof();
	cout << "Loading ProtonSpot Offset: " <<  ProtonSpot_offset<< endl;
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
  int passed_atime_cuts;

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
  C->SetBranchAddress("passed_atime_cuts",&passed_atime_cuts);

  C->SetBranchAddress("is_proton",&is_proton);
  C->SetBranchAddress("is_neutron",&is_neutron);
  C->SetBranchAddress("corrected_weight",&corrected_weight);


  



 /// 1d histo of hcal_x_exp with proton spot and HCal variables. 
  std::string proton_detected_hcal_x_exp_study_string= EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+ HCal_Energy_CutString + "&&"+  HCal_Shower_atime_CutString +"&&"+e_over_p_CutString +"&&"+ W2CutString+"&&("+ ProtonSpot_CutString +")";
  //TrackQualityCutString;
  //// Draw the 2D histogram
  C->Draw("hcal_x_exp>>proton_detected_hcal_x_exp_hist(300, -3, 3)",  proton_detected_hcal_x_exp_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *proton_detected_hcal_x_exp_hist= (TH1D*)gDirectory->Get("proton_detected_hcal_x_exp_hist");
  //TH1D *shifted_proton_detected_hcal_x_exp_hist = utilityHandler.ScaleAndShiftHistogram(proton_detected_hcal_x_exp_hist,1,ProtonSpot_offset);// hist, scale, shift
  TH1D *shifted_proton_detected_hcal_x_exp_hist =(TH1D*)proton_detected_hcal_x_exp_hist->Clone("shifted_proton_detected_hcal_x_exp_hist");
  if (shifted_proton_detected_hcal_x_exp_hist) {
    shifted_proton_detected_hcal_x_exp_hist->SetXTitle("detected proton hyp: hcal x expected");
  }
  

  std::string detected_title = proton_detected_hcal_x_exp_hist->GetTitle();
  // cout<<detected_title<<endl;
  
  /// 1d histo of hcal_x_exp  with no HCal variables. . 
  std::string proton_expected_hcal_x_exp_study_string =EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+e_over_p_CutString +"&&"+ W2CutString;
  //TrackQualityCutString+"&&"+ ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_x_exp>>proton_expected_hcal_x_exp_hist(300, -3, 3)",  proton_expected_hcal_x_exp_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *proton_expected_hcal_x_exp_hist= (TH1D*)gDirectory->Get("proton_expected_hcal_x_exp_hist");
  //TH1D *shifted_proton_expected_hcal_x_exp_hist = utilityHandler.ScaleAndShiftHistogram(proton_expected_hcal_x_exp_hist, 1, ProtonSpot_offset);
  TH1D *shifted_proton_expected_hcal_x_exp_hist = (TH1D*)proton_expected_hcal_x_exp_hist->Clone("shifted_proton_expected_hcal_x_exp_hist");
  if (shifted_proton_expected_hcal_x_exp_hist) {
    shifted_proton_expected_hcal_x_exp_hist->SetXTitle("proton hyp expected: hcal x expected");
  }
  std::string expected_title = shifted_proton_expected_hcal_x_exp_hist->GetTitle();
  //cout<<expected_title<<endl;


  // Get the detection eff by the ratio, bin-by-bin, by dividing hcal_x_exp histogram with HCal cuts (detected) by the hcal_x_exp histogram without hcal cuts (expected)
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH1D *HDE_hcal_x_exp_hist = (TH1D*)shifted_proton_detected_hcal_x_exp_hist->Clone("HDE_hcal_x_exp_hist"); 
  HDE_hcal_x_exp_hist->SetXTitle("proton hcal x expected");
  HDE_hcal_x_exp_hist->SetYTitle("Proton Detecton Eff");
  HDE_hcal_x_exp_hist ->Divide(shifted_proton_expected_hcal_x_exp_hist);
  TF1 *HDE_hcal_x_exp_fit = new TF1("HDE_hcal_x_exp_fit","[0]",hcal_x_exp_min,hcal_x_exp_max);
  HDE_hcal_x_exp_hist ->Fit(HDE_hcal_x_exp_fit,"Q R");
  //HDE_hcal_x_exp_hist ->Draw("E");

  shifted_proton_detected_hcal_x_exp_hist->Print("base");
  
  shifted_proton_expected_hcal_x_exp_hist->Print("base");
  


  

 /// 1d histo of hcal_y_exp with proton spot and HCal variables. 
  std::string proton_detected_hcal_y_exp_study_string= EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+ HCal_Energy_CutString + "&&"+  HCal_Shower_atime_CutString +"&&"+e_over_p_CutString +"&&"+ W2CutString+"&&("+ ProtonSpot_CutString +")";
  //TrackQualityCutString;
  //// Draw the 2D histogram
  C->Draw("hcal_y_exp>>proton_detected_hcal_y_exp_hist(300, -3, 3)",  proton_detected_hcal_y_exp_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *proton_detected_hcal_y_exp_hist= (TH1D*)gDirectory->Get("proton_detected_hcal_y_exp_hist");
  //TH1D *shifted_proton_detected_hcal_y_exp_hist = utilityHandler.ScaleAndShiftHistogram(proton_detected_hcal_y_exp_hist,1,ProtonSpot_offset);// hist, scale, shift
  TH1D *shifted_proton_detected_hcal_y_exp_hist =(TH1D*)proton_detected_hcal_y_exp_hist->Clone("shifted_proton_detected_hcal_y_exp_hist");
  if (shifted_proton_detected_hcal_y_exp_hist) {
    shifted_proton_detected_hcal_y_exp_hist->SetXTitle("detected proton hyp: hcal y expected");
  }
  

  std::string detected_title_y = proton_detected_hcal_y_exp_hist->GetTitle();
  // cout<<detected_title<<endl;
  
  /// 1d histo of hcal_y_exp  with no HCal variables. . 
  std::string proton_expected_hcal_y_exp_study_string =EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+e_over_p_CutString +"&&"+ W2CutString;
  //TrackQualityCutString+"&&"+ ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_y_exp>>proton_expected_hcal_y_exp_hist(300, -3, 3)",  proton_expected_hcal_y_exp_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *proton_expected_hcal_y_exp_hist= (TH1D*)gDirectory->Get("proton_expected_hcal_y_exp_hist");
  //TH1D *shifted_proton_expected_hcal_y_exp_hist = utilityHandler.ScaleAndShiftHistogram(proton_expected_hcal_y_exp_hist, 1, ProtonSpot_offset);
  TH1D *shifted_proton_expected_hcal_y_exp_hist = (TH1D*)proton_expected_hcal_y_exp_hist->Clone("shifted_proton_expected_hcal_y_exp_hist");
  if (shifted_proton_expected_hcal_y_exp_hist) {
    shifted_proton_expected_hcal_y_exp_hist->SetXTitle("proton hyp expected: hcal y expected");
  }
  std::string expected_title_y = shifted_proton_expected_hcal_y_exp_hist->GetTitle();
  //cout<<expected_title<<endl;




  // Get the detection eff by the ratio, bin-by-bin, by dividing hcal_y_exp histogram with HCal cuts (detected) by the hcal_y_exp histogram without hcal cuts (expected)
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH1D *HDE_hcal_y_exp_hist = (TH1D*)shifted_proton_detected_hcal_y_exp_hist->Clone("HDE_hcal_y_exp_hist"); 
  HDE_hcal_y_exp_hist->SetXTitle("proton hcal y expected");
  HDE_hcal_y_exp_hist->SetYTitle("Proton Detecton Eff");
  HDE_hcal_y_exp_hist ->Divide(shifted_proton_expected_hcal_y_exp_hist);
  TF1 *HDE_hcal_y_exp_fit = new TF1("HDE_hcal_y_exp_fit","[0]",hcal_y_exp_min,hcal_y_exp_max);
  HDE_hcal_y_exp_hist ->Fit(HDE_hcal_y_exp_fit,"Q R");
  //HDE_hcal_y_exp_hist ->Draw("E");

  shifted_proton_detected_hcal_y_exp_hist->Print("base");
  
  shifted_proton_expected_hcal_y_exp_hist->Print("base");
  

  

  TCanvas* HDE_hcal_y_exp_canvas = new TCanvas("HDE_hcal_y_exp_canvas", "HDE_hcal_y_exp_canvas", 1000, 600);
  HDE_hcal_y_exp_canvas->SetGrid();
  HDE_hcal_y_exp_canvas ->Divide(1,2);
  HDE_hcal_y_exp_canvas->cd(1);
  //gPad->SetLogy();
  proton_expected_hcal_y_exp_hist->SetTitle("hcal y expected");
  proton_expected_hcal_y_exp_hist->Draw("E");
  proton_detected_hcal_y_exp_hist->SetLineColor(kViolet);
  proton_detected_hcal_y_exp_hist->Draw("same E");
  TLegend *HDE_hcal_y_exp_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  HDE_hcal_y_exp_leg ->AddEntry(proton_expected_hcal_y_exp_hist,"No HCal Cuts","l");
  HDE_hcal_y_exp_leg ->AddEntry(proton_detected_hcal_y_exp_hist,"HCal Cuts and Proton Spot Cut","l");
  HDE_hcal_y_exp_leg->Draw();
  HDE_hcal_y_exp_canvas->cd(2);
  HDE_hcal_y_exp_hist->Draw("E");


  TCanvas* HDE_hcal_x_exp_canvas = new TCanvas("HDE_hcal_x_exp_canvas", "HDE_hcal_x_exp_canvas", 1000, 600);
  HDE_hcal_x_exp_canvas->SetGrid();
  HDE_hcal_x_exp_canvas ->Divide(1,2);
  HDE_hcal_x_exp_canvas->cd(1);
  //gPad->SetLogy();
  proton_expected_hcal_x_exp_hist->SetTitle("hcal x expected");
  proton_expected_hcal_x_exp_hist->Draw("E");
  proton_detected_hcal_x_exp_hist->SetLineColor(kViolet);
  proton_detected_hcal_x_exp_hist->Draw("same E");
  TLegend *HDE_hcal_x_exp_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  HDE_hcal_x_exp_leg ->AddEntry(proton_expected_hcal_x_exp_hist,"No HCal Cuts","l");
  HDE_hcal_x_exp_leg ->AddEntry(proton_detected_hcal_x_exp_hist,"HCal Cuts and Proton Spot cut","l");
  HDE_hcal_x_exp_leg->Draw();
  HDE_hcal_x_exp_canvas->cd(2);
  HDE_hcal_x_exp_hist->Draw("E");

 
 
  printParsedTitle(detected_title);
 
  fout ->Write();
}// end main


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



/// This expects the title to be in the form:
///  y_axis:x_axis {cut1&&cut2&&cut3....}
void printParsedTitle(const std::string& title) {
  
  cout<<title<<endl;
  
  // Find the position of the first '{' character
  size_t pos = title.find('{');
    
  // Extract the y_axis:x_axis part
  std::string axes = title.substr(0, pos);
  //cout<<"axis = "<<axes<<endl;
    
  // Extract the cuts part and remove '{' and '}'
  std::string cuts = title.substr(pos + 1, title.size() - pos - 2);
  //cout<<"cuts = "<<cuts<<endl;
   
  //cout<<"broken up cuts"<<endl;

  // Split the cuts into individual cut expressions
  std::vector<std::string> cutList;
  std::stringstream ss(cuts);
  std::string cut;
  while (std::getline(ss, cut, '&')) {
    // Remove leading and trailing whitespace
    cut.erase(0, cut.find_first_not_of(" \t"));
    cut.erase(cut.find_last_not_of(" \t") + 1);
        
    // Ensure the cut is not empty before processing
    if (!cut.empty() && cut.front() == '&') {
      cut.erase(cut.begin());
    }
        
    if (!cut.empty()) {
      cutList.push_back(cut);
      // cout<<cut<<endl;
    }
  }

    
  // Create a new canvas
  TCanvas* cuts_canvas = new TCanvas("cuts_canvas", "Parsed Histogram Title", 1000, 600);
    
  // Create a TLatex object to draw the text
  TLatex latex;
  latex.SetTextSize(0.03);  // Adjust text size
  latex.SetTextAlign(13);   // Align text to top left
    
  // Draw the axes part
  latex.DrawLatex(0.1, 0.9, axes.c_str());
    
  // Draw each cut expression on a new line
  double yPos = 0.8;  // Start position for the first cut
  for (const auto& cut : cutList) {
    latex.DrawLatex(0.1, yPos, cut.c_str());
    yPos -= 0.03;  // Move down for the next cut
  }
    
  // Update the canvas
  cuts_canvas->Update();
    
  // // Optionally, save the canvas as an image
  // // canvas->SaveAs(Form("%s/global_cuts_%s.pdf", output_dir.c_str(),outputname));
}
