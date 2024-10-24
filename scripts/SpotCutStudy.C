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



void SpotCutStudy(TString configfileinput="sbs4_30p_cuts"){ // main
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

  TChain *C = new TChain("P");

  //string configfilename = Form("../config/sbs%d.cfg",kine);
  TString configfilename = "../config/" + configfileinput + ".cfg";
  // string configfilename = "/w/halla-scshelf2102/sbs/msatnik/GMn/config/sbs4_30p.cfg";
  cout<<"reading from config file: "<<configfilename<<endl;


  // set location and name for output file 
  TString outputfilename = "../output/" + configfileinput + "_SpotStudy.root";
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



  




  // dy vs dx study
  std::string dy_study_string =  EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_dx:hcal_dy>>hcal_dx__hcal_dy(100, -2, 2, 200, -4, 4)", dy_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH2D *hcal_dx__hcal_dy= (TH2D*)gDirectory->Get("hcal_dx__hcal_dy");
  if (hcal_dx__hcal_dy) {
    hcal_dx__hcal_dy->SetXTitle("hcal_dy");
    hcal_dx__hcal_dy->SetYTitle("hcal_dx");
  }


  // seeing how the spot cuts look 
  // dy vs dx with proton spot cuts 
  std::string p_spot_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+ "&&" + ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_dx:hcal_dy>>hcal_dx__hcal_dy_protonspot(100, -2, 2, 200, -4, 4)", p_spot_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH2D *hcal_dx__hcal_dy_protonspot= (TH2D*)gDirectory->Get("hcal_dx__hcal_dy_protonspot");
  if (hcal_dx__hcal_dy_protonspot) {
    hcal_dx__hcal_dy_protonspot->SetXTitle("hcal_dy");
    hcal_dx__hcal_dy_protonspot->SetYTitle("hcal_dx");
  }

  // dy vs dx with neutron spot cuts 
  std::string n_spot_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+ "&&" + NeutronSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_dx:hcal_dy>>hcal_dx__hcal_dy_neutronspot(100, -2, 2, 200, -4, 4)", n_spot_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH2D *hcal_dx__hcal_dy_neutronspot= (TH2D*)gDirectory->Get("hcal_dx__hcal_dy_neutronspot");
  if (hcal_dx__hcal_dy_neutronspot) {
    hcal_dx__hcal_dy_neutronspot->SetXTitle("hcal_dy");
    hcal_dx__hcal_dy_neutronspot->SetYTitle("hcal_dx");
  }


  // Studying n/p over HCAL_Y 

  /// 1d histo of hcal_y 
  std::string hcal_y_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+ HCal_Energy_CutString + "&&"+  HCal_Shower_atime_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_y>>hcal_y_hist(200, -1, 1)",  hcal_y_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_y_hist= (TH1D*)gDirectory->Get("hcal_y_hist");
  if (hcal_y_hist) {
    hcal_y_hist->SetXTitle("hcal y");
  }

  /// 1d histo of hcal_y  with proton spot cut. 
  std::string proton_hcal_y_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+HCal_Energy_CutString  + "&&"+  HCal_Shower_atime_CutString +"&&" + ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_y>>hcal_y_hist_proton(200, -1, 1)",  proton_hcal_y_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *hcal_y_hist_proton= (TH1D*)gDirectory->Get("hcal_y_hist_proton");
  if (hcal_y_hist_proton) {
    hcal_y_hist_proton->SetXTitle("hcal y with proton spot cut");
  }

  /// 1d histo of hcal_y with neutron spot cut. 
  std::string neutron_hcal_y_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+HCal_Energy_CutString +"&&"+HCal_Shower_atime_CutString +"&&"+NeutronSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_y>>hcal_y_hist_neutron(200, -1, 1)",  neutron_hcal_y_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *hcal_y_hist_neutron= (TH1D*)gDirectory->Get("hcal_y_hist_neutron");
  if (hcal_y_hist_neutron) {
    hcal_y_hist_neutron->SetXTitle("hcal y with neutron spot cut");
  }


  // Get the neutron/proton ratio, bin-by-bin, by dividing the hcal_y histo with the neutron spot cut by the hcal_y histo with the proton spot cut.
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH1D *np_hcal_y_hist = (TH1D*)hcal_y_hist_neutron->Clone("np_hcal_y_hist"); 
  np_hcal_y_hist ->Divide(hcal_y_hist_proton);
  np_hcal_y_hist->SetXTitle("hcal y");
  np_hcal_y_hist->SetYTitle("n/p");
  np_hcal_y_hist ->Draw("E");


  // Studying n/p over hcal_y_exp 

  /// 1d histo of hcal_y_exp 
  std::string hcal_y_exp_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+ HCal_Energy_CutString + "&&"+  HCal_Shower_atime_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_y_exp>>hcal_y_exp_hist(200, -1, 1)",  hcal_y_exp_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_y_exp_hist= (TH1D*)gDirectory->Get("hcal_y_exp_hist");
  if (hcal_y_exp_hist) {
    hcal_y_exp_hist->SetXTitle("hcal y expected");
  }

  /// 1d histo of hcal_y_exp  with proton spot cut. 
  std::string proton_hcal_y_exp_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+HCal_Energy_CutString  + "&&"+  HCal_Shower_atime_CutString +"&&" + ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_y_exp>>hcal_y_exp_hist_proton(200, -1, 1)",  proton_hcal_y_exp_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *hcal_y_exp_hist_proton= (TH1D*)gDirectory->Get("hcal_y_exp_hist_proton");
  if (hcal_y_exp_hist_proton) {
    hcal_y_exp_hist_proton->SetXTitle("hcal y expected with proton spot cut");
  }

  /// 1d histo of hcal_y_exp with neutron spot cut. 
  std::string neutron_hcal_y_exp_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+HCal_Energy_CutString +"&&"+HCal_Shower_atime_CutString +"&&"+NeutronSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_y_exp>>hcal_y_exp_hist_neutron(200, -1, 1)",  neutron_hcal_y_exp_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *hcal_y_exp_hist_neutron= (TH1D*)gDirectory->Get("hcal_y_exp_hist_neutron");
  if (hcal_y_exp_hist_neutron) {
    hcal_y_exp_hist_neutron->SetXTitle("hcal y expected with neutron spot cut");
  }


  // Get the neutron/proton ratio, bin-by-bin, by dividing the hcal_y_exp histo with the neutron spot cut by the hcal_y_exp histo with the proton spot cut.
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH1D *np_hcal_y_exp_hist = (TH1D*)hcal_y_exp_hist_neutron->Clone("np_hcal_y_exp_hist"); 
  np_hcal_y_exp_hist->SetXTitle("hcal y expected");
  np_hcal_y_exp_hist->SetYTitle("n/p");
  np_hcal_y_exp_hist ->Divide(hcal_y_exp_hist_proton);
  TF1 *np_hcal_y_exp_fit = new TF1("np_hcal_y_exp_fit","[0]",hcal_y_exp_min,hcal_y_exp_max);
  np_hcal_y_exp_hist ->Fit(np_hcal_y_exp_fit,"Q R");
  np_hcal_y_exp_hist ->Draw("E");



  // Studying n/p over HCAL_X 

  /// 1d histo of hcal_x 
  std::string hcal_x_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+ HCal_Energy_CutString + "&&"+  HCal_Shower_atime_CutString +"&&" + W2CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_x>>hcal_x_hist(450, -3, 1.5)",  hcal_x_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_x_hist= (TH1D*)gDirectory->Get("hcal_x_hist");
  if (hcal_x_hist) {
    hcal_x_hist->SetXTitle("hcal x");
  }

  /// 1d histo of hcal_x  with proton spot cut. 
  std::string proton_hcal_x_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+HCal_Energy_CutString  + "&&"+  HCal_Shower_atime_CutString +"&&" +W2CutString+"&&"+ ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_x>>hcal_x_hist_proton(450, -3, 1.5)",  proton_hcal_x_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *hcal_x_hist_proton= (TH1D*)gDirectory->Get("hcal_x_hist_proton");
  if (hcal_x_hist_proton) {
    hcal_x_hist_proton->SetXTitle("hcal x with proton spot cut");
  }

  /// 1d histo of hcal_x with neutron spot cut. 
  std::string neutron_hcal_x_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+HCal_Energy_CutString +"&&"+HCal_Shower_atime_CutString +"&&" + W2CutString+"&&"+NeutronSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_x>>hcal_x_hist_neutron(450, -3, 1.5)",  neutron_hcal_x_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *hcal_x_hist_neutron= (TH1D*)gDirectory->Get("hcal_x_hist_neutron");
  if (hcal_x_hist_neutron) {
    hcal_x_hist_neutron->SetXTitle("hcal x with neutron spot cut");
  }


  // Get the neutron/proton ratio, bin-by-bin, by dividing the hcal_x histo with the neutron spot cut by the hcal_x histo with the proton spot cut.
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH1D *np_hcal_x_hist = (TH1D*)hcal_x_hist_neutron->Clone("np_hcal_x_hist"); 
  np_hcal_x_hist->SetXTitle("hcal x");
  np_hcal_x_hist->SetYTitle("n/p");
  np_hcal_x_hist ->Divide(hcal_x_hist_proton);
  np_hcal_x_hist ->Draw("E");




  // Studying n/p over HCAL_X Expected

  /// 1d histo of hcal_x_exp 
  std::string hcal_x_exp_study_string= EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+ HCal_Energy_CutString + "&&"+  HCal_Shower_atime_CutString +"&&"+e_over_p_CutString;
  //TrackQualityCutString;
  //// Draw the 2D histogram
  C->Draw("hcal_x_exp>>hcal_x_exp_hist(600, -3, 3)",  hcal_x_exp_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_x_exp_hist= (TH1D*)gDirectory->Get("hcal_x_exp_hist");
  if (hcal_x_exp_hist) {
    hcal_x_exp_hist->SetXTitle("hcal x expected");
  }

  /// 1d histo of hcal_x_exp  with proton spot cut. 
  std::string proton_hcal_x_exp_study_string =EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+HCal_Energy_CutString  + "&&"+  HCal_Shower_atime_CutString +"&&"+e_over_p_CutString+"&&"+ ProtonSpot_CutString;
  //TrackQualityCutString+"&&"+ ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_x_exp>>hcal_x_exp_hist_proton(600, -3, 3)",  proton_hcal_x_exp_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *hcal_x_exp_hist_proton= (TH1D*)gDirectory->Get("hcal_x_exp_hist_proton");
  if (hcal_x_exp_hist_proton) {
    hcal_x_exp_hist_proton->SetXTitle("hcal x expected with proton spot cut");
  }

  /// 1d histo of hcal_x_exp with neutron spot cut. 
  std::string neutron_hcal_x_exp_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+HCal_Energy_CutString +"&&"+HCal_Shower_atime_CutString +"&&"+e_over_p_CutString+"&&"+NeutronSpot_CutString;
  //TargetVertexCutString+"&&"+NeutronSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_x_exp>>hcal_x_exp_hist_neutron(600, -3, 3)",  neutron_hcal_x_exp_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *hcal_x_exp_hist_neutron= (TH1D*)gDirectory->Get("hcal_x_exp_hist_neutron");
  if (hcal_x_exp_hist_neutron) {
    hcal_x_exp_hist_neutron->SetXTitle("hcal x expected with neutron spot cut");
  }


  // Get the neutron/proton ratio, bin-by-bin, by dividing the hcal_x_exp histo with the neutron spot cut by the hcal_x_exp histo with the proton spot cut.
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH1D *np_hcal_x_exp_hist = (TH1D*)hcal_x_exp_hist_neutron->Clone("np_hcal_x_exp_hist"); 
  np_hcal_x_exp_hist->SetXTitle("hcal x expected");
  np_hcal_x_exp_hist->SetYTitle("n/p");
  np_hcal_x_exp_hist ->Divide(hcal_x_exp_hist_proton);
  TF1 *np_hcal_x_exp_fit = new TF1("np_hcal_x_exp_fit","[0]",hcal_x_exp_min,hcal_x_exp_max);
  np_hcal_x_exp_hist ->Fit(np_hcal_x_exp_fit,"Q R");
  np_hcal_x_exp_hist ->Draw("E");




  // Studying n/p over HCAL_X Exp vs HCal Y Exp 2d histo

  /// 2d histo of hcal_x_exp vs hcal_y_exp 
  std::string hcal_x_exp_hcal_y_exp_study_string= EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+ HCal_Energy_CutString + "&&"+  HCal_Shower_atime_CutString +"&&"+e_over_p_CutString;
  //TrackQualityCutString;
  //// Draw the 2D histogram
  C->Draw("hcal_x_exp:hcal_y_exp>>hcal_x_exp_hcal_y_exp_hist(20,-2,2,30, -3, 3)",  hcal_x_exp_hcal_y_exp_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_x_exp_hcal_y_exp_hist= (TH1D*)gDirectory->Get("hcal_x_exp_hcal_y_exp_hist");
  if (hcal_x_exp_hcal_y_exp_hist) {
    hcal_x_exp_hcal_y_exp_hist->SetXTitle("hcal y expected");
    hcal_x_exp_hcal_y_exp_hist->SetXTitle("hcal x expected");
  }

  /// 2d histo of hcal_x_exp vs hcal y exp  with proton spot cut. 
  std::string proton_hcal_x_exp_hcal_y_exp_study_string =EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+HCal_Energy_CutString  + "&&"+  HCal_Shower_atime_CutString +"&&"+e_over_p_CutString+"&&"+ ProtonSpot_CutString;
  //TrackQualityCutString+"&&"+ ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_x_exp:hcal_y_exp>>hcal_x_exp_hcal_y_exp_hist_proton(20,-2,2,30, -3, 3)",  proton_hcal_x_exp_hcal_y_exp_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_x_exp_hcal_y_exp_hist_proton= (TH1D*)gDirectory->Get("hcal_x_exp_hcal_y_exp_hist_proton");
  if (hcal_x_exp_hcal_y_exp_hist_proton) {
    hcal_x_exp_hcal_y_exp_hist_proton->SetXTitle("hcal y expected with proton spot cut");
    hcal_x_exp_hcal_y_exp_hist_proton->SetYTitle("hcal x expected with proton spot cut");
  }

  /// 1d histo of hcal_x_exp with neutron spot cut. 
  std::string neutron_hcal_x_exp_hcal_y_exp_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+HCal_Energy_CutString +"&&"+HCal_Shower_atime_CutString +"&&"+e_over_p_CutString+"&&"+NeutronSpot_CutString;
  //TargetVertexCutString+"&&"+NeutronSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_x_exp:hcal_y_exp>>hcal_x_exp_hcal_y_exp_hist_neutron(20,-2,2,30, -3, 3)",  neutron_hcal_x_exp_hcal_y_exp_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_x_exp_hcal_y_exp_hist_neutron= (TH1D*)gDirectory->Get("hcal_x_exp_hcal_y_exp_hist_neutron");
  if (hcal_x_exp_hcal_y_exp_hist_neutron) {
    hcal_x_exp_hcal_y_exp_hist_neutron->SetXTitle("hcal y expected with neutron spot cut");
    hcal_x_exp_hcal_y_exp_hist_neutron->SetYTitle("hcal x expected with neutron spot cut");
  }


  // Get the neutron/proton ratio, bin-by-bin, by dividing the hcal_x_exp histo with the neutron spot cut by the hcal_x_exp histo with the proton spot cut.
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH2D *np_hcal_x_exp_hcal_y_exp_hist = (TH2D*)hcal_x_exp_hcal_y_exp_hist_neutron->Clone("np_hcal_x_exp_hcal_y_exp_hist"); 
  np_hcal_x_exp_hcal_y_exp_hist ->SetXTitle("hcal y expected");
  np_hcal_x_exp_hcal_y_exp_hist  ->SetYTitle("hcal x expected");
  np_hcal_x_exp_hcal_y_exp_hist  ->Divide(hcal_x_exp_hcal_y_exp_hist_proton);

  



  // Studying the Nucleon Detection Eff  over HCAL_X Expected

  /// 1d histo of hcal_x_exp with proton and neutron spots and HCal variables. 
  std::string detected_hcal_x_exp_study_string= EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+ HCal_Energy_CutString + "&&"+  HCal_Shower_atime_CutString +"&&"+e_over_p_CutString +"&&"+ W2CutString+"&&("+ ProtonSpot_CutString +"||"+ NeutronSpot_CutString+")";
  //TrackQualityCutString;
  //// Draw the 2D histogram
  C->Draw("hcal_x_exp>>detected_hcal_x_exp_hist(600, -3, 3)",  detected_hcal_x_exp_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *detected_hcal_x_exp_hist= (TH1D*)gDirectory->Get("detected_hcal_x_exp_hist");
  if (detected_hcal_x_exp_hist) {
    detected_hcal_x_exp_hist->SetXTitle("detected: hcal x expected");
  }
  
  std::string detected_title = detected_hcal_x_exp_hist->GetTitle();
  // cout<<detected_title<<endl;
  
  /// 1d histo of hcal_x_exp  with no HCal variables. . 
  std::string expected_hcal_x_exp_study_string =EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+e_over_p_CutString +"&&"+ W2CutString;
  //TrackQualityCutString+"&&"+ ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_x_exp>>expected_hcal_x_exp_hist(600, -3, 3)",  expected_hcal_x_exp_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *expected_hcal_x_exp_hist= (TH1D*)gDirectory->Get("expected_hcal_x_exp_hist");
  if (expected_hcal_x_exp_hist) {
    expected_hcal_x_exp_hist->SetXTitle("expected: hcal x expected");
  }
  std::string expected_title = expected_hcal_x_exp_hist->GetTitle();
  //cout<<expected_title<<endl;


  // Get the detection eff by the ratio, bin-by-bin, by dividing hcal_x_exp histogram with HCal cuts (detected) by the hcal_x_exp histogram without hcal cuts (expected)
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH1D *HDE_hcal_x_exp_hist = (TH1D*)detected_hcal_x_exp_hist->Clone("HDE_hcal_x_exp_hist"); 
  HDE_hcal_x_exp_hist->SetXTitle("hcal x expected");
  HDE_hcal_x_exp_hist->SetYTitle("Nucleon Detecton Eff");
  HDE_hcal_x_exp_hist ->Divide(expected_hcal_x_exp_hist);
  TF1 *HDE_hcal_x_exp_fit = new TF1("HDE_hcal_x_exp_fit","[0]",hcal_x_exp_min,hcal_x_exp_max);
  HDE_hcal_x_exp_hist ->Fit(HDE_hcal_x_exp_fit,"Q R");
  HDE_hcal_x_exp_hist ->Draw("E");



  // Studying the Nucleon Detection Eff  over HCAL_Y Expected

  /// 1d histo of hcal_y_exp with proton and neutron spots and HCal variables. 
  std::string detected_hcal_y_exp_study_string= EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+ HCal_Energy_CutString + "&&"+  HCal_Shower_atime_CutString +"&&"+e_over_p_CutString +"&&"+ W2CutString+"&&("+ ProtonSpot_CutString +"||"+ NeutronSpot_CutString+")";
  //TrackQualityCutString;
  //// Draw the 2D histogram
  C->Draw("hcal_y_exp>>detected_hcal_y_exp_hist(400, -2, 2)",  detected_hcal_y_exp_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *detected_hcal_y_exp_hist= (TH1D*)gDirectory->Get("detected_hcal_y_exp_hist");
  if (detected_hcal_y_exp_hist) {
    detected_hcal_y_exp_hist->SetXTitle("detected: hcal y expected");
  }
  
  
  /// 1d histo of hcal_y_exp  with no HCal variables. . 
  std::string expected_hcal_y_exp_study_string =EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+e_over_p_CutString +"&&"+ W2CutString;
  //TrackQualityCutString+"&&"+ ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_y_exp>>expected_hcal_y_exp_hist(400, -2, 2)",  expected_hcal_y_exp_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *expected_hcal_y_exp_hist= (TH1D*)gDirectory->Get("expected_hcal_y_exp_hist");
  if (expected_hcal_y_exp_hist) {
    expected_hcal_y_exp_hist->SetXTitle("expected: hcal y expected");
  }
 

  // Get the detection eff by the ratio, bin-by-bin, by dividing hcal_y_exp histogram with HCal cuts (detected) by the hcal_y_exp histogram without hcal cuts (expected)
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH1D *HDE_hcal_y_exp_hist = (TH1D*)detected_hcal_y_exp_hist->Clone("HDE_hcal_y_exp_hist"); 
  HDE_hcal_y_exp_hist->SetXTitle("hcal y expected");
  HDE_hcal_y_exp_hist->SetYTitle("Nucleon Detecton Eff");
  HDE_hcal_y_exp_hist ->Divide(expected_hcal_y_exp_hist);
  TF1 *HDE_hcal_y_exp_fit = new TF1("HDE_hcal_y_exp_fit","[0]",hcal_y_exp_min,hcal_y_exp_max);
  HDE_hcal_y_exp_hist ->Fit(HDE_hcal_y_exp_fit,"Q R");
  HDE_hcal_y_exp_hist ->Draw("E");


  // Studying the Nucleon Detection Eff  over 2d HCAL X Expected vs HCAL Y expected 

  /// 2d histo of hcal_x_exp vs hcal_y_exp with proton and neutron spots and HCal variables. 
  std::string detected_hcal_x_exp_hcal_y_exp_study_string= EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+ HCal_Energy_CutString + "&&"+  HCal_Shower_atime_CutString +"&&"+e_over_p_CutString +"&&"+ W2CutString+"&&("+ ProtonSpot_CutString +"||"+ NeutronSpot_CutString+")";
  //// Draw the histogram
  C->Draw("hcal_x_exp:hcal_y_exp>>detected_hcal_x_exp_hcal_y_exp_hist(75,-3,3,50, -2, 2)",  detected_hcal_x_exp_hcal_y_exp_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH2D *detected_hcal_x_exp_hcal_y_exp_hist= (TH2D*)gDirectory->Get("detected_hcal_x_exp_hcal_y_exp_hist");
  if (detected_hcal_x_exp_hcal_y_exp_hist) {
    detected_hcal_x_exp_hcal_y_exp_hist->SetXTitle("detected: hcal y exp");
    detected_hcal_x_exp_hcal_y_exp_hist->SetYTitle("detected: hcal x exp");
  }

  /// 2d histo of hcal_x_exp vs hcal_y_exp with NO HCAL or Spot cuts 
  std::string expected_hcal_x_exp_hcal_y_exp_study_string= EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+e_over_p_CutString +"&&"+ W2CutString;
  //// Draw the histogram
  C->Draw("hcal_x_exp:hcal_y_exp>>expected_hcal_x_exp_hcal_y_exp_hist(75,-3,3,50, -2, 2)",  expected_hcal_x_exp_hcal_y_exp_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH2D *expected_hcal_x_exp_hcal_y_exp_hist= (TH2D*)gDirectory->Get("expected_hcal_x_exp_hcal_y_exp_hist");
  if (expected_hcal_x_exp_hcal_y_exp_hist) {
    expected_hcal_x_exp_hcal_y_exp_hist->SetXTitle("expected: hcal y exp");
    expected_hcal_x_exp_hcal_y_exp_hist->SetYTitle("expected: hcal x exp");
  }
  
  

  // Get the detection eff by the ratio, bin-by-bin, by dividing 2d histogram with HCal cuts (detected) by the 2d histogram without hcal cuts (expected)
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH2D *HDE_hcal_x_exp_hcal_y_exp_hist = (TH2D*)detected_hcal_x_exp_hcal_y_exp_hist->Clone("HDE_hcal_x_exp_hcal_y_exp_hist"); 
  HDE_hcal_x_exp_hcal_y_exp_hist->SetXTitle("hcal y expected");
  HDE_hcal_x_exp_hcal_y_exp_hist->SetYTitle("hcal x expected");
  HDE_hcal_x_exp_hcal_y_exp_hist ->Divide(expected_hcal_x_exp_hcal_y_exp_hist);

  
  
  // Studying n/p over W2 

  /// 1d histo of W2 
  std::string W2_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+ HCal_Energy_CutString + "&&"+  HCal_Shower_atime_CutString +"&&"+FidXCutString;
  //// Draw the 2D histogram
  C->Draw("W2>>W2_hist(300, 0, 3)",  W2_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *W2_hist= (TH1D*)gDirectory->Get("W2_hist");
  if (W2_hist) {
    W2_hist->SetXTitle("W^{2}");
  }

  /// 1d histo of W2 with proton spot cut. 
  std::string proton_W2_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+HCal_Energy_CutString  + "&&"+  HCal_Shower_atime_CutString + "&&"+FidXCutString+"&&" + ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("W2>>W2_hist_proton(300, 0, 3)",  proton_W2_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *W2_hist_proton= (TH1D*)gDirectory->Get("W2_hist_proton");
  if (W2_hist_proton) {
    W2_hist_proton->SetXTitle("W^{2} with proton spot cut");
  }

  /// 1d histo of W2 with neutron spot cut. 
  std::string neutron_W2_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+HCal_Energy_CutString +"&&"+HCal_Shower_atime_CutString + "&&"+FidXCutString+"&&" + NeutronSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("W2>>W2_hist_neutron(300, 0, 3)",  neutron_W2_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *W2_hist_neutron= (TH1D*)gDirectory->Get("W2_hist_neutron");
  if (W2_hist_neutron) {
    W2_hist_neutron->SetXTitle("W^{2} with neutron spot cut");
  }

  // Get the neutron/proton ratio, bin-by-bin, by dividing the W2 histo with the neutron spot cut by the W2 histo with the proton spot cut.
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH1D *np_W2_hist = (TH1D*)W2_hist_neutron->Clone("np_W2_hist"); 
  np_W2_hist->SetXTitle("W^{2}");
  np_W2_hist->SetYTitle("n/p");
  np_W2_hist ->Divide(W2_hist_proton);
  TF1 *np_W2_fit = new TF1("np_W2_fit","[0]",0.66,1.10);
  np_W2_hist ->Fit(np_W2_fit,"Q R");
  np_W2_hist ->Draw("E");



   
  

  // Studying n/p for e_over_p

  /// 1d histo of e_over_p
  std::string e_over_p_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+ HCal_Energy_CutString +"&&"+FidXCutString +"&&" + W2CutString+ "&&"+  HCal_Shower_atime_CutString ;
  //// Draw the 2D histogram
  C->Draw("e_over_p>>e_over_p_hist(100, 0.5, 1.5)",  e_over_p_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *e_over_p_hist= (TH1D*)gDirectory->Get("e_over_p_hist");
  if (e_over_p_hist) {
    e_over_p_hist->SetXTitle("e/p");
  }

  /// 1d histo of e_over_p with proton spot cut. 
  std::string proton_e_over_p_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+HCal_Energy_CutString  + "&&"+FidXCutString+"&&"+W2CutString+"&&"+  HCal_Shower_atime_CutString +"&&" + ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("e_over_p>>e_over_p_hist_proton(100, 0.5, 1.5)",  proton_e_over_p_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *e_over_p_hist_proton= (TH1D*)gDirectory->Get("e_over_p_hist_proton");
  if (e_over_p_hist_proton) {
    e_over_p_hist_proton->SetXTitle("e/p");
  }

  /// 1d histo of e_over_p  with neutron spot cut. 
  std::string neutron_e_over_p_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+HCal_Energy_CutString + "&&"+FidXCutString+"&&" +W2CutString+"&&"+  HCal_Shower_atime_CutString +"&&" + NeutronSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("e_over_p>>e_over_p_hist_neutron(100, 0.5, 1.5)",  neutron_e_over_p_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *e_over_p_hist_neutron= (TH1D*)gDirectory->Get("e_over_p_hist_neutron");
  if (e_over_p_hist_neutron) {
    e_over_p_hist_neutron->SetXTitle("e/p");
  }


  // Get the neutron/proton ratio, bin-by-bin, by dividing the e_over_p histo with the neutron spot cut by the e_over_p histo with the proton spot cut.
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH1D *np_e_over_p_hist = (TH1D*)e_over_p_hist_neutron->Clone("np_e_over_p_hist"); 
  np_e_over_p_hist->SetXTitle("e/p");
  np_e_over_p_hist->SetYTitle("n/p");
  np_e_over_p_hist ->Divide(e_over_p_hist_proton);
  TF1 *np_e_over_p_fit = new TF1("np_e_over_p_fit","[0]",0.78,1.18);
  np_e_over_p_hist ->Fit(np_e_over_p_fit,"Q R");
  np_e_over_p_hist ->Draw("E");


  // Studying n/p for hcal_dy

  /// 1d histo of hcal_dy
  std::string hcal_dy_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+FidXCutString + "&&"+  HCal_Energy_CutString +"&&"+HCal_Shower_atime_CutString ;
  //// Draw the 2D histogram
  C->Draw("hcal_dy>>hcal_dy_hist(100, -2, 2)",  hcal_dy_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_dy_hist= (TH1D*)gDirectory->Get("hcal_dy_hist");
  if (hcal_dy_hist) {
    hcal_dy_hist->SetXTitle("hcal_dy");
  }

  /// 1d histo of hcal_dy with proton spot cut. 
  std::string proton_hcal_dy_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString  + "&&"+FidXCutString+"&&"+HCal_Energy_CutString +"&&"+HCal_Shower_atime_CutString +"&&" + ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_dy>>hcal_dy_hist_proton(100, -2, 2)",  proton_hcal_dy_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *hcal_dy_hist_proton= (TH1D*)gDirectory->Get("hcal_dy_hist_proton");
  if (hcal_dy_hist_proton) {
    hcal_dy_hist_proton->SetXTitle("hcal_dy");
  }

  /// 1d histo of hcal_dy  with neutron spot cut. 
  std::string neutron_hcal_dy_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+ "&&"+FidXCutString+"&&" +  HCal_Energy_CutString +"&&"+HCal_Shower_atime_CutString +"&&" + NeutronSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_dy>>hcal_dy_hist_neutron(100, -2, 2)",  neutron_hcal_dy_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *hcal_dy_hist_neutron= (TH1D*)gDirectory->Get("hcal_dy_hist_neutron");
  if (hcal_dy_hist_neutron) {
    hcal_dy_hist_neutron->SetXTitle("hcal_dy");
  }

  // Get the neutron/proton ratio, bin-by-bin, by dividing the hcal_dy histo with the neutron spot cut by the hcal_dy histo with the proton spot cut.
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH1D *np_hcal_dy_hist = (TH1D*)hcal_dy_hist_neutron->Clone("np_hcal_dy_hist"); 
  np_hcal_dy_hist->SetXTitle("hcal_dy");
  np_hcal_dy_hist->SetYTitle("n/p");
  np_hcal_dy_hist ->Divide(hcal_dy_hist_proton);
  np_hcal_dy_hist ->Draw("E");


  // Studying n/p for vz

  /// 1d histo of vz
  std::string vz_study_string = EnergyCutString + "&&"+ TrackQualityCutString +"&&" + Optics_CutString+"&&"+FidXCutString + "&&"+  HCal_Energy_CutString +"&&"+HCal_Shower_atime_CutString ;
  //// Draw the 2D histogram
  C->Draw("bb_tr_vz>>vz_hist(200, -0.1, 0.1)",  vz_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *vz_hist= (TH1D*)gDirectory->Get("vz_hist");
  if (vz_hist) {
    vz_hist->SetXTitle("vz");
  }

  /// 1d histo of vz with proton spot cut. 
  std::string proton_vz_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+ Optics_CutString  + "&&"+FidXCutString+"&&"+HCal_Energy_CutString +"&&"+HCal_Shower_atime_CutString +"&&" + ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("bb_tr_vz>>vz_hist_proton(200, -0.1, 0.1)",  proton_vz_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *vz_hist_proton= (TH1D*)gDirectory->Get("vz_hist_proton");
  if (vz_hist_proton) {
    vz_hist_proton->SetXTitle("vz");
  }

  /// 1d histo of vz  with neutron spot cut. 
  std::string neutron_vz_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+ Optics_CutString+ "&&"+FidXCutString+"&&" +  HCal_Energy_CutString +"&&"+HCal_Shower_atime_CutString +"&&" + NeutronSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("bb_tr_vz>>vz_hist_neutron(200, -0.1, 0.1)",  neutron_vz_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *vz_hist_neutron= (TH1D*)gDirectory->Get("vz_hist_neutron");
  if (vz_hist_neutron) {
    vz_hist_neutron->SetXTitle("vz");
  }

  // Get the neutron/proton ratio, bin-by-bin, by dividing the vz histo with the neutron spot cut by the vz histo with the proton spot cut.
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH1D *np_vz_hist = (TH1D*)vz_hist_neutron->Clone("np_vz_hist"); 
  np_vz_hist->SetXTitle("vz");
  np_vz_hist->SetYTitle("n/p");
  np_vz_hist ->Divide(vz_hist_proton);
  np_vz_hist ->Draw("E");


  // Studying n/p for hcal_e

  /// 1d histo of hcal_e
  std::string hcal_e_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+FidXCutString +"&&" + W2CutString+ "&&"+  HCal_Shower_atime_CutString ;
  //// Draw the 2D histogram
  C->Draw("hcal_e>>hcal_e_hist(200, 0, 0.8)",  hcal_e_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *hcal_e_hist= (TH1D*)gDirectory->Get("hcal_e_hist");
  if (hcal_e_hist) {
    hcal_e_hist->SetXTitle("HCal Energy");
  }

  /// 1d histo of hcal_e with proton spot cut. 
  std::string proton_hcal_e_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString  + "&&"+FidXCutString+"&&"+W2CutString+"&&"+  HCal_Shower_atime_CutString +"&&" + ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_e>>hcal_e_hist_proton(200, 0, 0.8)",  proton_hcal_e_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *hcal_e_hist_proton= (TH1D*)gDirectory->Get("hcal_e_hist_proton");
  if (hcal_e_hist_proton) {
    hcal_e_hist_proton->SetXTitle("HCal Energy");
  }

  /// 1d histo of hcal_e  with neutron spot cut. 
  std::string neutron_hcal_e_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+ "&&"+FidXCutString+"&&" +W2CutString+"&&"+  HCal_Shower_atime_CutString +"&&" + NeutronSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_e>>hcal_e_hist_neutron(200, 0, 0.8)",  neutron_hcal_e_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *hcal_e_hist_neutron= (TH1D*)gDirectory->Get("hcal_e_hist_neutron");
  if (hcal_e_hist_neutron) {
    hcal_e_hist_neutron->SetXTitle("HCal Energy");
  }

  // Get the neutron/proton ratio, bin-by-bin, by dividing the hcal_e histo with the neutron spot cut by the hcal_e histo with the proton spot cut.
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH1D *np_hcal_e_hist = (TH1D*)hcal_e_hist_neutron->Clone("np_hcal_e_hist"); 
  np_hcal_e_hist->SetXTitle("HCal Energy");
  np_hcal_e_hist->SetYTitle("n/p");
  np_hcal_e_hist ->Divide(hcal_e_hist_proton);
  double max_hcal_e = np_hcal_e_hist->GetXaxis()->GetXmax();
  TF1 *np_hcal_e_fit = new TF1("np_hcal_e_fit","[0]",0.025,max_hcal_e);
  np_hcal_e_hist ->Fit(np_hcal_e_fit,"Q R");
  np_hcal_e_hist ->Draw("E");




  // Studying n/p for ps_e

  /// 1d histo of ps_e
  std::string ps_e_study_string =  TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+FidXCutString +"&&" + W2CutString+"&&"+HCal_Energy_CutString + "&&"+  HCal_Shower_atime_CutString ;
  //// Draw the 2D histogram
  C->Draw("bb_ps_e>>ps_e_hist(200, 0, 2)",  ps_e_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *ps_e_hist= (TH1D*)gDirectory->Get("ps_e_hist");
  if (ps_e_hist) {
    ps_e_hist->SetXTitle("Preshower Energy");
  }

  /// 1d histo of ps_e with proton spot cut. 
  std::string proton_ps_e_study_string = TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString  + "&&"+FidXCutString+"&&"+W2CutString+"&&"+HCal_Energy_CutString+"&&"+  HCal_Shower_atime_CutString +"&&" + ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("bb_ps_e>>ps_e_hist_proton(200, 0, 2)",  proton_ps_e_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *ps_e_hist_proton= (TH1D*)gDirectory->Get("ps_e_hist_proton");
  if (ps_e_hist_proton) {
    ps_e_hist_proton->SetXTitle("Preshower Energy");
  }

  /// 1d histo of ps_e  with neutron spot cut. 
  std::string neutron_ps_e_study_string =  TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+ "&&"+FidXCutString+"&&" +W2CutString+"&&"+HCal_Energy_CutString+"&&"+  HCal_Shower_atime_CutString +"&&" + NeutronSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("bb_ps_e>>ps_e_hist_neutron(200, 0, 2)",  neutron_ps_e_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *ps_e_hist_neutron= (TH1D*)gDirectory->Get("ps_e_hist_neutron");
  if (ps_e_hist_neutron) {
    ps_e_hist_neutron->SetXTitle("Preshower Energy");
  }

  // Get the neutron/proton ratio, bin-by-bin, by dividing the ps_e histo with the neutron spot cut by the ps_e histo with the proton spot cut.
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH1D *np_ps_e_hist = (TH1D*)ps_e_hist_neutron->Clone("np_ps_e_hist"); 
  np_ps_e_hist ->Divide(ps_e_hist_proton);
  double max_ps_e = np_ps_e_hist->GetXaxis()->GetXmax();
  TF1 *ps_e_fit = new TF1("ps_e_fit","[0]",0.2,max_ps_e);
  np_ps_e_hist ->Fit(ps_e_fit,"Q R");
  np_ps_e_hist->SetXTitle("Preshower Energy");
  np_ps_e_hist->SetYTitle("n/p");
  np_ps_e_hist ->Draw("E");



  // Studying n/p for ps_e + sh_e

  /// 1d histo of ps_e + sh_e
  std::string ps_sh_e_study_string =  TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+FidXCutString +"&&" + W2CutString+"&&"+HCal_Energy_CutString + "&&"+  HCal_Shower_atime_CutString ;
  //// Draw the 2D histogram
  C->Draw("bb_ps_e+bb_sh_e>>ps_sh_e_hist(600, 0, 6)",  ps_sh_e_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *ps_sh_e_hist= (TH1D*)gDirectory->Get("ps_sh_e_hist");
  if (ps_sh_e_hist) {
    ps_sh_e_hist->SetXTitle("Preshower + Shower Energy");
  }

  /// 1d histo of ps_e +sh_e with proton spot cut. 
  std::string proton_ps_sh_e_study_string = TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString  + "&&"+FidXCutString+"&&"+W2CutString+"&&"+HCal_Energy_CutString+"&&"+  HCal_Shower_atime_CutString +"&&" + ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("bb_ps_e+bb_sh_e>>ps_sh_e_hist_proton(600, 0, 6)",  proton_ps_sh_e_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *ps_sh_e_hist_proton= (TH1D*)gDirectory->Get("ps_sh_e_hist_proton");
  if (ps_sh_e_hist_proton) {
    ps_sh_e_hist_proton->SetXTitle("Preshower + Shower Energy");
  }

  /// 1d histo of ps_e+sh_e  with neutron spot cut. 
  std::string neutron_ps_sh_e_study_string =  TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+ "&&"+FidXCutString+"&&" +W2CutString+"&&"+HCal_Energy_CutString+"&&"+  HCal_Shower_atime_CutString +"&&" + NeutronSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("bb_ps_e+bb_sh_e>>ps_sh_e_hist_neutron(600, 0, 6)",  neutron_ps_sh_e_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *ps_sh_e_hist_neutron= (TH1D*)gDirectory->Get("ps_sh_e_hist_neutron");
  if (ps_sh_e_hist_neutron) {
    ps_sh_e_hist_neutron->SetXTitle("Preshower + Shower Energy");
  }

  // Get the neutron/proton ratio, bin-by-bin, by dividing the ps_e+sh_e histo with the neutron spot cut by the ps_e histo with the proton spot cut.
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH1D *np_ps_sh_e_hist = (TH1D*)ps_sh_e_hist_neutron->Clone("np_ps_sh_e_hist"); 
  np_ps_sh_e_hist ->Divide(ps_sh_e_hist_proton);
  double max_ps_sh_e = np_ps_sh_e_hist->GetXaxis()->GetXmax();
  TF1 *ps_sh_e_fit = new TF1("ps_sh_e_fit","[0]",1.7,max_ps_sh_e);
  np_ps_sh_e_hist ->Fit(ps_sh_e_fit,"Q R");
  np_ps_sh_e_hist->SetXTitle("Preshower + Shower Energy");
  np_ps_sh_e_hist->SetYTitle("n/p");
  np_ps_sh_e_hist->GetYaxis()->SetRangeUser(0,1);
  np_ps_sh_e_hist ->Draw("E");
 

  // Studying n/p for grinch_clus_size

  /// 1d histo of grinch_clus_size
  std::string grinch_clus_size_study_string = EnergyCutString+"&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+FidXCutString +"&&" + W2CutString+"&&"+HCal_Energy_CutString + "&&"+  HCal_Shower_atime_CutString ;
  //// Draw the 2D histogram
  C->Draw("bb_grinch_clus_size>>grinch_clus_size_hist(20, 0, 20)",  grinch_clus_size_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *grinch_clus_size_hist= (TH1D*)gDirectory->Get("grinch_clus_size_hist");
  if (grinch_clus_size_hist) {
    grinch_clus_size_hist->SetXTitle("GRINCH cluster size");
  }

  /// 1d histo of grinch_clus_size with proton spot cut. 
  std::string proton_grinch_clus_size_study_string =EnergyCutString+"&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString  + "&&"+FidXCutString+"&&"+W2CutString+"&&"+HCal_Energy_CutString+"&&"+  HCal_Shower_atime_CutString +"&&" + ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("bb_grinch_clus_size>>grinch_clus_size_hist_proton(20, 0, 20)",  proton_grinch_clus_size_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *grinch_clus_size_hist_proton= (TH1D*)gDirectory->Get("grinch_clus_size_hist_proton");
  if (grinch_clus_size_hist_proton) {
    grinch_clus_size_hist_proton->SetXTitle("GRINCH cluster size");
  }

  /// 1d histo of grinch_clus_size  with neutron spot cut. 
  std::string neutron_grinch_clus_size_study_string = EnergyCutString+"&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+ "&&"+FidXCutString+"&&" +W2CutString+"&&"+HCal_Energy_CutString+"&&"+  HCal_Shower_atime_CutString +"&&" + NeutronSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("bb_grinch_clus_size>>grinch_clus_size_hist_neutron(20, 0, 20)",  neutron_grinch_clus_size_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *grinch_clus_size_hist_neutron= (TH1D*)gDirectory->Get("grinch_clus_size_hist_neutron");
  if (grinch_clus_size_hist_neutron) {
    grinch_clus_size_hist_neutron->SetXTitle("GRINCH cluster size");
  }

  // Get the neutron/proton ratio, bin-by-bin, by dividing the grinch_clus_size histo with the neutron spot cut by the grinch_clus_size histo with the proton spot cut.
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH1D *np_grinch_clus_size_hist = (TH1D*)grinch_clus_size_hist_neutron->Clone("np_grinch_clus_size_hist"); 
  np_grinch_clus_size_hist ->Divide(grinch_clus_size_hist_proton);
  double max_grinch_clus_size = np_grinch_clus_size_hist->GetXaxis()->GetXmax();
  TF1 *grinch_clus_size_fit = new TF1("grinch_clus_size_fit","[0]",0,max_grinch_clus_size);
  np_grinch_clus_size_hist ->Fit(grinch_clus_size_fit,"Q R");
  np_grinch_clus_size_hist->SetXTitle("GRINCH cluster size");
  np_grinch_clus_size_hist->SetYTitle("n/p");
  np_grinch_clus_size_hist ->Draw("E");


  // Studying n/p for grinch_clus_adc

  /// 1d histo of grinch_clus_adc
  std::string grinch_clus_adc_study_string = EnergyCutString+"&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+FidXCutString +"&&" + W2CutString+"&&"+HCal_Energy_CutString + "&&"+  HCal_Shower_atime_CutString ;
  //// Draw the 2D histogram
  C->Draw("bb_grinch_clus_adc>>grinch_clus_adc_hist(500, 0, 500)",  grinch_clus_adc_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *grinch_clus_adc_hist= (TH1D*)gDirectory->Get("grinch_clus_adc_hist");
  if (grinch_clus_adc_hist) {
    grinch_clus_adc_hist->SetXTitle("GRINCH cluster ToT");
  }

  /// 1d histo of grinch_clus_adc with proton spot cut. 
  std::string proton_grinch_clus_adc_study_string =EnergyCutString+"&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString  + "&&"+FidXCutString+"&&"+W2CutString+"&&"+HCal_Energy_CutString+"&&"+  HCal_Shower_atime_CutString +"&&" + ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("bb_grinch_clus_adc>>grinch_clus_adc_hist_proton(500, 0, 500)",  proton_grinch_clus_adc_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *grinch_clus_adc_hist_proton= (TH1D*)gDirectory->Get("grinch_clus_adc_hist_proton");
  if (grinch_clus_adc_hist_proton) {
    grinch_clus_adc_hist_proton->SetXTitle("GRINCH cluster ToT");
  }

  /// 1d histo of grinch_clus_adc  with neutron spot cut. 
  std::string neutron_grinch_clus_adc_study_string = EnergyCutString+"&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+ "&&"+FidXCutString+"&&" +W2CutString+"&&"+HCal_Energy_CutString+"&&"+  HCal_Shower_atime_CutString +"&&" + NeutronSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("bb_grinch_clus_adc>>grinch_clus_adc_hist_neutron(500, 0, 500)",  neutron_grinch_clus_adc_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *grinch_clus_adc_hist_neutron= (TH1D*)gDirectory->Get("grinch_clus_adc_hist_neutron");
  if (grinch_clus_adc_hist_neutron) {
    grinch_clus_adc_hist_neutron->SetXTitle("GRINCH cluster ToT");
  }

  // Get the neutron/proton ratio, bin-by-bin, by dividing the grinch_clus_adc histo with the neutron spot cut by the grinch_clus_adc histo with the proton spot cut.
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH1D *np_grinch_clus_adc_hist = (TH1D*)grinch_clus_adc_hist_neutron->Clone("np_grinch_clus_adc_hist"); 
  np_grinch_clus_adc_hist ->Divide(grinch_clus_adc_hist_proton);
  double max_grinch_clus_adc = np_grinch_clus_adc_hist->GetXaxis()->GetXmax();
  TF1 *grinch_clus_adc_fit = new TF1("grinch_clus_adc_fit","[0]",0,max_grinch_clus_adc);
  np_grinch_clus_adc_hist ->Fit(grinch_clus_adc_fit,"Q R");
  np_grinch_clus_adc_hist->SetXTitle("GRINCH cluster ToT");
  np_grinch_clus_adc_hist->SetYTitle("n/p");
  np_grinch_clus_adc_hist ->Draw("E");


  // Studying Quasi-elastic Electron Detection Eff for GRINCH
  
  
  
  // Studying n/p for coin time

  /// 1d histo of coin time
  std::string coin_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+ HCal_Energy_CutString +"&&"+FidXCutString +"&&" + W2CutString ;
  //// Draw the 2D histogram
  C->Draw("hcal_sh_atime_diff>>coin_hist(300, -15, 15)",  coin_study_string.c_str(), "COLZ");
  // Retrieve and customize histogram
  TH1D *coin_hist= (TH1D*)gDirectory->Get("coin_hist");
  if (coin_hist) {
    coin_hist->SetXTitle("hcal time - shower time");
  }

  /// 1d histo of coin time with proton spot cut. 
  std::string proton_coin_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+HCal_Energy_CutString  + "&&"+FidXCutString+"&&"+W2CutString+"&&" + ProtonSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_sh_atime_diff>>coin_hist_proton(300, -15, 15)",  proton_coin_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *coin_hist_proton= (TH1D*)gDirectory->Get("coin_hist_proton");
  if (coin_hist_proton) {
    coin_hist_proton->SetXTitle("hcal time - shower time  with proton spot cut");
  }

  /// 1d histo of coin time  with neutron spot cut. 
  std::string neutron_coin_study_string = EnergyCutString + "&&"+ TrackQualityCutString + "&&"+TargetVertexCutString +"&&" + Optics_CutString+"&&"+HCal_Energy_CutString + "&&"+FidXCutString+"&&" +W2CutString+"&&" + NeutronSpot_CutString;
  //// Draw the 2D histogram
  C->Draw("hcal_sh_atime_diff>>coin_hist_neutron(300, -15, 15)",  neutron_coin_study_string.c_str(), "COLZ E");
  // Retrieve and customize histogram
  TH1D *coin_hist_neutron= (TH1D*)gDirectory->Get("coin_hist_neutron");
  if (coin_hist_neutron) {
    coin_hist_neutron->SetXTitle("hcal time - shower time with neutron spot cut");
  }


  // Get the neutron/proton ratio, bin-by-bin, by dividing the coin histo with the neutron spot cut by the coin histo with the proton spot cut.
  // Root knows how to do this bin-by-bin with the Divide() function.
  /// hist1/hist2  == hist1->Divide(hist2); 
  TH1D *np_coin_hist = (TH1D*)coin_hist_neutron->Clone("np_coin_hist"); 
  np_coin_hist->SetXTitle("coincidence time");
  np_coin_hist->SetYTitle("n/p");
  np_coin_hist ->Divide(coin_hist_proton);
  TF1 *np_coin_fit = new TF1("np_coin_fit","[0]",-10,10);
  np_coin_hist ->Fit(np_coin_fit,"Q R");
  np_coin_hist ->Draw("E");




  // making lines to draw on the 2d histos
  TLine *LineXi = new TLine(hcal_y_exp_min,hcal_x_exp_min, hcal_y_exp_max, hcal_x_exp_min); //horizontal line at hcal_x_exp_min from Yi to Yf.  geometry: (x1, y1, x2, y2)
  LineXi->SetLineWidth(2);
  // LineXi->SetLineStyle(2);
  LineXi ->SetLineColor(kRed);
  TLine *LineXf =  new TLine(hcal_y_exp_min,hcal_x_exp_max, hcal_y_exp_max, hcal_x_exp_max); //horizontal line at hcal_x_exp_max from Yi to Yf.
  LineXf->SetLineWidth(2);
  //LineXf->SetLineStyle(2);
  LineXf ->SetLineColor(kRed);
  TLine *LineXproton =  new TLine(hcal_y_exp_min, hcal_x_exp_min + dx_pn, hcal_y_exp_max, hcal_x_exp_min + dx_pn); //horizontal line at hcal_x_exp_min+dx_pn from Yi to Yf.
  LineXproton->SetLineWidth(2);
  LineXproton->SetLineStyle(2);
  LineXproton ->SetLineColor(kRed);
  TLine *LineYi = new TLine(hcal_y_exp_min,hcal_x_exp_min, hcal_y_exp_min, hcal_x_exp_max); //vert line at hcal_y_exp_min from Xi to Xf.
  LineYi->SetLineWidth(2);
  // LineYi->SetLineStyle(2);
  LineYi ->SetLineColor(kRed);
  TLine *LineYf = new TLine(hcal_y_exp_max,hcal_x_exp_min, hcal_y_exp_max, hcal_x_exp_max); //vert line at hcal_y_exp_max from Xi to Xf.
  LineYf->SetLineWidth(2);
  //LineYf->SetLineStyle(2);
  LineYf ->SetLineColor(kRed);

  

  TCanvas* W2_canvas = new TCanvas("W2_canvas", "W2_canvas", 1000, 600);
  W2_canvas ->Divide(1,2);
  W2_canvas->cd(1);
  W2_hist->SetTitle("W^{2}");
  W2_hist->Draw("E");
  W2_hist_proton->SetLineColor(kGreen);
  W2_hist_proton->Draw("same E");
  W2_hist_neutron->SetLineColor(kMagenta);
  W2_hist_neutron->Draw("same E");
  TLegend *W2_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  W2_leg ->AddEntry(W2_hist,"W^{2} w/ no spot cuts","l");
  W2_leg ->AddEntry(W2_hist_proton,"W^{2} w/ proton spot cut","l");
  W2_leg ->AddEntry(W2_hist_neutron,"W^{2} w/ neutron spot cut","l");
  W2_leg->Draw();
  W2_canvas->cd(2);
  np_W2_hist->Draw("E");

  TCanvas* hcal_x_canvas = new TCanvas("hcal_x_canvas", "hcal_x_canvas", 1000, 600);
  hcal_x_canvas ->Divide(1,2);
  hcal_x_canvas->cd(1);
  gPad->SetLogy();
  hcal_x_hist->SetTitle("hcal x");
  hcal_x_hist->Draw("E");
  hcal_x_hist_proton->SetLineColor(kGreen);
  hcal_x_hist_proton->Draw("same E");
  hcal_x_hist_neutron->SetLineColor(kMagenta);
  hcal_x_hist_neutron->Draw("same E");
  TLegend *hcal_x_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  hcal_x_leg ->AddEntry(hcal_x_hist,"No Spot Cuts","l");
  hcal_x_leg ->AddEntry(hcal_x_hist_proton,"Proton Spot Cut","l");
  hcal_x_leg ->AddEntry(hcal_x_hist_neutron,"Neutron Spot Cut","l");
  hcal_x_leg->Draw();
  hcal_x_canvas->cd(2);
  np_hcal_x_hist->Draw("E");

  TCanvas* hcal_x_exp_canvas = new TCanvas("hcal_x_exp_canvas", "hcal_x_exp_canvas", 1000, 600);
  hcal_x_exp_canvas ->Divide(1,2);
  hcal_x_exp_canvas->cd(1);
  gPad->SetLogy();
  hcal_x_exp_hist->SetTitle("hcal x expected");
  hcal_x_exp_hist->Draw("E");
  hcal_x_exp_hist_proton->SetLineColor(kGreen);
  hcal_x_exp_hist_proton->Draw("same E");
  hcal_x_exp_hist_neutron->SetLineColor(kMagenta);
  hcal_x_exp_hist_neutron->Draw("same E");
  TLegend *hcal_x_exp_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  hcal_x_exp_leg ->AddEntry(hcal_x_exp_hist,"No Spot Cuts","l");
  hcal_x_exp_leg ->AddEntry(hcal_x_exp_hist_proton,"Proton Spot Cut","l");
  hcal_x_exp_leg ->AddEntry(hcal_x_exp_hist_neutron,"Neutron Spot Cut","l");
  hcal_x_exp_leg->Draw();
  hcal_x_exp_canvas->cd(2);
  np_hcal_x_exp_hist->Draw("E");
  TLine* np_hcal_x_exp_min_line = new TLine(hcalposXi_mc,0, hcalposXi_mc, 1);
  np_hcal_x_exp_min_line->SetLineColor(kRed);  // Set line color (e.g., red)
  np_hcal_x_exp_min_line ->SetLineWidth(2);     // Set line width
  //np_hcal_x_exp_min_line->Draw("SAME");
  TLine* np_hcal_x_exp_max_line = new TLine(hcalposXf_mc, 0, hcalposXf_mc, 1 );
  np_hcal_x_exp_max_line->SetLineColor(kRed);  // Set line color (e.g., red)
  np_hcal_x_exp_max_line ->SetLineWidth(2);     // Set line width
  // np_hcal_x_exp_max_line->Draw("SAME");


  TCanvas* HDE_hcal_x_exp_canvas = new TCanvas("HDE_hcal_x_exp_canvas", "HDE_hcal_x_exp_canvas", 1000, 600);
  HDE_hcal_x_exp_canvas ->Divide(1,2);
  HDE_hcal_x_exp_canvas->cd(1);
  gPad->SetLogy();
  expected_hcal_x_exp_hist->SetTitle("hcal x expected");
  expected_hcal_x_exp_hist->Draw("E");
  detected_hcal_x_exp_hist->SetLineColor(kViolet);
  detected_hcal_x_exp_hist->Draw("same E");
  TLegend *HDE_hcal_x_exp_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  HDE_hcal_x_exp_leg ->AddEntry(expected_hcal_x_exp_hist,"No HCal Cuts","l");
  HDE_hcal_x_exp_leg ->AddEntry(detected_hcal_x_exp_hist,"HCal Cuts and Spot Cuts","l");
  HDE_hcal_x_exp_leg->Draw();
  HDE_hcal_x_exp_canvas->cd(2);
  HDE_hcal_x_exp_hist->Draw("E");

  

  
  TCanvas* hcal_y_canvas = new TCanvas("hcal_y_canvas", "hcal_y_canvas", 1000, 600);
  hcal_y_canvas ->Divide(1,2);
  hcal_y_canvas->cd(1);
  gPad->SetLogy();
  hcal_y_hist->SetTitle("hcal y");
  hcal_y_hist->Draw("E");
  hcal_y_hist_proton->SetLineColor(kGreen);
  hcal_y_hist_proton->Draw("same E");
  hcal_y_hist_neutron->SetLineColor(kMagenta);
  hcal_y_hist_neutron->Draw("same E");
  TLegend *hcal_y_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend boxd
  hcal_y_leg ->AddEntry(hcal_y_hist,"No Spot Cuts","l");
  hcal_y_leg ->AddEntry(hcal_y_hist_proton,"Proton Spot Cut","l");
  hcal_y_leg ->AddEntry(hcal_y_hist_neutron,"Neutron Spot Cut","l");
  hcal_y_leg->Draw();
  hcal_y_canvas->cd(2);
  np_hcal_y_hist->Draw("E");

  TCanvas* hcal_y_exp_canvas = new TCanvas("hcal_y_exp_canvas", "hcal_y_exp_canvas", 1000, 600);
  hcal_y_exp_canvas ->Divide(1,2);
  hcal_y_exp_canvas->cd(1);
  gPad->SetLogy();
  hcal_y_exp_hist->SetTitle("hcal y expected");
  hcal_y_exp_hist->Draw("E");
  hcal_y_exp_hist_proton->SetLineColor(kGreen);
  hcal_y_exp_hist_proton->Draw("same E");
  hcal_y_exp_hist_neutron->SetLineColor(kMagenta);
  hcal_y_exp_hist_neutron->Draw("same E");
  TLegend *hcal_y_exp_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  hcal_y_exp_leg ->AddEntry(hcal_y_exp_hist,"No Spot Cuts","l");
  hcal_y_exp_leg ->AddEntry(hcal_y_exp_hist_proton,"Proton Spot Cut","l");
  hcal_y_exp_leg ->AddEntry(hcal_y_exp_hist_neutron,"Neutron Spot Cut","l");
  hcal_y_exp_leg->Draw();
  hcal_y_exp_canvas->cd(2);
  np_hcal_y_exp_hist ->GetYaxis() ->SetRangeUser(0, 1);
  np_hcal_y_exp_hist->Draw("E");


  TCanvas* HDE_hcal_y_exp_canvas = new TCanvas("HDE_hcal_y_exp_canvas", "HDE_hcal_y_exp_canvas", 1000, 600);
  HDE_hcal_y_exp_canvas ->Divide(1,2);
  HDE_hcal_y_exp_canvas->cd(1);
  gPad->SetLogy();
  expected_hcal_y_exp_hist->SetTitle("hcal y expected");
  expected_hcal_y_exp_hist->Draw("E");
  detected_hcal_y_exp_hist->SetLineColor(kViolet);
  detected_hcal_y_exp_hist->Draw("same E");
  TLegend *HDE_hcal_y_exp_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  HDE_hcal_y_exp_leg ->AddEntry(expected_hcal_y_exp_hist,"No HCal Cuts","l");
  HDE_hcal_y_exp_leg ->AddEntry(detected_hcal_y_exp_hist,"HCal Cuts and Spot Cuts","l");
  HDE_hcal_y_exp_leg->Draw();
  HDE_hcal_y_exp_canvas->cd(2);
  HDE_hcal_y_exp_hist->Draw("E");

  
  TCanvas* hcal_dy_canvas = new TCanvas("hcal_dy_canvas", "hcal_dy_canvas", 1000, 600);
  hcal_dy_canvas ->Divide(1,2);
  hcal_dy_canvas->cd(1);
  gPad->SetLogy();
  hcal_dy_hist->SetTitle("hcal dy");
  hcal_dy_hist->Draw("E");
  hcal_dy_hist_proton->SetLineColor(kGreen);
  hcal_dy_hist_proton->Draw("same E");
  hcal_dy_hist_neutron->SetLineColor(kMagenta);
  hcal_dy_hist_neutron->Draw("same E");
  TLegend *hcal_dy_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  hcal_dy_leg ->AddEntry(hcal_dy_hist,"No Spot Cuts","l");
  hcal_dy_leg ->AddEntry(hcal_dy_hist_proton,"Proton Spot Cut","l");
  hcal_dy_leg ->AddEntry(hcal_dy_hist_neutron,"Neutron Spot Cut","l");
  hcal_dy_leg->Draw();
  hcal_dy_canvas->cd(2);
  np_hcal_dy_hist->Draw("E");


  TCanvas* coin_canvas = new TCanvas("coin_canvas", "coin_canvas", 1000, 600);
  coin_canvas ->Divide(1,2);
  coin_canvas->cd(1);
  coin_hist->SetTitle("HCal - SH time");
  coin_hist->Draw("E");
  coin_hist_proton->SetLineColor(kGreen);
  coin_hist_proton->Draw("same E");
  coin_hist_neutron->SetLineColor(kMagenta);
  coin_hist_neutron->Draw("same E");
  TLegend *coin_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  coin_leg ->AddEntry(coin_hist,"No Spot Cuts","l");
  coin_leg ->AddEntry(coin_hist_proton,"Proton Spot Cut","l");
  coin_leg ->AddEntry(coin_hist_neutron,"Neutron Spot Cut","l");
  coin_leg->Draw();
  coin_canvas->cd(2);
  np_coin_hist->Draw("E");

  TCanvas* ps_e_canvas = new TCanvas("ps_e_canvas", "ps_e_canvas", 1000, 600);
  ps_e_canvas ->Divide(1,2);
  ps_e_canvas->cd(1);
  ps_e_hist->SetTitle("Preshower Energy");
  ps_e_hist->Draw("E");
  ps_e_hist_proton->SetLineColor(kGreen);
  ps_e_hist_proton->Draw("same E");
  ps_e_hist_neutron->SetLineColor(kMagenta);
  ps_e_hist_neutron->Draw("same E");
  TLegend *ps_e_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  ps_e_leg ->AddEntry(ps_e_hist,"No Spot Cuts","l");
  ps_e_leg ->AddEntry(ps_e_hist_proton,"Proton Spot Cut","l");
  ps_e_leg ->AddEntry(ps_e_hist_neutron,"Neutron Spot Cut","l");
  ps_e_leg->Draw();
  ps_e_canvas->cd(2);
  np_ps_e_hist->Draw("E");

  TCanvas* ps_sh_e_canvas = new TCanvas("ps_sh_e_canvas", "ps_sh_e_canvas", 1000, 600);
  ps_sh_e_canvas ->Divide(1,2);
  ps_sh_e_canvas->cd(1);
  ps_sh_e_hist->SetTitle("Preshower + Shower Energy");
  ps_sh_e_hist->Draw("E");
  ps_sh_e_hist_proton->SetLineColor(kGreen);
  ps_sh_e_hist_proton->Draw("same E");
  ps_sh_e_hist_neutron->SetLineColor(kMagenta);
  ps_sh_e_hist_neutron->Draw("same E");
  TLegend *ps_sh_e_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  ps_sh_e_leg ->AddEntry(ps_sh_e_hist,"No Spot Cuts","l");
  ps_sh_e_leg ->AddEntry(ps_sh_e_hist_proton,"Proton Spot Cut","l");
  ps_sh_e_leg ->AddEntry(ps_sh_e_hist_neutron,"Neutron Spot Cut","l");
  ps_sh_e_leg->Draw();
  ps_sh_e_canvas->cd(2);
  np_ps_sh_e_hist->Draw("E");
  
  TCanvas* grinch_clus_size_canvas = new TCanvas("grinch_clus_size_canvas", "grinch_clus_size_canvas", 1000, 600);
  grinch_clus_size_canvas ->Divide(1,2);
  grinch_clus_size_canvas->cd(1);
  grinch_clus_size_hist->SetTitle("GRINCH cluster size");
  grinch_clus_size_hist->Draw("E");
  grinch_clus_size_hist_proton->SetLineColor(kGreen);
  grinch_clus_size_hist_proton->Draw("same E");
  grinch_clus_size_hist_neutron->SetLineColor(kMagenta);
  grinch_clus_size_hist_neutron->Draw("same E");
  TLegend *grinch_clus_size_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  grinch_clus_size_leg ->AddEntry(grinch_clus_size_hist,"No Spot Cuts","l");
  grinch_clus_size_leg ->AddEntry(grinch_clus_size_hist_proton,"Proton Spot Cut","l");
  grinch_clus_size_leg ->AddEntry(grinch_clus_size_hist_neutron,"Neutron Spot Cut","l");
  grinch_clus_size_leg->Draw();
  grinch_clus_size_canvas->cd(2);
  np_grinch_clus_size_hist->Draw("E");


  TCanvas* grinch_clus_adc_canvas = new TCanvas("grinch_clus_adc_canvas", "grinch_clus_adc_canvas", 1000, 600);
  grinch_clus_adc_canvas ->Divide(1,2);
  grinch_clus_adc_canvas->cd(1);
  grinch_clus_adc_hist->SetTitle("GRINCH cluster ToT");
  grinch_clus_adc_hist->Draw("E");
  grinch_clus_adc_hist_proton->SetLineColor(kGreen);
  grinch_clus_adc_hist_proton->Draw("same E");
  grinch_clus_adc_hist_neutron->SetLineColor(kMagenta);
  grinch_clus_adc_hist_neutron->Draw("same E");
  TLegend *grinch_clus_adc_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  grinch_clus_adc_leg ->AddEntry(grinch_clus_adc_hist,"No Spot Cuts","l");
  grinch_clus_adc_leg ->AddEntry(grinch_clus_adc_hist_proton,"Proton Spot Cut","l");
  grinch_clus_adc_leg ->AddEntry(grinch_clus_adc_hist_neutron,"Neutron Spot Cut","l");
  grinch_clus_adc_leg->Draw();
  grinch_clus_adc_canvas->cd(2);
  np_grinch_clus_adc_hist->Draw("E");

  TCanvas* hcal_e_canvas = new TCanvas("hcal_e_canvas", "hcal_e_canvas", 1000, 600);
  hcal_e_canvas ->Divide(1,2);
  hcal_e_canvas->cd(1);
  hcal_e_hist->SetTitle("HCal Energy");
  hcal_e_hist->Draw("E");
  hcal_e_hist_proton->SetLineColor(kGreen);
  hcal_e_hist_proton->Draw("same E");
  hcal_e_hist_neutron->SetLineColor(kMagenta);
  hcal_e_hist_neutron->Draw("same E");
  TLegend *hcal_e_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  hcal_e_leg ->AddEntry(hcal_e_hist,"No Spot Cuts","l");
  hcal_e_leg ->AddEntry(hcal_e_hist_proton,"Proton Spot Cut","l");
  hcal_e_leg ->AddEntry(hcal_e_hist_neutron,"Neutron Spot Cut","l");
  hcal_e_leg->Draw();
  hcal_e_canvas->cd(2);
  np_hcal_e_hist->Draw("E");

  TCanvas* e_over_p_canvas = new TCanvas("e_over_p_canvas", "e_over_p_canvas", 1000, 600);
  e_over_p_canvas ->Divide(1,2);
  e_over_p_canvas->cd(1);
  e_over_p_hist->SetTitle("E over P");
  e_over_p_hist->Draw("E");
  e_over_p_hist_proton->SetLineColor(kGreen);
  e_over_p_hist_proton->Draw("same E");
  e_over_p_hist_neutron->SetLineColor(kMagenta);
  e_over_p_hist_neutron->Draw("same E");
  TLegend *e_over_p_leg = new TLegend(0.7, 0.5, 0.9, 0.7);// (x1, y1, x2, y2) are the coordinates of the legend box
  e_over_p_leg ->AddEntry(e_over_p_hist,"No Spot Cuts","l");
  e_over_p_leg ->AddEntry(e_over_p_hist_proton,"Proton Spot Cut","l");
  e_over_p_leg ->AddEntry(e_over_p_hist_neutron,"Neutron Spot Cut","l");
  e_over_p_leg->Draw();
  e_over_p_canvas->cd(2);
  np_e_over_p_hist->Draw("E");


  
  
  TCanvas* results_canvas = new TCanvas("results_canvas", "results_canvas", 1000, 600);
  results_canvas ->Divide(2,2);
  results_canvas ->cd(1);
  np_W2_hist->GetYaxis() ->SetRangeUser(0, 1);
  np_W2_hist->SetTitle("n:p ratio across W^{2}");
  np_W2_hist->Draw("E");
  results_canvas ->cd(2);
  np_hcal_x_hist->GetYaxis() ->SetRangeUser(0, 1);
  np_hcal_x_hist->SetTitle("n:p ratio across hcal x");
  np_hcal_x_hist ->Draw("E");
  results_canvas->cd(3);
  np_hcal_y_hist->GetYaxis() ->SetRangeUser(0, 1);
  np_hcal_y_hist->SetTitle("n:p ratio across hcal y");
  np_hcal_y_hist->Draw("E");
  results_canvas->cd(4);
  //np_coin_hist->GetYaxis() ->SetRangeUser(0, 2);
  np_coin_hist->SetTitle("n:p ratio across Coincidence Time");
  np_coin_hist->Draw("E");

  TCanvas* results_canvas2 = new TCanvas("results_canvas2", "results_canvas2", 1000, 600);
  results_canvas2 ->Divide(2,2);
  adjustCanvas(results_canvas2);
  results_canvas2 ->cd(1);
  np_e_over_p_hist->GetYaxis() ->SetRangeUser(0, 1);
  np_e_over_p_hist->SetTitle("n:p ratio across e/p");
  np_e_over_p_hist->Draw("E");
  results_canvas2 ->cd(2);
  np_hcal_e_hist->GetYaxis() ->SetRangeUser(0, 1);
  np_hcal_e_hist->SetTitle("n:p ratio across HCal Energy");
  np_hcal_e_hist->Draw("E");
  results_canvas2 ->cd(3);
  np_ps_e_hist->GetYaxis() ->SetRangeUser(0, 1);
  np_ps_e_hist->SetTitle("n:p ratio across Preshower Energy");
  np_ps_e_hist->Draw("E");
  results_canvas2 ->cd(4);
  np_hcal_dy_hist->GetYaxis() ->SetRangeUser(0, 1);
  np_hcal_dy_hist->SetTitle("n:p ratio across hcal dy");
  np_hcal_dy_hist->Draw("E");

  TCanvas* results_canvas3 = new TCanvas("results_canvas3", "results_canvas3", 1000, 600);
  results_canvas3 ->Divide(2,2);
  results_canvas3 ->cd(1);
  np_hcal_y_exp_hist->GetYaxis() ->SetRangeUser(0, 1);
  np_hcal_y_exp_hist->SetTitle("n:p ratio across hcal y expected");
  np_hcal_y_exp_hist->Draw("E");
  // hcal_x_exp_hcal_y_exp_hist_neutron->SetTitle("netutron");
  // hcal_x_exp_hcal_y_exp_hist_neutron ->Draw("colz");
  results_canvas3 ->cd(2);
  np_hcal_x_exp_hist->GetYaxis() ->SetRangeUser(0, 1);
  np_hcal_x_exp_hist->SetTitle("n:p ratio across hcal x Expected");
  np_hcal_x_exp_hist->Draw("E");
  //hcal_x_exp_hcal_y_exp_hist_proton->SetTitle("proton");
  // hcal_x_exp_hcal_y_exp_hist_proton ->Draw("colz");
  results_canvas3 ->cd(3);
  np_hcal_x_exp_hcal_y_exp_hist->SetTitle("n:p ratio across hcal x exp and hcal y exp");
  np_hcal_x_exp_hcal_y_exp_hist->Draw("colz");
  
  
  TCanvas* results_canvas4 = new TCanvas("results_canvas4", "results_canvas4", 1000, 600);
  results_canvas4 ->Divide(2,2);
  results_canvas4 ->cd(1);
  HDE_hcal_y_exp_hist->GetYaxis() ->SetRangeUser(0, 1);
  HDE_hcal_y_exp_hist->SetTitle("Nucleon Detection Efficiency across hcal_y_exp");
  HDE_hcal_y_exp_hist->Draw("E");
  results_canvas4 ->cd(2);
  HDE_hcal_x_exp_hist->GetYaxis() ->SetRangeUser(0, 1);
  HDE_hcal_x_exp_hist->SetTitle("Nucleon Detection Efficiency across hcal_x_exp");
  HDE_hcal_x_exp_hist->Draw("E");
  results_canvas4->cd(3);
  HDE_hcal_x_exp_hcal_y_exp_hist->SetTitle("Nucleon Detection Efficiency across hcal_x_exp and hcal_y_exp");
  HDE_hcal_x_exp_hcal_y_exp_hist->Draw("Colz");
  LineXi->Draw();
  LineXf->Draw();
  LineYi->Draw();
  LineYf->Draw();

  TCanvas* results_canvas5 = new TCanvas("results_canvas5", "results_canvas5", 1000, 600);
  results_canvas5 ->Divide(2,2);
  results_canvas5 ->cd(1);
  np_vz_hist->GetYaxis() ->SetRangeUser(0, 1);
  np_vz_hist->SetTitle("n:p ratio across track vz");
  np_vz_hist->Draw("E");

 
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
