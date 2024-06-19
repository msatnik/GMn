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
using namespace std;


//global variables
const int maxTracks = 16; 
const double PI = TMath::Pi();
const double M_e = 0.0051;
const double M_p = 0.938272;
const double M_n = 0.939565;

//const double hcalheight = -0.2897; // Need for Pass0 or Pass1 data
//const double hcalheight = 0; // should put in config file maybe

 // HCAL 
  // Dimensions
  static const Int_t hcalchan = 288;
  static const Int_t hcalcol = 12;
  static const Int_t hcalrow = 24;
  static const Double_t hcal_hrange = 1.85928;    //m, total range in horizontal direction of HCal (end-to-end)
  static const Double_t hcal_vrange = 3.81;       //m, total range in vertical direction of HCal (end-to-end)
  static const Double_t hcaladc_binw = 4.;        //ns, width of each ADC bin
  static const Double_t hcalblk_w_p0 = 0.15;      //m, width of a HCAL block, pass0/1
  static const Double_t hcalblk_h_p0 = 0.15;      //m, height of a HCAL block, pass0/1 
  static const Double_t hcalblk_w = 0.1524;       //m, width of a HCAL block
  static const Double_t hcalblk_h = 0.1524;       //m, height of a HCAL block
  static const Double_t hcalblk_div_h = 0.15494;  //m, horizontal center-to-center dist.
  static const Double_t hcalblk_div_v = 0.15875;  //m, vertical center-to-center dist.
  static const Double_t hcalblk_div_hyp = 0.15875;//m, division corner-to-cornter dist.
  static const Double_t hcalblk_gap_h = 0.00254;  //m, horiz. gap bet. two blocks
  static const Double_t hcalblk_gap_v = 0.00635;  //m, vert. gap bet. two blocks
  // Positions (mc)
  static const Double_t hcalposXi_mc = -2.655;    //m, distance from beam center to top of HCal w/75cm offset
  static const Double_t hcalposXf_mc = 1.155;     //m, distance from beam center to bottom of HCal w/75cm offset
  static const Double_t hcalposYi_mc = -0.92964;  //m, distance from beam center to opposite-beam side of HCal
  static const Double_t hcalposYf_mc = 0.92964;   //m, distance from beam center to beam side of HCal
  // Pass0/1 (no block spacing)
  static const Double_t hcalposXi_p0 = -2.16014;  //m, distance from beam center to top of HCal w/75cm offset
  static const Double_t hcalposXf_p0 = 1.43826;   //m, distance from beam center to bottom of HCal w/75cm offset
  static const Double_t hcalposYi_p0 = -0.9;      //m, distance from beam center to opposite-beam side of HCal
  static const Double_t hcalposYf_p0 = 0.9;       //m, distance from beam center to beam side of HCal
  // Positions (data fits)
  static const Double_t hcalposXi = -2.268095;    //m, distance from beam center to top of HCal (obsolete)
  static const Double_t hcalposXf = 1.538095;     //m, distance from beam center to bottom of HCal (obsolete)
  static const Double_t hcalposYi = -0.931545;    //m, distance from beam center to opposite-beam side of HCal (obsolete)
  static const Double_t hcalposYf = 0.931545;     //m, distance from beam center to beam side of HCal (obsolete)
  // Global
  //static const Double_t hcalvoff = -0.2897;       //m, height of the center of hcal above beam (m)
  //static const Double_t hcalvoff = -0.3735;       //m, height of the center of hcal above beam (m) (sbs8)
  //static const Double_t hcalvoff = -0.75;       //m, height of the center of hcal above beam (m)
  //static const Double_t hcalvoff = 0.0;         //m, height of the center of hcal above beam (m) after pass2 corr


//function declarations
std::vector<Double_t> hcalaa_data (int exblkN_x=1, int exblkN_y=1);
bool hcalaaON (Double_t hcalx, Double_t hcaly, std::vector<Double_t> hcalaa);
std::vector<Double_t> hcalfid (Double_t dxsig_p, Double_t dxsig_n, Double_t dysig, std::vector<Double_t> hcalaa);
bool hcalfidIN (Double_t hcalx_exp, Double_t hcaly_exp, Double_t dx_pn, vector<Double_t> fid);
std::vector<Double_t> hcalaa_mc (int exblkN_x=1, int exblkN_y=1);


void dxdy_simc_v2(int entries_input = -1, int kine = 4 ){//main
  
  gStyle->SetNumberContours(255); 


  TChain *C = new TChain("T");
 
  //SBS 4
  double E_e = 3.728;
  double BB_d = 1.7988;
  double BB_th = 36.0* TMath::DegToRad();
  double HCal_d =11.0;
  double HCal_th = 31.9 * TMath::DegToRad();
  double W2_mean = 1.00;
  double W2_sig = 0.24;
  //const double hcalheight = -0.2897; // Need for Pass0 or Pass1 data
  double hcalheight = 0; // should put in config file maybe
  double dxsig_p = 0.194; // sbs4 30% data 
  double dxsig_n = 0.21; // sbs4 30% data 
  double dysig = 0.3;
  double dx_pn = 0.639;


  cout<<"kine: "<<kine<<endl;
  //cout<<"config file"<<configfileinput<<endl;

  //// set location to config file
   //string configfilename = Form("../config/sbs%d.cfg",kine);
   //string configfilename = Form("../config/%s.cfg",configfileinput);
   string configfilename = "/w/halla-scshelf2102/sbs/msatnik/GMn/config/sbs4_simc_30p_deen.cfg";

   // Declare outfile
   // TFile *fout = new TFile( Form("../output/sbs%d.root",kine), "RECREATE" );
  TFile *fout = new TFile("../output/sbs4_simc_30p_deep.root", "RECREATE" );
  //TFile *fout = new TFile(outputfilename,"RECREATE");


  cout<<"testing reading in from a config file"<<endl;

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
      if( skey == "HCal_d" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_d = sval.Atof();
	cout << "Loading HCal distance: " << HCal_d << endl;
      }
      if( skey == "HCal_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_th = sval.Atof();	
	cout << "Loading HCal angle: " << HCal_th << endl;
	HCal_th = HCal_th * TMath::DegToRad();
      }
      if( skey == "hcalheight" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	hcalheight = sval.Atof();	
	cout << "Loading HCal 'height': " << hcalheight<< endl;
      }
      if( skey == "BB_d" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	BB_d = sval.Atof();
	cout << "Loading BB distance: " << BB_d << endl;
      }
      if( skey == "BB_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	BB_th = sval.Atof();	
	cout << "Loading BB angle: " << BB_th << endl;
	BB_th = BB_th * TMath::DegToRad();
      }
      if( skey == "W2_mean" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W2_mean = sval.Atof();
	cout << "Loading W2 mean cut: " << W2_mean << endl;
      }
      if( skey == "W2_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W2_sig = sval.Atof();
	cout << "Loading W2 sigma cut: " << W2_sig << endl;
      }
    }
    delete tokens;
  }

 

  //TCut globalcut = "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&sbs.hcal.e>0.025&&bb.ps.e+bb.sh.e>1.7";

  cout<<"Applying Global Cut to rootfiles. May take a few minutes...."<<endl;
  TEventList *elist = new TEventList("elist","Elastic Event List");
  C->Draw(">>elist",globalcut);



  double BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  double BBtr_vz[maxTracks];
  double BBtgt_x[maxTracks], BBtgt_y[maxTracks], BBtgt_th[maxTracks], BBtgt_ph[maxTracks];
  double BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;
  double HCALx, HCALy, HCALe, ekineW2;

  double mc_sigma, mc_nucl;
  double simc_weight;

  C ->SetBranchStatus("*",0);

  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.p",1);
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("bb.ps.x",1);
  C->SetBranchStatus("bb.ps.y",1);
  C->SetBranchStatus("bb.sh.e",1);
  C->SetBranchStatus("bb.sh.x",1);
  C->SetBranchStatus("bb.sh.y",1);
  C->SetBranchStatus("e.kine.W2",1);

  // MC branches
  C->SetBranchStatus("MC.mc_sigma",1);
  C->SetBranchStatus("MC.mc_nucl",1);
  C->SetBranchStatus("MC.simc_Weight",1);


  C->SetBranchAddress("sbs.hcal.x", &HCALx);
  C->SetBranchAddress("sbs.hcal.y", &HCALy);
  C->SetBranchAddress("sbs.hcal.e", &HCALe);
  C->SetBranchAddress("bb.tr.n", &BBtr_n);
  C->SetBranchAddress("bb.tr.px", BBtr_px);
  C->SetBranchAddress("bb.tr.py", BBtr_py);
  C->SetBranchAddress("bb.tr.pz", BBtr_pz);
  C->SetBranchAddress("bb.tr.p", BBtr_p);
  C->SetBranchAddress("bb.tr.vz", BBtr_vz);
  C->SetBranchAddress("bb.ps.e", &BBps_e);
  C->SetBranchAddress("bb.ps.x", &BBps_x);
  C->SetBranchAddress("bb.ps.y", &BBps_y);
  C->SetBranchAddress("bb.sh.e", &BBsh_e);
  C->SetBranchAddress("bb.sh.x", &BBsh_x);
  C->SetBranchAddress("bb.sh.y", &BBsh_y);
  C->SetBranchAddress("e.kine.W2", &ekineW2);

  // MC branches
  C->SetBranchAddress("MC.mc_sigma", &mc_sigma);
  C->SetBranchAddress("MC.mc_nucl", &mc_nucl);
  C->SetBranchAddress("MC.simc_Weight", &simc_weight);
 



  TH1D *h_dx_HCAL = new TH1D("h_dx_HCAL ", " ; x_{HCAL} - x_{exp} (m)  ", 250,-4,3);
  TH1D *h_dy_HCAL = new TH1D("h_dy_HCAL ", " ; y_{HCAL} - y_{exp} (m)  ", 250,-2,2);
  TH2D *h_dxdy_HCAL = new TH2D("h_dxdy_HCAL","HCAL_{actual} - HCAL_{expected} positions  ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -2.5, 2.5 );
  TH1D* h_W2 = new TH1D("h_W2", "W2",100,0,2);
  TH1D *h_dx_HCAL_p = new TH1D("h_dx_HCAL_p ", " ; x_{HCAL} - x_{exp} (m)  for proton only", 250,-4,3);
  TH1D *h_dx_HCAL_n = new TH1D("h_dx_HCAL_n ", " ; x_{HCAL} - x_{exp} (m)  for neutron only", 250,-4,3);
 TH1D *h_dx_HCAL_noFid = new TH1D("h_dx_HCAL_noFid", " ; x_{HCAL} - x_{exp} (m)  ", 250,-4,3);
  TH1D *h_dy_HCAL_noFid = new TH1D("h_dy_HCAL_noFid", " ; y_{HCAL} - y_{exp} (m)  ", 250,-2,2);
  TH2D *h_dxdy_HCAL_noFid = new TH2D("h_dxdy_HCAL_noFid","HCAL_{actual} - HCAL_{expected} positions  ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -2.5, 2.5 );
  TH1D* h_W2_noFid = new TH1D("h_W2_noFid", "W2",100,0,2);

  TH1D *h_dx_HCAL_noAA = new TH1D("h_dx_HCAL_noAA", " ; x_{HCAL} - x_{exp} (m)  ", 250,-4,3);
  TH1D *h_dy_HCAL_noAA = new TH1D("h_dy_HCAL_noAA", " ; y_{HCAL} - y_{exp} (m)  ", 250,-2,2);
  TH2D *h_dxdy_HCAL_noAA = new TH2D("h_dxdy_HCAL_noAA","HCAL_{actual} - HCAL_{expected} positions  ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -2.5, 2.5 );
  TH1D* h_W2_noAA = new TH1D("h_W2_noAA", "W2",100,0,2);
 

  //TH1D* h_mc_sigma = new TH1D("h_mc_sigma","mc_sigma",


//// define hcal active area
  std::vector<Double_t> hcalaa = hcalaa_mc(1, 1);//
  cout<< "HCal acitve area cut params: "<< hcalaa[0] << ", " << hcalaa[1] << ", " << hcalaa[2] << ", " << hcalaa[3] <<endl;

  //// define hcal fiducial cut 
  std::vector<Double_t> hcalfid_mc =hcalfid(dxsig_p, dxsig_n, dysig, hcalaa);
  cout<< "HCal fiducial cut params: "<< hcalfid_mc[0] << ", " << hcalfid_mc[1] << ", " << hcalfid_mc[2] << ", " << hcalfid_mc[3] <<endl;


  Long64_t Nevents = elist ->GetN();

  Int_t max = 0;
    if (entries_input == -1){
      max = Nevents;
    }
    else{
      max = entries_input;
    }
    if (max > Nevents){ max = Nevents;}
    cout<<"max = "<<max<<endl;

 cout<<endl << "Opened up TChain with nentries: " << C->GetEntries() << ". After globalcut: " << Nevents << "." << endl << endl;



  cout<< "# Events = "<<Nevents<<endl;
  for( Long64_t nevent = 1; nevent <max; nevent++){

    if (nevent%50000==0) cout << " Entry = " << nevent << endl;
    
    C->GetEntry(elist ->GetEntry(nevent)); 
    
    double etheta = acos( BBtr_pz[0]/BBtr_p[0]);
    double ephi = atan2(BBtr_py[0],BBtr_px[0]);
    
    TVector3 vertex(0,0,BBtr_vz[0]);
    TLorentzVector Pbeam(0,0,E_e,E_e);
    TLorentzVector kprime(BBtr_px[0], BBtr_py[0], BBtr_pz[0], BBtr_p[0]);
    TLorentzVector Ptarg(0, 0, 0, M_p);

    TLorentzVector q = Pbeam - kprime;
    TLorentzVector PgammaN = Ptarg +q; //should go through and write this out. Momentum of virtual photon

    double pel = E_e/ (1. +E_e/M_p*(1.-cos(etheta)));//momentum of elastically scattered electron 
    double nu = E_e -BBtr_p[0]; //kinetic energy of the elasticlly scattered electron 
    double pp = sqrt(pow(nu,2)+2 *M_p*nu); 
    double phinucleon = ephi + PI; //coplanar 
    double thetanucleon = acos((E_e - BBtr_pz[0])/pp);

    TVector3 pNhat( sin(thetanucleon) *cos(phinucleon),sin(thetanucleon)*sin(phinucleon),cos(thetanucleon));
    
    TVector3 HCAL_zaxis(sin(-HCal_th),0,cos(-HCal_th));
    TVector3 HCAL_xaxis(0,-1,0);
    TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();

    TVector3 HCAL_origin = HCal_d * HCAL_zaxis + hcalheight * HCAL_xaxis;

    double sintersect = ( HCAL_origin - vertex).Dot( HCAL_zaxis) / (pNhat.Dot( HCAL_zaxis));
    
    TVector3 HCAL_intersect = vertex + sintersect * pNhat; 

    double yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis);
    double xexpect_HCAL =  (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis);

    double Q2 = 2*E_e *BBtr_p[0]*( 1-BBtr_pz[0]/BBtr_p[0]);
    
    double W = PgammaN.M();
    double W2 = ekineW2;

    double dx = HCALx - xexpect_HCAL;
    double dy = HCALy - yexpect_HCAL;

    if (abs(ekineW2- W2_mean)>W2_sig ) continue; 

    h_dx_HCAL_noAA ->Fill(dx, simc_weight);
    h_dy_HCAL_noAA ->Fill(dy, simc_weight);
    h_dxdy_HCAL_noAA ->Fill(dy,dx);
    h_W2_noAA ->Fill(ekineW2);

      //// check if it's on the hcal active area
    if ( hcalaaON( HCALx , HCALy , hcalaa))
      {
	h_dx_HCAL_noFid ->Fill(dx, simc_weight);
	h_dy_HCAL_noFid ->Fill(dy, simc_weight);
	h_dxdy_HCAL_noFid ->Fill(dy,dx);
	h_W2_noFid->Fill(ekineW2);

	//// check if it's within the fiducial cuts 
	if(  hcalfidIN(xexpect_HCAL, xexpect_HCAL, dx_pn, hcalfid_mc ) ) 
	  {
	    h_dx_HCAL ->Fill(dx, simc_weight);
	    h_dy_HCAL ->Fill(dy, simc_weight);
	    h_dxdy_HCAL ->Fill(dy,dx);
	    h_W2->Fill(ekineW2);
	  }//end fid cut
      }//end AA cut


  }


  fout->Write();

}// end main


 // Establish hcal active area excluding N blks from edge, Pass0/1 DB
  std::vector<Double_t> hcalaa_data(int exblkN_x=1, int exblkN_y=1) {
    std::vector<Double_t> hcalaa;
    Double_t hcalaaXi = hcalposXi_p0 + exblkN_x*hcalblk_w_p0;
    Double_t hcalaaXf = hcalposXf_p0 - exblkN_x*hcalblk_w_p0;
    Double_t hcalaaYi = hcalposYi_p0 + exblkN_y*hcalblk_h_p0;
    Double_t hcalaaYf = hcalposYf_p0 - exblkN_y*hcalblk_h_p0;
    hcalaa.push_back( hcalaaXi ); 
    hcalaa.push_back( hcalaaXf );
    hcalaa.push_back( hcalaaYi );
    hcalaa.push_back( hcalaaYf );
    return hcalaa;
  }

 // Check position per event against hcal active area (TRUE if detection on active area)
  bool hcalaaON (Double_t hcalx, Double_t hcaly, std::vector<Double_t> hcalaa) {
    bool on = false;
    // active area dimensions
    Double_t hcalx_t = hcalaa[0];
    Double_t hcalx_b = hcalaa[1];
    Double_t hcaly_r = hcalaa[2];
    Double_t hcaly_l = hcalaa[3];
    on = hcaly>hcaly_r && hcaly<hcaly_l && hcalx>hcalx_t && hcalx<hcalx_b;
    return on;
  } 


// Overload for most cases where dy margin negligable and dx top/bottom configurable by p/n
  std::vector<Double_t> hcalfid (Double_t dxsig_p, Double_t dxsig_n, Double_t dysig, std::vector<Double_t> hcalaa) {
    std::vector<Double_t> fid;
    Double_t hcalx_t = hcalaa[0] + dxsig_p;  // top margin (relevant for proton)
    Double_t hcalx_b = hcalaa[1] - dxsig_n;  // bottom margin (relevant for neutron)
    Double_t hcaly_r = hcalaa[2] + dysig;    // right margin
    Double_t hcaly_l = hcalaa[3] - dysig;    // left margin
    fid.push_back( hcalx_t ); 
    fid.push_back( hcalx_b );
    fid.push_back( hcaly_r );
    fid.push_back( hcaly_l );
    return fid;
  }


  // Check position by event and verify in hcal fiducial area
  bool hcalfidIN (Double_t hcalx_exp, Double_t hcaly_exp, Double_t dx_pn, vector<Double_t> fid) {
    Double_t hcalx_t = fid[0];
    Double_t hcalx_b = fid[1];
    Double_t hcaly_r = fid[2];
    Double_t hcaly_l = fid[3];

    Double_t hcalx_exp_p = hcalx_exp - dx_pn;      //define the exp pos of a proton from obs dx peak diff

    bool infid = hcaly_exp>hcaly_r && hcaly_exp<hcaly_l &&      //dy same for protons and neutrons
	         hcalx_exp>hcalx_t && hcalx_exp<hcalx_b &&      //dx for neutrons
	         hcalx_exp_p>hcalx_t && hcalx_exp_p<hcalx_b;	//dx for protons							     
    return infid;
  } 


 // Establish hcal active area excluding N blks from edge, MC DB
  std::vector<Double_t> hcalaa_mc (int exblkN_x=1, int exblkN_y=1) {
    std::vector<Double_t> hcalaa;
    Double_t hcalaaXi = hcalposXi_mc + exblkN_x*hcalblk_div_v;
    Double_t hcalaaXf = hcalposXf_mc - exblkN_x*hcalblk_div_v;
    Double_t hcalaaYi = hcalposYi_mc + exblkN_y*hcalblk_div_h;
    Double_t hcalaaYf = hcalposYf_mc - exblkN_y*hcalblk_div_h;
    hcalaa.push_back( hcalaaXi ); 
    hcalaa.push_back( hcalaaXf );
    hcalaa.push_back( hcalaaYi );
    hcalaa.push_back( hcalaaYf );
    return hcalaa;
  }


