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

const double hcalheight = -0.2897;
//function declarations


void dxdydraft(const char *configfilename = "setup_gmn.cfg",const char *outputfilename="gmn_out.root"){//main
  
  gStyle->SetNumberContours(255); 

  // TChain *fchain = new TChain("T");
  // while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ) {
  //    if( !currentline.BeginsWith("#") ){
  //     cout << " add file : " << currentline << endl;
  //     fchain->Add(currentline);
  //   }   
  // } 

  TChain *C = new TChain("T");
  // C->Add("/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS4/LD2/rootfiles/e1209019_fullreplay_11595*");
 //  C->Add("/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS4/LD2/rootfiles/e1209019_fullreplay_11589*");
 //  C->Add("/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS4/LD2/rootfiles/e1209019_fullreplay_11590*");
 //  C->Add("/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS4/LD2/rootfiles/e1209019_fullreplay_11592*");
 // C->Add("/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS4/LD2/rootfiles/e1209019_fullreplay_11593*");
 // C->Add("/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS4/LD2/rootfiles/e1209019_fullreplay_11594*");
 // C->Add("/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS4/LD2/rootfiles/e1209019_fullreplay_11595*");
 // C->Add("/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS4/LD2/rootfiles/e1209019_fullreplay_11597*");
  //C->Add("/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS4/LD2/rootfiles/e1209019_fullreplay_1159*");
  //SBS 4
  double E_e = 3.728;
  double BB_d = 1.7988;
  double BB_th = 36.0* TMath::DegToRad();
  double HCal_d =11.0;
  double HCal_th = 31.9 * TMath::DegToRad();
  double W2_mean = 1.00;
  double W2_sig = 0.24;

   cout<<"kine: "<<kine<<endl;

  string configfilename = Form("config/sbs%d.cfg",kine);

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

  cout<<"Applying Global Cut to rootfiles. May take a few minutes....";
  TEventList *elist = new TEventList("elist","Elastic Event List");
  C->Draw(">>elist",globalcut);

  double BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  double BBtr_vz[maxTracks];
  double BBtgt_x[maxTracks], BBtgt_y[maxTracks], BBtgt_th[maxTracks], BBtgt_ph[maxTracks];
  double BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;
  double HCALx, HCALy, HCALe, ekineW2;


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



  TFile *fout = new TFile(outputfilename,"RECREATE");

  TH1D *hdx_HCAL = new TH1D("hdx_HCAL ", " ; x_{HCAL} - x_{exp} (m)  ", 250,-4,3);
  TH1D *hdy_HCAL = new TH1D("hdy_HCAL ", " ; y_{HCAL} - y_{exp} (m)  ", 250,-2,2);
  TH2D *hdxdy_HCAL = new TH2D("hdxdy_HCAL","HCAL_{actual} - HCAL_{expected} positions  ;y_{HCAL}-y_{expect} (m); x_{HCAL}-x_{expect} (m)", 250, -1.25, 1.25, 250, -2.5, 2.5 );

  TH1D* h_W2 = new TH1D("h_W2", "W2",100,0,2);

  Long64_t Nevents = elist ->GetN();
  

  cout<< "# Events = "<<Nevents<<endl;
  for( Long64_t nevent = 1; nevent <Nevents; nevent++){

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

    hdx_HCAL ->Fill(dx);
    hdy_HCAL ->Fill(dy);
    hdxdy_HCAL ->Fill(dy,dx);
    h_W2 ->Fill(ekineW2);
    
  }



  fout->Write();

}// end main







