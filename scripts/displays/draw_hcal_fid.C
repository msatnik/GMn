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


//// default values to be overwritten 
TString DataFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_tight_2Dhistos.root";
TString ProtonFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_tight_2Dhistos.root";
TString NeutronFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_tight_2Dhistos.root";
TString InelasticFileString = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_inel_2Dhistos.root";
double dx_pn = 0.738;
double sigma_p =0.2549;
double sigma_n = 0.2523;
double sigma_dy =0.2813;

// SBS4 0p 
TString DataFileString_sbs4_0p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_2Dhistos_sept26.root";
TString ProtonFileString_sbs4_0p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_2Dhistos_sept26.root";
TString NeutronFileString_sbs4_0p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_2Dhistos_sept26.root";
TString InelasticFileString_sbs4_0p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_inel_2Dhistos.root";
double dx_pn_sbs4_0p = 0;
double sigma_p_sbs4_0p =0.255;
double sigma_n_sbs4_0p = 0.252;
double sigma_dy_sbs4_0p =0.28; 


// SBS4 30p 
TString DataFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_2Dhistos.root";
TString ProtonFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deep_2Dhistos.root";
TString NeutronFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_deen_2Dhistos.root";
TString InelasticFileString_sbs4_30p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_30p_cuts_inel_2Dhistos.root";
double dx_pn_sbs4_30p = 0.738;
double sigma_p_sbs4_30p =0.255;
double sigma_n_sbs4_30p = 0.252;
double sigma_dy_sbs4_30p =0.28; 


// SBS4 50p 
TString DataFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_2Dhistos.root";
TString ProtonFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deep_2Dhistos.root";
TString NeutronFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_deen_2Dhistos.root";
TString InelasticFileString_sbs4_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs4_50p_cuts_inel_2Dhistos.root";
double dx_pn_sbs4_50p = 1.21;
double sigma_p_sbs4_50p =0.272;
double sigma_n_sbs4_50p = 0.273;
double sigma_dy_sbs4_50p =0.28; 


// SBS8 70p
TString DataFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_2Dhistos_sept11.root";
TString ProtonFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deep_2Dhistos.root";
TString NeutronFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_deen_2Dhistos.root";
TString InelasticFileString_sbs8_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_70p_cuts_inel_2Dhistos.root";
double dx_pn_sbs8_70p = 0.867;
double sigma_p_sbs8_70p =0.210;
double sigma_n_sbs8_70p = 0.202;
double sigma_dy_sbs8_70p =0.244;

// SBS8 50p
TString DataFileString_sbs8_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_50p_cuts_2Dhistos_sept24.root";
TString ProtonFileString_sbs8_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_50p_cuts_deep_2Dhistos_sept24.root";
TString NeutronFileString_sbs8_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_50p_cuts_deen_2Dhistos_sept24.root";
TString InelasticFileString_sbs8_50p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_50p_cuts_inel_2Dhistos.root";// not real
double dx_pn_sbs8_50p = 0.615;
double sigma_p_sbs8_50p =0.210;
double sigma_n_sbs8_50p = 0.200;
double sigma_dy_sbs8_50p =0.245;

// SBS8 100p
TString DataFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_2Dhistos_sept24.root";
TString ProtonFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_deep_2Dhistos_sept24.root";
TString NeutronFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_deen_2Dhistos_sept24.root";
TString InelasticFileString_sbs8_100p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs8_100p_cuts_inel_2Dhistos.root";// not real
double dx_pn_sbs8_100p = 1.23;
double sigma_p_sbs8_100p =0.210;
double sigma_n_sbs8_100p = 0.20;
double sigma_dy_sbs8_100p =0.245;


// SBS9 70p
TString DataFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_2Dhistos_sept24.root";
TString ProtonFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_deep_2Dhistos_sept24_z.root";
TString NeutronFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_deen_2Dhistos_sept24_z.root";
TString InelasticFileString_sbs9_70p = "/w/halla-scshelf2102/sbs/msatnik/GMn/output/sbs9_70p_cuts_inel_2Dhistos.root";
double dx_pn_sbs9_70p = 0.91;
double sigma_p_sbs9_70p =0.232;
double sigma_n_sbs9_70p = 0.227;
double sigma_dy_sbs9_70p =0.228;


//fuctions

void draw_hcal_fid(TString KineString="sbs4_30p"){// main
  gStyle->SetOptFit(11111);
  gStyle->SetCanvasPreferGL(1);
  gStyle -> SetOptStat(0);
  gStyle ->SetEndErrorSize(0);
  gStyle->SetNumberContours(255);

  

  Utility utilityHandler;

  if (KineString == "sbs4_30p"){
    DataFileString = DataFileString_sbs4_30p;
    ProtonFileString =ProtonFileString_sbs4_30p;
    NeutronFileString = NeutronFileString_sbs4_30p;
    InelasticFileString =  InelasticFileString_sbs4_30p;
    dx_pn = dx_pn_sbs4_30p;
    sigma_p = sigma_p_sbs4_30p;
    sigma_n = sigma_n_sbs4_30p;
    sigma_dy = sigma_dy_sbs4_30p; 
  } if (KineString == "sbs4_0p"){
    DataFileString = DataFileString_sbs4_0p;
    ProtonFileString =ProtonFileString_sbs4_0p;
    NeutronFileString = NeutronFileString_sbs4_0p;
    InelasticFileString =  InelasticFileString_sbs4_0p;
    dx_pn = dx_pn_sbs4_0p;
    sigma_p = sigma_p_sbs4_0p;
    sigma_n = sigma_n_sbs4_0p;
    sigma_dy = sigma_dy_sbs4_0p; 
  }else if (KineString == "sbs4_50p"){
    DataFileString = DataFileString_sbs4_50p;
    ProtonFileString =ProtonFileString_sbs4_50p;
    NeutronFileString = NeutronFileString_sbs4_50p;
    InelasticFileString =  InelasticFileString_sbs4_50p;
    dx_pn = dx_pn_sbs4_50p;
    sigma_p = sigma_p_sbs4_50p;
    sigma_n = sigma_n_sbs4_50p;
    sigma_dy = sigma_dy_sbs4_50p; 
  }else if (KineString == "sbs8_70p"){
    DataFileString = DataFileString_sbs8_70p;
    ProtonFileString =ProtonFileString_sbs8_70p;
    NeutronFileString = NeutronFileString_sbs8_70p;
    InelasticFileString =  InelasticFileString_sbs8_70p;
    dx_pn = dx_pn_sbs8_70p;
    sigma_p = sigma_p_sbs8_70p;
    sigma_n = sigma_n_sbs8_70p;
    sigma_dy = sigma_dy_sbs8_70p; 
  }else if (KineString == "sbs8_50p"){
    DataFileString = DataFileString_sbs8_50p;
    ProtonFileString =ProtonFileString_sbs8_50p;
    NeutronFileString = NeutronFileString_sbs8_50p;
    InelasticFileString =  InelasticFileString_sbs8_50p;
    dx_pn = dx_pn_sbs8_50p;
    sigma_p = sigma_p_sbs8_50p;
    sigma_n = sigma_n_sbs8_50p;
    sigma_dy = sigma_dy_sbs8_50p; 
  }else if (KineString == "sbs8_100p"){
    DataFileString = DataFileString_sbs8_100p;
    ProtonFileString =ProtonFileString_sbs8_100p;
    NeutronFileString = NeutronFileString_sbs8_100p;
    InelasticFileString =  InelasticFileString_sbs8_100p;
    dx_pn = dx_pn_sbs8_100p;
    sigma_p = sigma_p_sbs8_100p;
    sigma_n = sigma_n_sbs8_100p;
    sigma_dy = sigma_dy_sbs8_100p; 
  }else if (KineString == "sbs9_70p"){
    DataFileString = DataFileString_sbs9_70p;
    ProtonFileString =ProtonFileString_sbs9_70p;
    NeutronFileString = NeutronFileString_sbs9_70p;
    InelasticFileString =  InelasticFileString_sbs9_70p;
    dx_pn = dx_pn_sbs9_70p;
    sigma_p = sigma_p_sbs9_70p;
    sigma_n = sigma_n_sbs9_70p;
    sigma_dy = sigma_dy_sbs9_70p;  
  }else {
    std::cout<<"KineString: "<<KineString<<endl;
    std::cout<<"Error with kinematic setting"<<std::endl;
    return;
  }

  
  cout<<endl;
  cout<<"What we're running: "<<KineString<<endl;
  cout<<DataFileString<<endl;
  cout<<ProtonFileString<<endl;
  cout<<NeutronFileString<<endl;
  //cout<<InelasticFileString<<endl;
  cout<<"dx_pn = "<<dx_pn<<", sigma_p = "<<sigma_p<<", sigma_n = "<<sigma_n<<", sigma_dy = "<<sigma_dy<<endl;
  cout<<endl;

  // load files 
  TFile *f1 = TFile::Open(DataFileString); // data
  TFile *f2 = TFile::Open(ProtonFileString); //proton
  TFile *f3 = TFile::Open(NeutronFileString); //neutron
  //TFile *f4 = TFile::Open(InelasticFileString); // inelastic sim

  if (!f1 || !f2 || !f3) {
    std::cout<<"Error with loading one of the rootfiles"<<std::endl;
    return;
  }

  // if(!f4){
  //   std::cout<<"Error loading the inelastic root file"<<std::endl;
  //   return;
  // }

  // set location and name for output file 
  TString outputfilelocation="../../output/hcal_fid_plots/"+KineString;
  TString outputfilename = outputfilelocation +"/"+ KineString+".root";
  // Declare outfile
  TFile *fout = new TFile(outputfilename,"RECREATE");
  cout<<"writing to file: "<< outputfilename <<endl;

  // Load Histograms
  TH2D *hist_data = (TH2D*)f1->Get("hcal_x_exp__hcal_y_exp2");
  TH2D *hist_proton_sim = (TH2D*)f2->Get("hcal_x_exp__hcal_y_exp");
  TH2D *hist_neutron_sim = (TH2D*)f3->Get("hcal_x_exp__hcal_y_exp");
  //TH2D *hist_inel = (TH2D*)f4->Get("hcal_x_exp__hcal_y_exp");


  TH2D *hist_dxdy_data = (TH2D*)f1->Get("hcal_dx__hcal_dy");


  if (!hist_data||!hist_proton_sim || !hist_neutron_sim)
    {
      std::cout<<"One of the histograms is empty"<<std::endl;
      return;
    }

  // if (!hist_inel)
  //   {
  //     std::cout<<"Error with inelastic histogram"<<std::endl;
  //     return; 
  //   }
  
  hist_data->GetXaxis()->SetRangeUser(-1.25,1.25);
  hist_proton_sim->GetXaxis()->SetRangeUser(-1.25,1.25);
  hist_neutron_sim->GetXaxis()->SetRangeUser(-1.25,1.25);
  //hist_inel->GetXaxis()->SetRangeUser(-1.25,1.25);

  // make hcal position boundary vector of lines: Xi, Xf, Yi, Yf
  std::vector<double> hcal_pos_vector = utilityHandler.hcalpos();
  std::vector<TLine*> hcal_pos_box = utilityHandler.CreateBox(hcal_pos_vector , kGreen, 2);
  // make hcal active area vector of lines: Xi, Xf, Yi, Yf
  std::vector<double> hcal_active_area_vector = utilityHandler.hcalaa_mc();
  std::vector<TLine*> hcal_aa_box = utilityHandler.CreateBox( hcal_active_area_vector , kMagenta, 2);


  std::cout<<"HCal Active Area: "<<std::endl;
  for (auto coord : hcal_active_area_vector)
    {
      std::cout<<coord<<", ";
    }
  std::cout<<std::endl<<std::endl;

  // make hcal fiducal vector of lines: Xi, Xf, Yi, Yf
  std::vector<std::vector<double>> hcal_fiducal_vector_of_vectors;
  std::vector<std::vector<double>> hcal_fiducal_shift_vector_of_vectors;
  for (int i = 0; i<20; i++){
    double nsigma = i*0.2;
    std::vector<double> hcalfid = utilityHandler.hcalfid(sigma_p, sigma_n,sigma_dy,hcal_active_area_vector,nsigma,nsigma,nsigma);
    hcal_fiducal_vector_of_vectors.push_back(hcalfid);
    std::vector<double>hcalfid_shift = hcalfid;
    hcalfid_shift[0] = hcalfid[0] + dx_pn;
    hcal_fiducal_shift_vector_of_vectors.push_back(hcalfid_shift);
    std::cout<<"nsigma: "<<nsigma<<std::endl;
    for(int j = 0; j<4; j++){
      std::cout<<hcalfid[j]<<", "<<std::endl;
    }
    std::cout<<"possible x-exp, y-exp cut"<<std::endl;;
    for(int j = 0; j<4; j++){
      std::cout<<hcalfid_shift[j]<<", "<<std::endl;
    }
    
    if(hcalfid_shift[0]>hcalfid_shift[1])
      {
	std::cout<<"x fiducail cut is too large"<<std::endl;
      }
    if(hcalfid[2]>hcalfid[3])
      {
	std::cout<<"y fiducail cut is too large"<<std::endl;
      }
    std::cout<<std::endl;
  }

  

  std::vector<double> hcal_fid_1=utilityHandler.hcalfid(sigma_p, sigma_n,sigma_dy,hcal_active_area_vector,1,1,1);
  std::vector<TLine*> box_fid1 = utilityHandler.CreateBox(hcal_fid_1, kRed, 2);

   
  std::vector<double> hcal_fid_0=utilityHandler.hcalfid(sigma_p, sigma_n,sigma_dy,hcal_active_area_vector,0,0,0);
  std::vector<TLine*> box_fid0 = utilityHandler.CreateBox(hcal_fid_0, kRed, 2);

  // shifting Xi  by the proton-neutron seperation. (for neutron hyp)
  std::vector<double> hcal_fid_1_shift_xi = hcal_fid_1;
  hcal_fid_1_shift_xi[0] = hcal_fid_1[0]+dx_pn;
  std::vector<TLine*> hcal_shift_xi_box_1 = utilityHandler.CreateBox( hcal_fid_1_shift_xi , kOrange, 2);
  std::vector<double> hcal_fid_0_shift_xi = hcal_fid_0;
  hcal_fid_0_shift_xi[0] = hcal_fid_0[0]+dx_pn;
  std::vector<TLine*> hcal_shift_xi_box_0 = utilityHandler.CreateBox( hcal_fid_0_shift_xi , kOrange, 2);
        
  //shifting Xf by the proton-neutron seperation. (for proton hyp) 
  std::vector<double> hcal_fid_1_shift_xf = hcal_fid_1;
  hcal_fid_1_shift_xf[1] = hcal_fid_1[1]-dx_pn;
  std::vector<TLine*> hcal_shift_xf_box_1 = utilityHandler.CreateBox( hcal_fid_1_shift_xf , kOrange, 2);
  std::vector<double> hcal_fid_0_shift_xf = hcal_fid_0;
  hcal_fid_0_shift_xf[1] = hcal_fid_0[1]-dx_pn;
  std::vector<TLine*> hcal_shift_xf_box_0 = utilityHandler.CreateBox( hcal_fid_0_shift_xf , kOrange, 2);
  
  
  TH2* hist_data_proton_hyp = utilityHandler.Shift2DHistogramY(hist_data, -dx_pn);
  TH2* hist_proton_sim_proton_hyp = utilityHandler.Shift2DHistogramY(hist_proton_sim, -dx_pn);
  TH2* hist_neutron_sim_proton_hyp = utilityHandler.Shift2DHistogramY(hist_neutron_sim, -dx_pn);


  
 
  

  TCanvas *data_canvas =  new TCanvas("data_canvas", "data_canvas", 800, 600);
  data_canvas ->Divide(2,1);
  data_canvas->cd(1);
  utilityHandler.adjustPadVirtual(gPad,0.10, 0.15, 0.10, 0.10);
  // neutron hyp
  hist_data->Draw("colz");
  data_canvas->Update();
  // draw fid box with Xi shift for pn seperation
  for (auto* line : hcal_shift_xi_box_1) {
    line->Draw("same");
    data_canvas->Update();
  }
  // draw fid box
  for (auto* line : box_fid1) {
    line->Draw("same");
    data_canvas->Update();
  } 
  // draw pos box
  for (auto* line : hcal_pos_box) {
    line->Draw("same");
    data_canvas->Update();
  }
  // draw pos box
  for (auto* line : hcal_aa_box) {
    line->Draw("same");
    data_canvas->Update();
  }
  // proton hyp
  data_canvas->cd(2);
  utilityHandler.adjustPadVirtual(gPad, 0.10, 0.15, 0.10, 0.10);
  hist_data_proton_hyp->Draw("colz");
  data_canvas->Update();
  // draw fid box with Xf shift for pn seperation
  for (auto* line : hcal_shift_xf_box_1) {
    line->Draw("same");
    data_canvas->Update();
  }
  // draw fid box
  for (auto* line : box_fid1) {
    line->Draw("same");
    data_canvas->Update();
  } 
  // draw pos box
  for (auto* line : hcal_pos_box) {
    line->Draw("same");
    data_canvas->Update();
  }
  // draw pos box
  for (auto* line : hcal_aa_box) {
    line->Draw("same");
    data_canvas->Update();
  }


  TCanvas *sim_canvas =  new TCanvas("sim_canvas", "sim_canvas", 800, 600);
  utilityHandler.adjustCanvas(sim_canvas);
  sim_canvas ->Divide(2,1);
  sim_canvas->cd(1);
  utilityHandler.adjustPadVirtual(gPad,0.10, 0.15, 0.10, 0.10);
  // neutron hyp
  hist_proton_sim->Draw("colz");
  sim_canvas->Update();
  // draw fid box with Xi shift for pn seperation
  for (auto* line : hcal_shift_xi_box_1) {
    line->Draw("same");
    sim_canvas->Update();
  }
  // draw fid box
  for (auto* line : box_fid1) {
    line->Draw("same");
    sim_canvas->Update();
  } 
  // draw pos box
  for (auto* line : hcal_pos_box) {
    line->Draw("same");
    sim_canvas->Update();
  }
  // draw pos box
  for (auto* line : hcal_aa_box) {
    line->Draw("same");
    sim_canvas->Update();
  }
  // proton hyp
  sim_canvas->cd(2);
  utilityHandler.adjustPadVirtual(gPad,0.10, 0.15, 0.10, 0.10);
  hist_proton_sim_proton_hyp->Draw("colz");
  sim_canvas->Update();
  // draw fid box with Xf shift for pn seperation
  for (auto* line : hcal_shift_xf_box_1) {
    line->Draw("same");
    sim_canvas->Update();
  }
  // draw fid box
  for (auto* line : box_fid1) {
    line->Draw("same");
    sim_canvas->Update();
  } 
  // draw pos box
  for (auto* line : hcal_pos_box) {
    line->Draw("same");
    sim_canvas->Update();
  }
  // draw pos box
  for (auto* line : hcal_aa_box) {
    line->Draw("same");
    sim_canvas->Update();
  }

    
    
  TCanvas *proton_hyp_canvas = new TCanvas("proton_hyp_canvas", "proton_hyp_canvas", 800, 600);
  hist_data_proton_hyp->Draw("colz");
  proton_hyp_canvas->Update();
  // draw fid box with Xf shift for pn seperation
  for (auto* line : hcal_shift_xf_box_1) {
    line->Draw("same");
    proton_hyp_canvas->Update();
  }
  // draw fid box
  for (auto* line : box_fid1) {
    line->Draw("same");
    proton_hyp_canvas->Update();
  } 
  // draw pos box
  for (auto* line : hcal_pos_box) {
    line->Draw("same");
    proton_hyp_canvas->Update();
  }
  // draw pos box
  for (auto* line : hcal_aa_box) {
    line->Draw("same");
    proton_hyp_canvas->Update();
  }

  TH2D* cut_test = utilityHandler.CutXYRangeTH2D(hist_data, hcal_fid_1_shift_xi[2], hcal_fid_1_shift_xi[3], hcal_fid_1_shift_xi[0],hcal_fid_1_shift_xi[1]);
 
  
  TCanvas *neutron_hyp_canvas = new TCanvas("neutron_hyp_canvas", "neutron_hyp_canvas", 800, 600);
  //hist_data->Draw("colz");
  cut_test->Draw("colz");
  neutron_hyp_canvas->Update();
  // draw fid box with Xi shift for pn seperation
  for (auto* line : hcal_shift_xi_box_1) {
    line->Draw("same");
    neutron_hyp_canvas->Update();
  }
  // draw fid box
  for (auto* line : box_fid1) {
    line->Draw("same");
    neutron_hyp_canvas->Update();
  } 
  // draw pos box
  for (auto* line : hcal_pos_box) {
    line->Draw("same");
    neutron_hyp_canvas->Update();
  }
  // draw pos box
  for (auto* line : hcal_aa_box) {
    line->Draw("same");
    neutron_hyp_canvas->Update();
  }


  TCanvas *spotcut_canvas = new TCanvas("spotcut_canvas", "spotcut_canvas", 800, 600);
  
  TEllipse* proton_spot = utilityHandler.CreateEllipse(0, -dx_pn, sigma_dy, sigma_p, 0, kRed, 2);

   TEllipse* neutron_spot = utilityHandler.CreateEllipse(0, 0, sigma_dy, sigma_n, 0, kRed, 2);
  hist_dxdy_data->Draw("colz");
  proton_spot->Draw("same");
  neutron_spot->Draw("same");
  spotcut_canvas->Update();

  fout->Write();

  //f1->Close();
  // f2->Close();
  //f3->Close();
  
}//end main
