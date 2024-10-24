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

// Global Params
// physical constants 
const double M_p = 0.93827208816; // Proton mass in GeV
const double M_n = 0.93956542052; // Neutron mass in GeV
const double Pi = acos(-1);
const double alpha = 1.0 / 137.0; // Fine-structure constant
const double GD_const = 0.71; // Dipole form factor range factor
const double mu_p = 2.79284735; //Proton magnetic moment
const double mu_n = -1.9130427; //Neutron anomalous magnetic moment
const double M_pi = 0.13957; // Charged pion mass in GeV. Need for Ye fit. 

const double ye_tcut = 4* M_pi * M_pi;
const double ye_tnot = -0.7; // Ye parameterization constant t0
const double GE_over_GD_errslope = 0.0184; //via interpolated fit to Arrington07 list
const double GM_over_GDmup_errslope = 0.0025; //via interpolated fit to Arrington07 list

// Function Headers 
double calculate_Q2(double E, double theta, double M);
double calculate_tau(double Q2, double M);
double calculate_epsilon(double tau, double theta);

// Fits to world data
double GEp_Kelly(double tau);
double GMp_Kelly(double tau);
double GMn_Kelly(double tau);
double calculate_Kelly(double tau, const std::vector<double>& a_coef, const std::vector<double>& b_coef);
double GEn_Riordan(double tau, double GD);
void GEp_Arrington07(double &GEp, double &GEp_err, double tau, double Q2, double GD);
void GMp_Arrington07(double &GMp, double &GMp_err, double tau, double Q2, double GD);
void GEp_Arrington07_borntpe(double &GEp, double & GEp_err, double tau, double Q2, double GD);
void GMp_Arrington07_borntpe(double &GMp, double &GMp_err, double tau, double Q2, double GD);

double calculateYeZ(double Q2);
double calculateYeX(double Q2);
double calculate_Ye(const std::vector<double>& a_coef, double yez);
double calculate_Ye_err(const std::vector<double>& b_coef, double yex, double GD);

double GEn_Ye(double yez);
double GEn_Ye_err(double yex, double GD);


double calculate_reduced_cross_section(double GM, double GE , double tau, double epsilon);
double calculate_reduced_cross_section_err(double GM, double GM_err, double GE, double GE_err , double tau, double epsilon);
double calculate_GMn(double Rsf, double RCSR_simc, double GEp, double GMp, double GEn, double epsilon_p, double epsilon_n, double tau_p, double tau_n);
double calculate_GMn_err(double Rsf, double Rsf_err, double RCSR_simc, double GMn, double GEp, double GEp_err, double GMp, double GMp_err, double GEn, double GEn_err, double epsilon_p, double epsilon_n, double tau_p, double tau_n);


std::vector<double> GEp_Ye_coef = {
  0.239163298067,
  -1.109858574410,
  1.444380813060,
  0.479569465603,
  -2.286894741870,
  1.126632984980,
  1.250619843540,
  -3.631020471590,
  4.082217023790,
  0.504097346499,
  -5.085120460510,
  3.967742543950,
  -0.981529071103};

std::vector<double> GMp_over_mup_Ye_coef = {
  0.264142994136,
  -1.095306122120,
  1.218553781780,
  0.661136493537,
  -1.405678925030,
  -1.356418438880,
  1.447029155340,
  4.235669735900,
  -5.334045653410,
  -2.916300520960,
  8.707403067570,
  -5.706999943750,
  1.280814375890};

std::vector<double> GEn_Ye_coef = {
  0.048919981379,
  -0.064525053912,
  -0.240825897382,
  0.392108744873,
  0.300445258602,
  -0.661888687179,
  -0.175639769687,
  0.624691724461,
  -0.077684299367,
  -0.236003975259,
  0.090401973470};


std::vector<double> GMn_over_mun_Ye_coef = {
  0.257758326959,
  -1.079540642058,
  1.182183812195,
  0.711015085833,
  -1.348080936796,
  -1.662444025208,
  2.624354426029,
  1.751234494568,
  -4.922300878888,
  3.197892727312,
  -0.712072389946};
  

  
std::vector<double> GEp_Ye_err_coef = {
  -1.97750308,
  -4.46566998*pow(10,-1),
  2.94508717*pow(10,-1),
  1.54467525,
  9.05268347*pow(10,-1),
  -6.00008111*pow(10,-1),
  -1.10732394,
  -9.85982716*pow(10,-2),
  4.63035988*pow(10,-1),
  1.37729116*pow(10,-1),
  -7.82991627*pow(10,-2),
  -3.63056932*pow(10,-2),
  2.64219326*pow(10,-3),
  3.13261383*pow(10,-3),
  3.89593858*pow(10,-4)};

std::vector<double> GMp_over_mup_Ye_err_coef = {
  -1.76549673,
  1.67218457*pow(10,-1),
  -1.20542733,
  -4.72244127*pow(10,-1),
  1.41548871,
  6.61320779*pow(10,-1),
  -8.16422909*pow(10,-1),
  -3.73804477*pow(10,-1),
  2.62223992*pow(10,-1),
  1.28886639*pow(10,-1),
  -3.90901510*pow(10,-2),
  -2.44995181*pow(10,-2),
  8.34270064*pow(10,-4),
  1.88226433*pow(10,-3),
  2.43073327*pow(10,-4)};

std::vector<double> GEn_Ye_err_coef = {
  -2.07194073,
  1.13809127,
  1.01431277,
  -3.13301380*pow(10,-1),
  -2.73293676*pow(10,-1),
  2.57350595*pow(10,-1),
  -2.06042113*pow(10,-1),
  -1.68497332*pow(10,-1),
  1.37784515*pow(10,-1),
  7.57591964*pow(10,-2),
  -2.67511301*pow(10,-2),
  -1.72573088*pow(10,-2),
  7.03581500*pow(10,-4),
  1.47962095*pow(10,-3),
  1.97375221*pow(10,-4)};
  
  
std::vector<double>  GMn_over_mun_Ye_err_coef= {
  -2.06920873,
  6.431564*pow(10,-2),
  -3.55593786*pow(10,-1),
  4.1489766*pow(10,-1),
  1.95746824,
  2.705257*pow(10,-1),
  -1.52685784,
  -4.43527359*pow(10,-1),
  5.16884065*pow(10,-1),
  2.07915837*pow(10,-1),
  -7.48665703*pow(10,-2),
  -4.25411431 *pow(10,-2),
  1.54965016 *pow(10,-3),
  3.25322279*pow(10,-3),
  4.20819518*pow(10,-4)};
  
  
// defaults to be overwritten
double E_beam = 4.0268; // Beam energy in GeV/c2 
double theta_BB_deg = 49.0; // Central anlge of BigBite
double theta_BB = theta_BB_deg * TMath::DegToRad();// convert to radians
double Rsf =  1.06448; // maybe real
double Rsf_err = 0.106448; // arbitrary 10% error


// SBS 4 30%
double E_beam_sbs4_30p = 3.7393; // Beam energy in GeV/c2  3.728
double theta_BB_deg_sbs4_30p = 36.0; // Central anlge of BigBite 
double Rsf_sbs4_30p=  0.9655;
double Rsf_err_sbs4_30p =  0.09655;// arbitray 10% error 
//double Rsf_err_sbs4_30p = 0.006;

// // Hard-code sbs4 50% for the moment
double E_beam_sbs4_50p = 3.7393; // Beam energy in GeV/c2 
double theta_BB_deg_sbs4_50p = 36.0; // Central anlge of BigBite 
double Rsf_sbs4_50p = 1.0485;
double Rsf_err_sbs4_50p = 0.019;

/// SBS 9 70%
double E_beam_sbs9_70p = 4.0268; // Beam energy in GeV/c2 
double theta_BB_deg_sbs9_70p = 49.0; // Central anlge of BigBite 
double Rsf_sbs9_70p =  1.06448; // maybe real
double Rsf_err_sbs9_70p = 0.106448; // arbitrary 10% error


// SBS 8 70%
double E_beam_sbs8_70p = 5.9826; // Beam energy in GeV/c2 
double theta_BB_deg_sbs8_70p = 26.5; // Central anlge of BigBite 
double Rsf_sbs8_70p =  1.089; //
double Rsf_err_sbs8_70p = 0.1089;// Arbitrary 10% error 
//double Rsf_err = 0.01; // 


void getGMn(TString KineString="sbs4_30p"){// MAIN

  // // Hard-code sbs4 30% for the moment
  // cout<<"Running SBS4 30%"<<endl;
  // double E_beam = 3.7393; // Beam energy in GeV/c2  3.728
  // double theta_BB_deg = 36.0; // Central anlge of BigBite 
  // double theta_BB = theta_BB_deg * TMath::DegToRad();// convert to radians
  // double Rsf =  0.9655;
  // double Rsf_err =  0.09655;// arbitray 10% error 
  // //double Rsf_err = 0.01;


  // // Hard-code sbs4 50% for the moment
  // cout<<"Running SBS4 50%"<<endl;
  // double E_beam = 3.7393; // Beam energy in GeV/c2 
  // double theta_BB_deg = 36.0; // Central anlge of BigBite 
  // double theta_BB = theta_BB_deg * TMath::DegToRad();// convert to radians
  // double Rsf = 1.0485;
  // double Rsf_err = 0.019;
  


  // // Hard-code sbs8 70% for the moment
  // cout<<"Running SBS8 70%"<<endl;
  // double E_beam = 5.9826; // Beam energy in GeV/c2 
  // double theta_BB_deg = 26.5; // Central anlge of BigBite 
  // double theta_BB = theta_BB_deg * TMath::DegToRad();// convert to radians
  // double Rsf =  1.089; //
  // double Rsf_err = 0.1089;// Arbitrary 10% error 
  // //double Rsf_err = 0.01; // 
 
  
  // // Hard-code sbs9 70% for the moment
  // cout<<"Running SBS9 70%"<<endl;
  // double E_beam = 4.0268; // Beam energy in GeV/c2 
  // double theta_BB_deg = 49.0; // Central anlge of BigBite 
  // double theta_BB = theta_BB_deg * TMath::DegToRad();// convert to radians
  // double Rsf =  1.06448; // maybe real
  // double Rsf_err = 0.106448; // arbitrary 10% error
  // //double Rsf_err = 0.004; // maybe real

  
  if (KineString == "sbs4_30p"){
    E_beam = E_beam_sbs4_30p;
    theta_BB_deg = theta_BB_deg_sbs4_30p;
    Rsf = Rsf_sbs4_30p;
    Rsf_err = Rsf_err_sbs4_30p;
  }else if (KineString == "sbs4_50p"){
    E_beam = E_beam_sbs4_50p;
    theta_BB_deg = theta_BB_deg_sbs4_50p;
    Rsf = Rsf_sbs4_50p;
    Rsf_err = Rsf_err_sbs4_50p;
  }else if (KineString == "sbs8_70p"){
    E_beam = E_beam_sbs8_70p;
    theta_BB_deg = theta_BB_deg_sbs8_70p;
    Rsf = Rsf_sbs8_70p;
    Rsf_err = Rsf_err_sbs8_70p;
  }
  // else if (KineString == "sbs8_50p"){ 
  // }else if (KineString == "sbs8_100p"){
  // }
  else if (KineString == "sbs9_70p"){
    E_beam = E_beam_sbs9_70p;
    theta_BB_deg = theta_BB_deg_sbs9_70p;
    Rsf = Rsf_sbs9_70p;
    Rsf_err = Rsf_err_sbs9_70p;
  }else {
    std::cout<<"Error with kinematic setting"<<std::endl;
    return;
  }


  double theta_BB = theta_BB_deg * TMath::DegToRad();// convert to radians

  cout<<endl;
  cout<<"What we're running: "<<KineString<<endl;

  cout << "E_beam = "<< E_beam <<endl;
  cout<< "theta BB = "<<theta_BB_deg<<" deg, or "<<theta_BB<<" in radians "<<endl;
  cout<<"user input Rsf = "<<Rsf<<" +/- "<<Rsf_err<<endl;
  cout<<endl;
  
  double M_avg =(M_n + M_p) / 2; // Average mass of the proton and neutron
  double Q2_cent  =  calculate_Q2(E_beam, theta_BB, M_avg); // calculate the central Q2 usign the average mass

  //cout<<"M_n = "<<M_n<<", M_p = "<<M_p<<" , M_avg = " << M_avg<<endl;
  cout<< "Q2 = "<< Q2_cent<<endl;
 
  //Q2_cent = 3.0326; //GET RID OF THIS TEST
  
  double GD = pow(1+Q2_cent/0.71,-2); // dipole form factor

  double tau_p = calculate_tau(Q2_cent, M_p);
  double tau_n = calculate_tau(Q2_cent, M_n);

  double epsilon_p = calculate_epsilon(tau_p, theta_BB);
  double epsilon_n = calculate_epsilon(tau_n, theta_BB);

  cout<<"GD = "<<GD<<endl;
  cout<< "tau_n = "<< tau_n <<", tau_p = "<<tau_p<<endl;
  cout<< "epsilon_n = "<< epsilon_n <<", epsilon_p = "<<epsilon_p<<endl;
  cout<<endl;

  // Calculated the form factors used by simc
  /// simc uses tau_p in the code, and Kelly uses tau_p in the paper
  double GEp_simc = GEp_Kelly(tau_p);
  double GMp_simc = GMp_Kelly(tau_p);
  double GEn_simc = GEn_Riordan(tau_p,GD);
  /// Oh sebastian has a  1.05 times the whole thing when looking at the world data... why?
  // These are the parameters Sebastian has in his code for GEn Riordan. I should try and verify them.
  double GMn_simc = GMn_Kelly(tau_p);

  cout<<"Form Factors used in simc"<<endl;
  cout<<"JJ Kelly: "<<endl;
  cout<<"    GEp_simc = "<<GEp_simc<<endl;
  cout<<"    GEp_simc/GD = "<<GEp_simc/GD<<endl; 
  cout<<"    GMp_simc  = "<< GMp_simc <<endl;
  cout<<"    GMp_simc/(GD*mu_p) = "<<GMp_simc/(GD*mu_p)<<endl;
  cout<<"    GMn_simc = "<<GMn_simc<<endl;
  cout<<"    GMn_simc/(GD*mu_n) = "<<GMn_simc/(GD*mu_n)<<endl;
  cout<<"    mu_p*GEp/GMp = "<< mu_p*GEp_simc/GMp_simc<<endl;
  cout<<"Riordan:"<<endl;
  cout<<"    GEn_simc = "<<GEn_simc<<endl;
  cout<<"    GEn_simc/GD = "<<GEn_simc/GD<<endl;
  cout<<endl;

  double reduced_cross_section_n = calculate_reduced_cross_section(GMn_simc, GEn_simc, tau_n, epsilon_n);
  double reduced_cross_section_p = calculate_reduced_cross_section(GMp_simc,GEp_simc,tau_p,epsilon_p);

  //// Calculate the reduced cross section ratio found in simc
  double RCSR_simc = reduced_cross_section_n / reduced_cross_section_p;
  //// Get the corrected reduced cross section ratio by mulitplying it by Rsf
  double corrected_RCSR_simc = Rsf * RCSR_simc;
  
  cout<<"Calculating the reduced cross sections used in simc"<<endl;
  cout<< "Neutron reduced cross section = "<< reduced_cross_section_n<<endl;
  cout<<"Proton reduced cross section = "<<reduced_cross_section_p<<endl;
  cout<<"Reduced cross section ratio: n/p = "<< RCSR_simc<<endl;
  cout<< "Corrected simc Reduced Cross Section Ratio: "<<endl;
  cout<<"Rsf * (simc reduced cross section ratio) = "<<Rsf << " * "<< RCSR_simc <<" = "<<corrected_RCSR_simc<<endl;
  cout<<endl;


  //// Calculate form factor values from the most modern fits to world data.
  //// Consider this model dependent from here out.
  double GEp, GEp_err;
  double GMp, GMp_err;
  double GEn, GEn_err;
  
  double GEp_borntpe, GEp_borntpe_err;
  double GMp_borntpe, GMp_borntpe_err;


  /// Calculate GEn from Ye fit
  cout<<"Ye Calculation"<<endl;
    
  double Ye_z = calculateYeZ(Q2_cent);
  double Ye_x  = calculateYeX(Q2_cent);

  cout<< "Ye_z = "<< Ye_z<< ", Ye_x =  "<<Ye_x<<endl;

  GEn = calculate_Ye(GEn_Ye_coef ,Ye_z);
  GEn_err = calculate_Ye_err(GEn_Ye_err_coef, Ye_x,GD);
  cout <<"GEn  = "<<GEn <<" +/- "<<GEn_err<<endl;
  cout <<"GEn/GD  = "<<GEn/GD <<" +/- "<<GEn_err/GD<<endl;

  int using_ye = 1;

  if (using_ye){
  
  
    double GMp_over_mup_ye = calculate_Ye(GMp_over_mup_Ye_coef ,Ye_z);
    GMp = mu_p * GMp_over_mup_ye;
    double GMp_over_mup_ye_err = calculate_Ye_err(GMp_over_mup_Ye_err_coef, Ye_x,GD);
    GMp_err =mu_p*GMp_over_mup_ye_err;
    GEp =  calculate_Ye(GEp_Ye_coef ,Ye_z);
    GEp_err = calculate_Ye_err(GEp_Ye_err_coef, Ye_x,GD);
    
    cout<<"GMp = "<<GMp<<" +/- "<<GMp_err<<endl;
    cout <<"GMp/(mu_p*GD)  = "<<GMp_over_mup_ye/GD <<" +/- "<<GMp_over_mup_ye_err/GD<<endl;
    cout<<"GEp = "<<GEp <<" +/= "<<GEp_err<<endl;
    cout<<"GEp/GD = "<<GEp/GD <<" +/= "<<GEp_err/GD<<endl;
  
    
    cout<<endl;
  
  } else {
  
    GEp_Arrington07(GEp, GEp_err, tau_p, Q2_cent, GD); // reference params for GEp, GEp_err
    GMp_Arrington07(GMp, GMp_err, tau_p, Q2_cent, GD); // reference params for GMp, GMp_err
    GEp_Arrington07_borntpe(GEp_borntpe, GEp_borntpe_err, tau_p, Q2_cent, GD); // reference params for GEp_borntpe, GEp_borntpe_err
    GMp_Arrington07_borntpe(GMp_borntpe, GMp_borntpe_err, tau_p, Q2_cent,GD); // reference params for GMp_borntpe, GMp_borntpe_err


   

    cout<<"Arrington07 Calculations"<<endl;
    cout<<"With born tpe corrections"<<endl;
    cout<<"    GEp  = "<< GEp<<" +/- "<<GEp_err<<endl;
    cout<<"    GEp/GD = "<<GEp/GD<<" +/- "<<GEp_err/GD<<endl;
    cout<<"    GMp  = "<< GMp<< " +/- "<<GMp_err<<endl;
    cout<<"    GMp/(mu_p*GD) = "<<GMp/(mu_p*GD)<< " +/- " << GMp_err/(mu_p*GD)<<endl;
    cout<<"    mu_p*GEp/GMp = "<< mu_p*GEp/GMp<<endl;
    cout<<"without born tpe corrections. Needs work. Confused. "<<endl;
    cout<<"    GEp born tpe = "<< GEp_borntpe<<" +/- "<<GEp_borntpe_err<<endl;
    cout<<"    GEp/GD = "<<GEp_borntpe/GD<<" +/- "<<GEp_borntpe_err/GD<<endl;
    cout<<"    GMp born tpe = "<< GMp_borntpe<<" +/- "<< GMp_borntpe_err<<endl;
    cout<<"    GMp/(mu_p*GD) = "<< GMp_borntpe/(mu_p*GD)<<" +/- "<< GMp_borntpe_err/(mu_p*GD)<<endl;

    // not sure how to use the born tpe values or exactly where they come from. 


  }

  {
  
    GEp_Arrington07(GEp, GEp_err, tau_p, Q2_cent, GD); // reference params for GEp, GEp_err
    GMp_Arrington07(GMp, GMp_err, tau_p, Q2_cent, GD); // reference params for GMp, GMp_err
    GEp_Arrington07_borntpe(GEp_borntpe, GEp_borntpe_err, tau_p, Q2_cent, GD); // reference params for GEp_borntpe, GEp_borntpe_err
    GMp_Arrington07_borntpe(GMp_borntpe, GMp_borntpe_err, tau_p, Q2_cent,GD); // reference params for GMp_borntpe, GMp_borntpe_err


   

    // cout<<"Arrington07 Calculations"<<endl;
    // cout<<"With born tpe corrections"<<endl;
    // cout<<"    GEp  = "<< GEp<<" +/- "<<GEp_err<<endl;
    // cout<<"    GEp/GD = "<<GEp/GD<<" +/- "<<GEp_err/GD<<endl;
    // cout<<"    GMp  = "<< GMp<< " +/- "<<GMp_err<<endl;
    // cout<<"    GMp/(mu_p*GD) = "<<GMp/(mu_p*GD)<< " +/- " << GMp_err/(mu_p*GD)<<endl;
    // cout<<"    mu_p*GEp/GMp = "<< mu_p*GEp/GMp<<endl;
    // cout<<"without born tpe corrections. Needs work. Confused. "<<endl;
    // cout<<"    GEp born tpe = "<< GEp_borntpe<<" +/- "<<GEp_borntpe_err<<endl;
    // cout<<"    GEp/GD = "<<GEp_borntpe/GD<<" +/- "<<GEp_borntpe_err/GD<<endl;
    // cout<<"    GMp born tpe = "<< GMp_borntpe<<" +/- "<< GMp_borntpe_err<<endl;
    // cout<<"    GMp/(mu_p*GD) = "<< GMp_borntpe/(mu_p*GD)<<" +/- "<< GMp_borntpe_err/(mu_p*GD)<<endl;

    // // not sure how to use the born tpe values or exactly where they come from. 


  }


  double GMn_over_mun_ye = calculate_Ye(GMn_over_mun_Ye_coef ,Ye_z);
  double GMn_ye = mu_n * GMn_over_mun_ye;
  double GMn_over_mun_ye_err = calculate_Ye_err(GMn_over_mun_Ye_err_coef, Ye_x,GD);
  double GMn_ye_err =mu_n*GMn_over_mun_ye_err;
  cout<<endl;
  cout<< "GMn from Ye Fit: "<<endl;
  cout<<"GMn = "<<GMn_ye<<" +/- "<<GMn_ye_err<<endl;
  cout <<"GMn/(mu_n*GD)  = "<<GMn_over_mun_ye/GD <<" +/- "<<GMn_over_mun_ye_err/GD<<endl;
  cout<<endl;
  
  

  
  //// Calculating the reduced cross section for protons from just world data.
  double RCSp = calculate_reduced_cross_section(GMp, GEp, tau_p, epsilon_p);
  double RCSp_err = calculate_reduced_cross_section_err(GMp, GMp_err, GEp, GEp_err, tau_p, epsilon_p);

  cout<<endl;
  cout<< "Proton Reduced Cross Section just from world data fit ="<<endl;
  cout<< RCSp <<" +/- "<<RCSp_err <<endl;
  cout<<endl;

 


  //// calculate GMn 
  double GMn = calculate_GMn(Rsf,RCSR_simc,GEp, GMp, GEn, epsilon_p, epsilon_n,tau_p,tau_n);
  double GMn_err = calculate_GMn_err(Rsf, Rsf_err, RCSR_simc, GMn, GEp, GEp_err, GMp, GMp_err, GEn, GEn_err,  epsilon_p, epsilon_n, tau_p, tau_n);

   
  double GMn_over_GDmun = abs( GMn / (GD*mu_n) );
  double GMn_over_GDmun_err = abs( GMn_err / (GD*mu_n));

  cout <<"---------------------------------------------------------------------------"<<endl;
  cout<< "Q2 = "<< Q2_cent<<endl;
  cout<<"Rsf = "<<Rsf<<" +/- "<<Rsf_err<<endl;
  cout<< "GMn = "<<GMn<<" +/- "<<GMn_err<<endl;
  cout<< "GMn/(mu_n*GD) = "<<GMn_over_GDmun <<" +/- "<<GMn_over_GDmun_err<<endl;
  cout<<"----------------------------------------------------------------------------"<<endl;
  cout<<endl;

  // cout<<"trying out the values without the born tpe correction"<<endl;
  // double GMn_borntpe = calculate_GMn(Rsf,RCSR_simc,GEp_borntpe, GMp_borntpe, GEn, epsilon_p, epsilon_n,tau_p,tau_n);
  //  double GMn_borntpe_err = calculate_GMn_err(Rsf, Rsf_err, RCSR_simc, GMn_borntpe, GEp_borntpe, GEp_borntpe_err, GMp_borntpe, GMp_borntpe_err, GEn, GEn_err,  epsilon_p, epsilon_n, tau_p, tau_n);   
  // double GMn_over_GDmun_borntpe = abs( GMn_borntpe / (GD*mu_n) );
  // double GMn_over_GDmun_borntpe_err = abs( GMn_borntpe_err / (GD*mu_n));
  // cout<< "GMn_borntpe = "<<GMn_borntpe<<" +/- "<<GMn_borntpe_err<<endl;
  // cout<< "GMn_borntpe/(mu_n*GD) = "<<GMn_over_GDmun_borntpe <<" +/- "<<GMn_over_GDmun_borntpe_err<<endl;
  // cout<<endl;
  
}// end MAIN







// Calculate Q2 for a given beam energy, BB angle, and nucleon mass 
// Solve for E' using q2 = -4 E E' sin^2 (theta/2) and q2 = -2 nu M and nu = E - E'
// Q2 = -q2
double calculate_Q2(double E, double theta, double M){
  double E_prime = E / (1+ (2*E / M) *pow(sin(theta/2),2) );
  cout <<"E_prime = "<< E_prime<<endl;
  double Q2 = 4 * E * E_prime * pow(sin(theta/2),2);
  return Q2; 
}

// calculate tau for a given Q2 and nucleon mass
double calculate_tau(double Q2, double M){
  return Q2 / (4 * M * M);
}

// calculate epsilon for a given tau and angle 
double calculate_epsilon(double tau, double theta){
  return 1 / ( 1+2*(1+tau) * pow( tan(theta/2),2)  );
}


//// Various functions to calculate Ye parameters
//// Calculate z value for the z expansion. 
double calculateYeZ(double Q2){
  double numerator = sqrt(ye_tcut+Q2) - sqrt(ye_tcut-ye_tnot);
  double denominator = sqrt(ye_tcut+Q2) + sqrt(ye_tcut-ye_tnot);
  return numerator / denominator;
}

/// calculate the ye x param needed for error calculation 
double calculateYeX(double Q2){
  double yex;
  if (Q2 == 1)
    yex = 0.00000001;
  else
    yex = log10(Q2);
  return yex;
}

/// Calculate GEn from z expansion using coefficients 
double GEn_Ye(double yez){
  double GEn = 0.048919981 +
    -0.064525054 * yez +
    -0.240825897 * pow(yez, 2) +
    0.392108745 * pow(yez, 3) +
    0.300445259 * pow(yez, 4) +
    -0.661888687 * pow(yez, 5) +
    -0.17563977 * pow(yez, 6) +
    0.624691724 * pow(yez, 7) +
    -0.077684299 * pow(yez, 8) +
    -0.236003975 * pow(yez, 9) +
    0.090401973 * pow(yez, 10);
  return GEn;
}



//// calculate the error on GEn from yex log10(Q2)
double GEn_Ye_err(double yex, double GD){
  double GEn_err_exp = -2.07194073 +
    1.13809127 * yex +
    1.01431277 * pow(yex, 2) +
    -0.31330138 * pow(yex, 3) +
    -0.273293676 * pow(yex, 4) +
    0.257350595 * pow(yex, 5) +
    -0.206042113 * pow(yex, 6) +
    -0.168497322 * pow(yex, 7) +
    0.137784515 * pow(yex, 8) +
    0.075759196 * pow(yex, 9) +
    -0.02675113 * pow(yex, 10) +
    -0.017525731 * pow(yex, 11) +
    0.000703582 * pow(yex, 12) +
    0.001479621 * pow(yex, 13) +
    0.000197375 * pow(yex, 14);

  double GEn_err = pow(10, GEn_err_exp) * GD;
  return GEn_err;
}


// Ripped directly from simulation
// g4sbs/src/G$SBSEventGen.cc line 790
// this is labled as "our fit" in the code, so I don't know where it comes from. 
double GEn_Riordan(double tau, double GD){
  double GEn =GD* (1.520*tau + 2.629*tau*tau + 3.055*tau*tau*tau)/(1.0+5.222*tau+0.040*tau*tau+11.438*tau*tau*tau);
  return GEn;
}

  

double calculate_Kelly(double tau, const std::vector<double>& a_coef, const std::vector<double>& b_coef){
  double numerator = 0;
  for (int i = 0; i<a_coef.size();i++)
    {
      numerator += a_coef[i]*pow(tau,i);
    }
  double denominator = 1;
  for (int i = 0; i<b_coef.size(); i++)
    {
      denominator += b_coef[i]*pow(tau,i);
    }
  return numerator/denominator;
}


// Calculate Ye form factor from the coefs and the ye z.
// G = GEp, GMp/mu_p, GEn, or GMn/mu_n
double calculate_Ye(const std::vector<double>& a_coef, double yez)
{
  double G = 0;
  for (int i = 0; i<a_coef.size();i++)
    {
      G += a_coef[i]*pow(yez, i);
    }
  return G;
}


// Calculate Ye Error  form factor from the coefs and the ye x.
// G_err = error(GEp), error( GMp/mu_p) , error( GEn) , or error(GMn/mu_n)
double calculate_Ye_err(const std::vector<double>& b_coef, double yex, double GD)
{
  double exponent = 0;
  for (int i = 0; i<b_coef.size();i++)
    {
      exponent += b_coef[i]*pow(yex, i);
    }
  double G_err = pow(10,exponent) * GD;
  return G_err;
}



// double GEn_Riordan_test(double tau, double GMn_Kelly, double GD)
// {
//   // if I'm reading this paper correctly, we modify the Kelly GMn value to get GEn.
//   double GEn = mu_n;
  
// }

// GEp calculatation from  JJ Kelly fit 
double GEp_Kelly(double tau){
  double GEp = (1.0 - 0.24 * tau) / (1.0 + 10.98 * tau + 12.82 * pow(tau,2) + 21.97 * pow(tau,3));
  return GEp;
}


// GMp calculation from JJ Kelly fit.
// The fit parameters are for GMp / mu_p
/// ripped directly from g4sbs/src/G4SBSEventGen.cc 
double GMp_Kelly(double tau){
  double GMp = 2.79 * (1 + 0.12 * tau) / (1.0 + 10.97 * tau + 18.86 * pow(tau,2)+ 6.55 * pow(tau,3) );
  /// in G4BS, they hard-coded mu_p as 2.79 in this calculation
  return GMp;
}

/// GMn calculation from JJ Kelly fit.
/// the fit parameters are for GMn/mu_n
/// ripped directly from g4sbs/src/G4SBSEventGen.cc 
double GMn_Kelly(double tau){
  double GMn = -1.913 * (1.0 + 2.33 * tau) / (1.0 + 14.72 * tau + 24.2 * pow(tau,2) + 84.1 * pow(tau,3));
  // in G4SBS, they hard-coded mu_n as -1.913 in this calculation
  return GMn;
}


// GEp calculation from Arrington07 fit. No borne tpe correction. 
void GEp_Arrington07(double &GEp, double &GEp_err, double tau, double Q2, double GD){
  double numerator = 1+ 3.439*tau - 1.602* pow(tau,2) + 0.068 * pow(tau,3);
  double denominator = 1 + 15.055 * tau + 48.061*pow(tau,2) + 99.304*pow(tau,3) + 0.012* pow(tau,4) + 8.65 * pow(tau,5);
  GEp = numerator/denominator; 
  GEp_err = GE_over_GD_errslope * Q2 * GD;
  return;
}


// GMp calculation from Arrington07 fit. Incluedes borntpe correction  
void GMp_Arrington07(double &GMp, double &GMp_err, double tau, double Q2, double GD){
  //// for GMp/mu_p
  double numerator = (1 - 1.465*tau + 1.26*pow(tau,2) + 0.262*pow(tau,3) );
  double denominator =  (1 + 9.627*tau + 0.*pow(tau,2) + 0.*pow(tau,3) + 11.179*pow(tau,4) + 13.245*pow(tau,5));
  double GMp_over_mup = numerator/denominator;
  GMp = mu_p * numerator/denominator;
  GMp_err = GM_over_GDmup_errslope * Q2 * GD *mu_p;
  return;
}

/// GEp calculation from Arrington07 fit. I belive this excludes born TPE corrections 
void GEp_Arrington07_borntpe(double &GEp, double & GEp_err, double tau, double Q2, double GD){
  double numerator = (1 - 1.651*tau + 1.287*pow(tau,2) - 0.185*pow(tau,3));
  double denominator =  (1 + 9.531*tau + 0.591*pow(tau,2) + 0.0*pow(tau,3) + 0.0*pow(tau,4) + 4.994*pow(tau,5));
  GEp = numerator/denominator;
  GEp_err = GE_over_GD_errslope * Q2 * GD;
  return;
}

//// GMp calculation from Arrington07 fit. Excludes the TPE correction 
void GMp_Arrington07_borntpe(double &GMp, double &GMp_err, double tau, double Q2, double GD){
  double numerator = (1 - 2.151*tau + 4.261*pow(tau,2) + 0.159*pow(tau,3));
  double denominator = (1 + 8.647*tau + 0.001*pow(2,tau) + 5.245*pow(3,tau) + 82.817*pow(4,tau) + 14.191*pow(5,tau));
  double GMp_over_mup = numerator/denominator;
  GMp = mu_p * numerator/denominator;
  GMp_err = GM_over_GDmup_errslope * Q2 * GD * mu_p;
  return; 
}


/// Calculate the reduced cross section 
double calculate_reduced_cross_section(double GM, double GE , double tau, double epsilon){
  return GM*GM + (epsilon/tau)*GE*GE;
}

// Calculate the reduced cross section err
double calculate_reduced_cross_section_err(double GM, double GM_err, double GE, double GE_err , double tau, double epsilon){
  double GM_partial = 2 * GM;
  double GE_partial  = 2 * (epsilon/tau)* GE;

  double variance = pow(GM_partial,2)* pow(GM_err,2) + pow(GE_partial,2)*pow(GE_err,2);
  double stdev  =  sqrt(variance);
  return stdev; 
}



// Calculate GMn
double calculate_GMn(double Rsf, double RCSR_simc, double GEp, double GMp, double GEn, double epsilon_p, double epsilon_n, double tau_p, double tau_n){
  double GMn_squared = Rsf * RCSR_simc*(GMp*GMp + epsilon_p / tau_p * GEp * GEp) - epsilon_n / tau_n * GEn * GEn;
  return sqrt(GMn_squared);
}

//// Calculate GMn error
//// Assuming uncorrelated error
double calculate_GMn_err(double Rsf, double Rsf_err, double RCSR_simc, double GMn, double GEp, double GEp_err, double GMp, double GMp_err, double GEn, double GEn_err, double epsilon_p, double epsilon_n, double tau_p, double tau_n){

  double Rsf_partial  = RCSR_simc * (GMp * GMp +(epsilon_n / tau_n) * GEp * GEp)   / (GMn*2);
  double GMp_partial = 2 * Rsf * RCSR_simc * GMp / (GMn*2);
  double GEp_partial  = 2 * Rsf * RCSR_simc * (epsilon_p / tau_p) * GEp / (GMn*2);
  double GEn_partial  = 2 * (epsilon_n/tau_n)*GEn / (GMn*2);

  double variance  = pow(Rsf_partial,2)*pow(Rsf_err,2) + pow(GMp_partial,2) * pow(GMp_err,2) +pow(GEp_partial,2) *pow(GEp_err,2) +  pow(GEn_partial,2)*pow(GEn_err,2) ;
  double stdev = sqrt(variance);
  return stdev;
}
