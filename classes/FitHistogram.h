#ifndef FITHISTOGRAM_H
#define FITHISTOGRAM_H

#include <TH1.h>
#include <TF1.h>

class FitHistogram {
 private:

 public:
  TH1D *hist_p;  // proton simulation histogram
  TH1D *hist_n; // neutron simulation histogram
  TF1 *fitFunc; // Custom fit function
  double xMin; // min x on fit range
  double xMax; // max x on fit range
  double scale_p; // proton scale factor from fit
  double scale_p_err; //
  double scale_n; // neutron scale factor from fit
  double scale_n_err;
  double shift_p;// how much the neutron histogram should be shifted from fit
  double shift_p_err;
  double shift_n;// how much the neutron histogram should be shifted from fit
  double shift_n_err;
  double R; // ratio n/p from fit
  double R_err; //
  double BGscale; // scale for the provided bg polynomial used in mc_p_n_BG_R
  double BGscale_err;//
  double polyorder; // order for the polynomial background used in mc_p_n_poly_R
  double ChiSq; // chi^2 from fit
  double NDF; // number of degrees of freedom from fit.
  std::vector<double> poly_result; // background from poly fit (used in mc_p_n_poly2_R and mc_p_n_poly_R)
  std::vector<double> poly_result_err;
  std::vector<double> polyCoefficients;  // Polynomial coefficients (passed to mc_p_n_BG_R from user)
    
  FitHistogram(TH1D *hist_p, TH1D *hist_n, double xmin, double xmax);  // Constructor
  ~FitHistogram();  // Destructor

  
  Double_t mc_p_n_poly2_R(Double_t *x, Double_t *par);  // Fit function for data-mc comparison
  Double_t mc_p_n_poly_R(Double_t *x, Double_t *par);  // Fit function for data-mc comparison
  Double_t mc_p_n_BG_R(Double_t *x, Double_t *par);     // Fit function with polynomial scaling

  void setPolynomialCoefficients(const std::vector<double>& coeffs);  // Set polynomial coefficients
  void setPolyOrder(double order);// set the order for the background fit for mc_p_n_poly
   
  void fitDataPoly(TH1D *h_data, std::vector<double> initialParams ={1,1,0,0,1,1,-1} ); // Fit data using mc_p_n_poly_R
  void fitDataWithBG(TH1D *h_data,std::vector<double> initialParams ={1,1,0,0,1,1,-1});  // Fit data using mc_p_n_BG_R
  
};

#endif
