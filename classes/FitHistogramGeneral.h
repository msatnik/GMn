#ifndef FITHISTOGRAMGENERAL_H
#define FITHISTOGRAMGENERAL_H

#include <TH1.h>
#include <TF1.h>

class FitHistogramGeneral {
private:

public:
  TH1D *hist1;  //  Histogram 1 
  TH1D *hist2; // Histogram 2
  TF1 *fitFunc; // Custom fit function
  double xMin; // min x on fit range
  double xMax; // max x on fit range
  double scale1; // inelastic scale factor from fit
  double scale1_err; //
  double scale2; // elastic scale factor from fit
  double scale2_err;
  double ChiSq; // chi^2 from fit
  double NDF; // number of degrees of freedom from fit.

    
  FitHistogramGeneral(TH1D *hist_p, TH1D *hist_n, double xmin, double xmax);  // Constructor
  ~FitHistogramGeneral();  // Destructor
  
  void Reset();
  void Set(TH1D *h_1, TH1D *h_2, double xmin, double xmax);
  

  
  Double_t hist1_hist2_sum(Double_t *x, Double_t *par);  // Fit function for data-mc comparison

  void fitData(TH1D *h_data);  // Fit data using mc_p_n_poly_R
 
};

#endif
