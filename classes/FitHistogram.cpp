#include "FitHistogram.h"
#include <TCanvas.h>
#include <TLegend.h>

//// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//// !!!!!!!!!!!!!! WARNING WARNING WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!
//// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//// The histograms that you fit using these functions (or at least the mc_p_n_poly2_R fit function) CAN NOT BE DRAWN. IT WILL CRASH. Make a clone of the histogram and a clone of the fit and draw those. To be extra safe, explictly make your own c1 canvas after you make the fit and draw something else on it so it doesn't accidentally draw it and crash. 


//// This class is intended to help with data - Monte Carlo comparison for GMn extraction
//// It allows you to fit the dx distribution from data to scaled version of montecarlo plus background,
//// This should remove the necessity of the simulation variables being globals

const int n_fitpoints = 300000;

// Constructor
FitHistogram::FitHistogram(TH1D *h_p, TH1D *h_n, double xmin, double xmax) : hist_p(h_p), hist_n(h_n), xMin(xmin), xMax(xmax), fitFunc(nullptr),scale_p(1),scale_p_err(0),R(1),R_err(0), ChiSq(0),NDF(0) {
  if (!hist_p || !hist_n) {
    std::cerr << "Error: One or both histograms are null in FitHistogram constructor!" << std::endl;
  }
  //std::cout << "In constructor" << std::endl;
}

// Destructor
FitHistogram::~FitHistogram() {
  if (fitFunc) {
    delete fitFunc;
    fitFunc = nullptr;
  }
  //std::cout << "In destructor" << std::endl;
}


// Define the fit function
Double_t FitHistogram::mc_p_n_poly_R(Double_t *x, Double_t *par) {
  
  // Check if histograms are valid
  if (!hist_p || !hist_n) {
    std::cerr << "Error: One of the histograms is null in mc_p_n_poly_R!" << std::endl;
    return 0;
  }
  // // Additional checks to confirm that hist_p and hist_n are not corrupted
  // std::cout << "hist_p: " << (hist_p ? "Valid" : "Null") << std::endl;
  // std::cout << "hist_n: " << (hist_n ? "Valid" : "Null") << std::endl;
  
  Double_t val = 0.0;

  // Get x value
  Double_t xx = x[0];

  // Retrieve parameters
  Double_t scale_p= par[0];
  Double_t R = par[1];
  Double_t shift_p = par[2];
  Double_t shift_n = par[3];

  std::vector<double> polyCoeffs={0};
  polyCoeffs.resize(polyorder+1,0); // Initialize to zero
    
  // Fill polynomial coefficients
  for (Int_t i = 0; i <=polyorder; ++i) {
    polyCoeffs[i] = par[i + 4];
  }

  double n = hist_n->Interpolate(xx-shift_n);
  double p = hist_p->Interpolate(xx-shift_p);
    
  // Calculate value using combination of histograms and polynomial background
  val = scale_p*(p + R*n); 

  // Add polynomial background
  for (Int_t i = 0; i <= polyorder; ++i) {
    val += polyCoeffs[i] * TMath::Power(xx, i);
  }

  return val;
}// end mc_p_n_poly_R

// Define the fit function with polynomial scaling (mc_p_n_BG_R)
Double_t FitHistogram::mc_p_n_BG_R(Double_t *x, Double_t *par) {
  // Check if histograms are valid
  if (!hist_p || !hist_n) {
    std::cerr << "Error: One of the histograms is null in mc_p_n_BG_R!" << std::endl;
    return 0;
  }

  Double_t xx = x[0];
  Double_t scale_p = par[0];  // Proton scale factor
  Double_t R = par[1];        // Neutron to proton ratio
  Double_t shift_p = par[2];
  Double_t shift_n = par[3];
  Double_t bg_scale = par[4]; // Scale factor for the polynomial background

  // Get proton and neutron values after shifting
  double n = hist_n->Interpolate(xx - shift_n);
  double p = hist_p->Interpolate(xx - shift_p);

  // Calculate the value using the scaled polynomial background
  Double_t val = scale_p * (p + R * n);

  // Scale and add the polynomial background
  for (size_t i = 0; i < polyCoefficients.size(); ++i) {
    val += bg_scale * polyCoefficients[i] * TMath::Power(xx,static_cast<Int_t>(i) );
  }

  return val;
}// end mc_p_n_BG_R

// Set polynomial coefficients
void FitHistogram::setPolynomialCoefficients(const std::vector<double>& coeffs) {
  polyCoefficients = coeffs;
}

// Set polynomial coefficients
void FitHistogram::setPolyOrder(double order) {
  polyorder = order;
}


// Method to fit data
void FitHistogram::fitDataPoly(TH1D *h_data, std::vector<double> initialParams ={1,1,0,0,1,1,-1} ) {

  if (!hist_p || !hist_n || !h_data) {
    std::cerr << "Error: One or both histograms are null in fitData!" << std::endl;
  }

  bool FitWithNoBackground = false;

  if (polyorder == -1)
    {
      polyorder = 0;
      FitWithNoBackground = true;
    }
  

  const int numpars = 4 + polyorder +1;
  
  // this utilizes a "lambda function". Google it to learn more
  fitFunc = new TF1("fitFuncPoly", [this](Double_t *x, Double_t *par) -> Double_t {
    return this->mc_p_n_poly_R(x, par);
  }, xMin, xMax, numpars);


  // Set initial parameter guesses and limits 
  //fitFunc->SetParameters(1,1,-0.05,-0.05,1,1,-1);
  fitFunc->SetParameters(&initialParams[0]);
  fitFunc->SetParLimits(0, 0, 10000);// par 0 limit, scale_p positive
  fitFunc->SetParLimits(1,0,2);// par 1 limit, R between 0 and 2 
  fitFunc->SetParLimits(2, -0.10, 0.10); // par 2 limit +- 10cm (proton shift)
  fitFunc->SetParLimits(3, -0.10, 0.10); // par 3 limit +- 10cm  (neutron shift)
  if(polyorder ==2){
    fitFunc->SetParLimits(6,-10000000,-0.0000000001); // x^2 term negative to force downward concavity
  }
  if(FitWithNoBackground){
    fitFunc->FixParameter(4, 0);// par 4 fix to zero. (pol0)
  }
  fitFunc->SetNpx(n_fitpoints);
  

  // Fit the data
  h_data->Fit(fitFunc,"R Q");

  // Save fit results
  scale_p = fitFunc ->GetParameter(0);
  scale_p_err = fitFunc ->GetParError(0);
  R = fitFunc ->GetParameter(1);
  R_err = fitFunc ->GetParError(1);
  shift_p = fitFunc ->GetParameter(2);
  shift_p_err = fitFunc ->GetParError(2); 
  shift_n = fitFunc ->GetParameter(3);
  shift_n_err = fitFunc ->GetParError(3);
  ChiSq = fitFunc->GetChisquare();
  NDF = fitFunc->GetNDF();

  // calculate scale_n from results
  scale_n = R*scale_p;
  scale_n_err = scale_n * sqrt( pow( (R_err / R), 2) + pow( (scale_p_err / scale_p),2) );
  // adding in quadrature isn't exactly right because they are correlated but it's fine because we don't really need to use it.

  poly_result.clear();

  // save poly results 
  for (int i =0 ; i <= polyorder ;i++)
    {
      poly_result.push_back( fitFunc->GetParameter(4+i));
      poly_result_err.push_back( fitFunc->GetParError(4+i));
    } 
}// end FitDataPoly



// Method to fit data
void FitHistogram::fitDataWithBG(TH1D *h_data,std::vector<double> initialParams ={1,1,0,0,1,1,-1}) {

  if (!hist_p || !hist_n || !h_data) {
    std::cerr << "Error: One or both histograms are null in fitData!" << std::endl;
  }

  // this utilizes a "lambda function". Google it to learn more
  fitFunc = new TF1("fitFuncBG", [this](Double_t *x, Double_t *par) -> Double_t {
    return this->mc_p_n_BG_R(x, par);
  }, xMin, xMax, 5);


  // Set initial parameter guesses and limits 
  //fitFunc->SetParameters(1,1,0,0,1,1,-1);
  fitFunc->SetParameters(&initialParams[0]);
  fitFunc->SetParLimits(0, 0, 10000);// par 0 limit, scale_p positive
  fitFunc->SetParLimits(1,0,2);// par 1 limit, R between 0 and 2 
  fitFunc->SetParLimits(2, -0.10, 0.10); // par 2 limit +- 10cm (proton shift)
  fitFunc->SetParLimits(3, -0.10, 0.10); // par 3 limit +- 10cm  (neutron shift)
  fitFunc->SetParLimits(4, 0, 10000);// par 4 limit, background scale positive
  fitFunc->SetNpx(n_fitpoints);

  // Fit the data
  h_data->Fit(fitFunc,"R Q");

  // Save fit results
  scale_p = fitFunc ->GetParameter(0);
  scale_p_err = fitFunc ->GetParError(0);
  R = fitFunc ->GetParameter(1);
  R_err = fitFunc ->GetParError(1);
  shift_p = fitFunc ->GetParameter(2);
  shift_p_err = fitFunc ->GetParError(2); 
  shift_n = fitFunc ->GetParameter(3);
  shift_n_err = fitFunc ->GetParError(3);
  BGscale = fitFunc->GetParameter(4);
  BGscale_err = fitFunc->GetParError(4);
  ChiSq = fitFunc->GetChisquare();
  NDF = fitFunc->GetNDF();

  // calculate scale_n from results
  scale_n = R*scale_p;
  scale_n_err = scale_n * sqrt( pow( (R_err / R), 2) + pow( (scale_p_err / scale_p),2) );
  // adding in quadrature isn't exactly right because they are correlated but it's fine because we don't really need to use it.
  
}// end FitDataWithBG



