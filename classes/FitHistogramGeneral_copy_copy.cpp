#include "FitHistogramGeneral.h"
#include <TCanvas.h>
#include <TLegend.h>

//// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//// !!!!!!!!!!!!!! WARNING WARNING WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!
//// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//// The histograms that you fit using these functions CANNOT BE DRAWN. IT WILL CRASH. Make a clone of the histogram and a clone of the fit and draw those. To be extra safe, explictly make your own c1 canvas after you make the fit and draw something else on it so it doesn't accidentally draw it and crash. 

//// This class is intended to help with data - Monte Carlo comparison for GMn extraction
//// It fits a histogram to two other histograms that are scaled. (No background). Inteded for comparing inelastic MC + elastic MC to data for W2. 

//// This should remove the necessity of the histgrams you are using for the fit having to be global variables 

// Constructor
FitHistogramGeneral::FitHistogramGeneral(TH1D *h_1, TH1D *h_2, double xmin, double xmax) 
  : hist1(h_1), hist2(h_2), xMin(xmin), xMax(xmax), fitFunc(nullptr), scale1(0), scale1_err(0), 
    scale2(0), scale2_err(0), ChiSq(0), NDF(0) {

  if (!hist1) {
    std::cerr << "Error: hist1 is null in FitHistogramGeneral constructor!" << std::endl;
  }
  if (h_1 && !h_1->GetNbinsX()) {
    std::cerr << "Error: hist1 has no bins in FitHistogramGeneral constructor!" << std::endl;
    hist1 = nullptr;
  }

  if (!hist2) {
    std::cerr << "(may be okay given situation): hist2 is null in FitHistogramGeneral constructor!" << std::endl;
  }
  if (h_2 && !h_2->GetNbinsX()) {
    std::cerr << "Error: hist2 has no bins in FitHistogramGeneral constructor!" << std::endl;
    hist2 = nullptr;
  }
}// end constructor 


// Destructor
FitHistogramGeneral::~FitHistogramGeneral() {
  if (fitFunc) {
    delete fitFunc;
    fitFunc = nullptr;
  }
}// end Destructor



void FitHistogramGeneral::Reset() {
  if (fitFunc) {
    delete fitFunc;
    fitFunc = nullptr;
  }

  hist1 = nullptr;
  hist2 = nullptr;
  xMin = 0;
  xMax = 0;
  scale1 = 0;
  scale1_err = 0;
  scale2 = 0;
  scale2_err = 0;
  ChiSq = 0;
  NDF = 0;
}

void FitHistogramGeneral::Set(TH1D *h_1, TH1D *h_2, double xmin, double xmax) {
  if (fitFunc) {
    delete fitFunc;
    fitFunc = nullptr;
  }

  // Check if the histograms are valid
  if (!h_1) {
    std::cerr << "Error: hist1 is null in Set!" << std::endl;
    hist1 = nullptr;
  } else if (h_1->GetNbinsX() == 0) {
    std::cerr << "Error: hist1 has no bins in Set!" << std::endl;
    hist1 = nullptr;
  } else {
    hist1 = h_1;
  }

  if (!h_2) {
    std::cerr << "(may be okay): hist2 is null in Set!" << std::endl;
    hist2 = nullptr;
  } else if (h_2->GetNbinsX() == 0) {
    std::cerr << "Error: hist2 has no bins in Set!" << std::endl;
    hist2 = nullptr;
  } else {
    hist2 = h_2;
  }

  xMin = xmin;
  xMax = xmax;
  scale1 = 0;
  scale1_err = 0;
  scale2 = 0;
  scale2_err = 0;
  ChiSq = 0;
  NDF = 0;
}


void FitHistogramGeneral::fitData(TH1D *h_data) {
  if (!hist1 || !h_data) {
    std::cerr << "Error: hist1 or h_data is null in fitData!" << std::endl;
    return;
  }

  // Additional checks for histogram validity
  if (h_data->GetNbinsX() == 0) {
    std::cerr << "Error: h_data has no bins in fitData!" << std::endl;
    return;
  }

  if (hist1->GetNbinsX() == 0) {
    std::cerr << "Error: hist1 has no bins in fitData!" << std::endl;
    return;
  }

  if (hist2 && hist2->GetNbinsX() == 0) {
    std::cerr << "Warning: hist2 has no bins in fitData!" << std::endl;
    hist2 = nullptr;  // Treat as null
  }

  int nParams = (hist2) ? 2 : 1;  // Use 2 parameters if hist2 is valid, otherwise use only 1

  // Create fit function with a lambda to use hist1_hist2_sum
  fitFunc = new TF1("fitFunc", [this](Double_t *x, Double_t *par) -> Double_t {
    return this->hist1_hist2_sum(x, par);
  }, xMin, xMax, nParams);

  // Set initial parameter guesses and limits
  fitFunc->SetParameter(0, 1);  // Initial guess for scale1
  fitFunc->SetParLimits(0, 0, 100000);  // par 0 limit, scale1 positive
  
  if (hist2) {
    fitFunc->SetParameter(1, 1);  // Initial guess for scale2
    fitFunc->SetParLimits(1, 0, 100000);  // par 1 limit, scale2 positive
  }

  // Fit the data
  h_data->Fit(fitFunc, "R Q");

  // Save fit results
  scale1 = fitFunc->GetParameter(0);
  scale1_err = fitFunc->GetParError(0);
  
  if (hist2) {
    scale2 = fitFunc->GetParameter(1);
    scale2_err = fitFunc->GetParError(1);
  } else {
    scale2 = 0;
    scale2_err = 0;
  }

  ChiSq = fitFunc->GetChisquare();
  NDF = fitFunc->GetNDF();
}
