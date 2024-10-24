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
FitHistogramGeneral::FitHistogramGeneral(TH1D *h_1, TH1D *h_2, double xmin, double xmax) : hist1(h_1), hist2(h_2), xMin(xmin), xMax(xmax), fitFunc(nullptr),scale1(0),scale1_err(0), scale2(0),scale2_err(0),ChiSq(0),NDF(0) {
  if (!hist1) {
    std::cerr << "Error: hist1 is null in FitHistogramGeneral constructor!" << std::endl;
  }
  if (!hist2) {
    std::cerr << "(may be okay given situation): hist2 is null in FitHistogramGeneral constructor!" << std::endl;
  }
  //std::cout << "In constructor" << std::endl;
}

// Destructor
FitHistogramGeneral::~FitHistogramGeneral() {
  if (fitFunc) {
    delete fitFunc;
    fitFunc = nullptr;
  }
  //std::cout << "In destructor" << std::endl;
}


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
}// end Reset())

void FitHistogramGeneral::Set(TH1D *h_1, TH1D *h_2, double xmin, double xmax) {
   if (fitFunc) {
    delete fitFunc;
    fitFunc = nullptr;
  }
  hist1 = h_1;
  hist2 = h_2;
  xMin = xmin;
  xMax = xmax;
  scale1 = 0;
  scale1_err = 0;
  scale2 = 0;
  scale2_err = 0;
  ChiSq = 0;
  NDF = 0;
}// end Set


Double_t FitHistogramGeneral::hist1_hist2_sum(Double_t *x, Double_t *par) {
  
  // Check if histograms are valid
  if (!hist1) {
    std::cerr << "Error: hist1 is null in hist1_hist2_sum!" << std::endl;
    return 0;
  }

  Double_t val = 0.0;
  Double_t xx = x[0];
  Double_t scale1 = par[0];
  Double_t scale2 = (hist2) ? par[1] : 0;  // If hist2 is null, set scale2 to 0

  int bin1 = hist1->FindBin(xx);
  double y1 = hist1->GetBinContent(bin1);
  double y2 = (hist2) ? hist2->GetBinContent(hist2->FindBin(xx)) : 0;  // Only use hist2 if it's valid

  // Calculate value using hist1 and (optionally) hist2
  val = scale1 * y1 + scale2 * y2;

  return val;
}


// // Define the fit function
// Double_t FitHistogramGeneral::hist1_hist2_sum(Double_t *x, Double_t *par) {
  
//   // Check if histograms are valid
//   if (!hist1 || !hist2) {
//     std::cerr << "Error: One of the histograms is null in elas_inel_sum!" << std::endl;
//     return 0;
//   }
//   // // Additional checks to confirm that hist1 and hist2 are not corrupted
//   // std::cout << "hist1: " << (hist1 ? "Valid" : "Null") << std::endl;
//   // std::cout << "hist2: " << (hist2 ? "Valid" : "Null") << std::endl;
  
//   Double_t val = 0.0;

//   // Get x value
//   Double_t xx = x[0];

//   // Retrieve parameters
//   Double_t scale1= par[0];
//   Double_t scale2 = par[1];

//   //double y1 = hist1->Interpolate(xx);
//   //double y2 = hist2->Interpolate(xx);

//   int bin1 = hist1->FindBin(xx);
//   int bin2 = hist2->FindBin(xx);
//   double y1 = hist1->GetBinContent(bin1);
//   double y2 = hist2->GetBinContent(bin2);
    
//   // Calculate value using combination of histograms and letting them scale up
//   val = scale2*y2 + scale1*y1;

//   return val;
// }// end hist1_hist2_sum



// Fit the data. This allows for histogram 2 to be null. 
void FitHistogramGeneral::fitData(TH1D *h_data) {
  if (!hist1 || !h_data) {
    std::cerr << "Error: hist1 or h_data is null in fitData!" << std::endl;
    return;
  }

  int nParams = (hist2) ? 2 : 1;  // Use 2 parameters if hist2 is valid, otherwise use only 1

  fitFunc = new TF1("fitFunc", [this](Double_t *x, Double_t *par) -> Double_t {
    return this->hist1_hist2_sum(x, par);
  }, xMin, xMax, nParams);
  fitFunc->SetNpx(10000000);
  // Set initial parameter guesses and limits
  fitFunc->SetParameter(0, 1);
  fitFunc->SetParLimits(0, 0, 100000);  // par 0 limit, scale1 positive
  
  if (hist2) {
    fitFunc->SetParameter(1, 1);
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
}// end fitData



// // Method to fit data
// void FitHistogramGeneral::fitData(TH1D *h_data) {

//   if (!hist1 || !hist2 || !h_data) {
//     std::cerr << "Error: One or both histograms are null in fitData!" << std::endl;
//   }

//   // this utilizes a "lambda function". Google it to learn more
//   fitFunc = new TF1("fitFunc", [this](Double_t *x, Double_t *par) -> Double_t {
//     return this->hist1_hist2_sum(x, par);
//   }, xMin, xMax, 2);


//   // Set initial parameter guesses
//   fitFunc->SetParameters(1,1);
//   fitFunc->SetParLimits(0, 0, 100000);// par 0 limit, scale1 positive
//   fitFunc->SetParLimits(1, 0, 100000);// par 1 limit, scale2 positive
//   fitFunc->SetNpx(500);

//   // Fit the data
//   h_data->Fit(fitFunc,"R Q");

//   std::cout << "Integral of hist1: " << hist1->Integral() << std::endl;
//   std::cout << "Integral of hist2: " << hist2->Integral() << std::endl;
//   std::cout << "Integral of h_data: " << h_data->Integral() << std::endl;
  
//   // Save fit results
//   scale1 = fitFunc ->GetParameter(0);
//   scale1_err = fitFunc ->GetParError(0);
//   scale2 = fitFunc ->GetParameter(1);
//   scale2_err = fitFunc ->GetParError(1);

//   // std::cout<<"fit results in class: "<< std::endl;
//   // std::cout<<"scale1 = "<<scale1<<" +/- "<<scale1_err<<endl;
//   // std::cout<<"scale2 = "<<scale2<<" +/- "<<scale2_err<<endl;
  
//   ChiSq = fitFunc->GetChisquare();
//   NDF = fitFunc->GetNDF();

// }// end FitData





