#include "Utility.h"
#include <iostream>

// Constructor
Utility::Utility() {
  // Constructor body (can be empty)
}

// Destructor
Utility::~Utility() {
  // Destructor body (can be empty)
}


// Method to scale and shift the histogram. Accounting for bins that may be out of the original histogram's range 
TH1D* Utility::ScaleAndShiftHistogram(TH1D *hist, double scale, double shift) {
  TH1D *scaledShiftedHist = (TH1D*)hist->Clone("scaledShiftedHist");
  // what if we shifted by only an integer value of bins
  // Loop over bins and apply scaling and shifting
  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
    Double_t binCenter = hist->GetBinCenter(i);
    Double_t shiftedBinCenter = binCenter - shift;  // Apply the shift
    Double_t Val = hist->Interpolate(shiftedBinCenter);  // Interpolate the value
    
    // Check if the shifted bin is within the histogram range
    if (shiftedBinCenter >= hist->GetXaxis()->GetXmin() && shiftedBinCenter <= hist->GetXaxis()->GetXmax()) {
      scaledShiftedHist->SetBinContent(i, scale * Val);  // Set scaled value

      // Propagate errors using the interpolated value
      Double_t Error = hist->GetBinError(hist->FindBin(shiftedBinCenter));  // Interpolated error
      scaledShiftedHist->SetBinError(i, scale * Error);  // Scale the error
    } else {
      // If the bin is outside the range, set content and error to 0
      scaledShiftedHist->SetBinContent(i, 0);
      scaledShiftedHist->SetBinError(i, 0);
    }
  }
  
  return scaledShiftedHist;
}// end scale and shift histogram 2


TH2* Utility::Shift2DHistogramY(const TH2* hist, double yShift) {
  if (!hist) {
    std::cerr << "Error: Input histogram is null!" << std::endl;
    return nullptr;
  }

  // Clone the original histogram to create a new shifted histogram
  TH2* shiftedHist = (TH2*) hist->Clone();
  shiftedHist->SetName(Form("%s_shifted", hist->GetName())); // Optionally give it a new name
  shiftedHist->Reset(); // Clear the contents to fill the shifted values

  int nBinsX = hist->GetNbinsX();
  int nBinsY = hist->GetNbinsY();

  double binContent, binError;

  for (int i = 1; i <= nBinsX; ++i) {
    for (int j = 1; j <= nBinsY; ++j) {
      double yValue = hist->GetYaxis()->GetBinCenter(j);
      double newYValue = yValue + yShift;

      // Find the bin corresponding to the shifted y-value
      int newYBin = hist->GetYaxis()->FindBin(newYValue);

      // Ensure the new bin is within the range
      if (newYBin >= 1 && newYBin <= nBinsY) {
	binContent = hist->GetBinContent(i, j);
	binError = hist->GetBinError(i, j);
	shiftedHist->SetBinContent(i, newYBin, binContent);
	shiftedHist->SetBinError(i, newYBin, binError);
      }
    }
  }

  return shiftedHist; // Return the new shifted histogram
}// end Shift2DHistogramY


void Utility::customizeGraph(TGraphErrors *graph, int markerStyle, int markerColor, double markersize) {
  int numPoints = graph->GetN();

  // Set marker style and color for each point
  for (int i = 0; i < numPoints; ++i) {
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerColor(markerColor);
    graph->SetMarkerSize(markersize);
  }
}// end customizeGraph

void Utility::customizeGraphMore(TGraphErrors *graph, int markerStyle, int markerColor, double markerSize, 
				 std::string graphTitle ="", std::string xAxisLabel="", std::string yAxisLabel="",
				 double TitleOffsetX = 1.4, double TitleOffsetY = 2, 
				 double LabelOffsetX = 0.01, double LabelOffsetY = 0.01) {

  // Set marker style, color, and size
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerColor(markerColor);
  graph->SetMarkerSize(markerSize);

  // Set graph title and axis labels
  graph->SetTitle(graphTitle.c_str());
  graph->GetXaxis()->SetTitle(xAxisLabel.c_str());
  graph->GetYaxis()->SetTitle(yAxisLabel.c_str());

  // Adjust axis title offsets to provide more space
  graph->GetXaxis()->SetTitleOffset(TitleOffsetX); // Adjust as needed
  graph->GetYaxis()->SetTitleOffset(TitleOffsetY); // Adjust as needed


  // Adjust axis label offsets to provide more space
  graph->GetXaxis()->SetLabelOffset(LabelOffsetX); // Adjust as needed
  graph->GetYaxis()->SetLabelOffset(LabelOffsetY); // Adjust as needed
}// end customizegraphmore

void Utility::SliceAndProjectHistogram_xMinxMax(TH2D* hist2D, const std::vector<double>& xMinimum,const std::vector<double>& xMaximum, std::vector<TH1D*>& histVector, std::string xAxisName, std::string yAxisName, std::string type) {
  // Clear the vector to ensure it's empty before filling
  histVector.clear();

  if (xMinimum.size()!=xMaximum.size())
    {
      cout<<"error in the xMin and xMax ranges. Not same size."<<endl;
    }

  cout<< "xMinimum.size() "<<xMinimum.size() <<" , xMaximum.size() "<<xMaximum.size()<<endl;
  
  // loop over the slices 
  for (size_t i = 0; i < xMinimum.size() ; ++i) { //-1
    double xMin = xMinimum[i] ;
    int binMin = hist2D ->GetXaxis()->FindBin(xMin)+1;

    double xMax = xMaximum[i];
    int binMax = hist2D ->GetXaxis()->FindBin(xMax)-1;

    // Define the name and title for the TH1D histogram
    // Format the histogram name to display 4 decimal places
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(4) << "_"<<xMin << "_to_" << xMax;
    std::string histName = yAxisName+"__"+xAxisName + stream.str() + "_" +type;
    TH1D *projY = hist2D->ProjectionY(histName.c_str(), binMin, binMax);
    histVector.push_back(projY); //

    // cout<<endl;
    // cout<<"slice id: "<<i<<endl;
    // cout<<"xMin = "<<xMin<<" , xMax = "<<xMax<<endl;
    // cout<< "FindBin(xMin) = "<<hist2D ->GetXaxis()->FindBin(xMin)<<", using binMin= "<<binMin<< " to be exclusive"<<endl;
    // cout<< "FindBin(xMax) = "<<hist2D ->GetXaxis()->FindBin(xMax)<<", using binMax= "<<binMax<< " to be exclusive"<<endl;

  }// end loop over slices
}// end SliceAndProjectHistogram_xMinxMax




void Utility::SliceAndProjectHistogram_xMinxMax_inclusiveMin(TH2D* hist2D, const std::vector<double>& xMinimum,const std::vector<double>& xMaximum, std::vector<TH1D*>& histVector, std::string xAxisName, std::string yAxisName, std::string type) {
  // Clear the vector to ensure it's empty before filling
  histVector.clear();

  if (xMinimum.size()!=xMaximum.size())
    {
      cout<<"error in the xMin and xMax ranges. Not same size."<<endl;
    }

  cout<< "xMinimum.size() "<<xMinimum.size() <<" , xMaximum.size() "<<xMaximum.size()<<endl;
  
  // loop over the slices 
  for (size_t i = 0; i < xMinimum.size() ; ++i) { //-1
    double xMin = xMinimum[i] ;
    int binMin = hist2D ->GetXaxis()->FindBin(xMin); // inclusive on xMin

    double xMax = xMaximum[i];
    int binMax = hist2D ->GetXaxis()->FindBin(xMax)-1;// exclusive on xMax
    // Define the name and title for the TH1D histogram
    // Format the histogram name to display 4 decimal places
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(4) << "_"<<xMin << "_to_" << xMax;
    std::string histName = yAxisName+"__"+xAxisName + stream.str() + "_" +type;
    TH1D *projY = hist2D->ProjectionY(histName.c_str(), binMin, binMax);
    histVector.push_back(projY); //

    // cout<<endl;
    // cout<<"slice id: "<<i<<endl;
    // cout<<"xMin = "<<xMin<<" , xMax = "<<xMax<<endl;
    // cout<< "FindBin(xMin) = "<<hist2D ->GetXaxis()->FindBin(xMin)<<", using binMin= "<<binMin<< " to be inclusive"<<endl;
    // cout<< "FindBin(xMax) = "<<hist2D ->GetXaxis()->FindBin(xMax)<<", using binMax= "<<binMax<< " to be exclusive"<<endl;

  }// end loop over slices
}// end SliceAndProjectHistogram_xMinxMax_inclusiveMin


// useful for coin stability
// making a vector of Y projections of the same size as the xMinimum and xMaximum vectors. 
void Utility::ProjectAndCloneHistogram(TH2D* hist2D, const std::vector<double>& xMinimum,const std::vector<double>& xMaximum, std::vector<TH1D*>& histVector, std::string xAxisName, std::string yAxisName, std::string type) {
  // Clear the vector to ensure it's empty before filling
  histVector.clear();

  // We're not actually going use the xMinimum and xMaximum beyond getting the vector sie
  if (xMinimum.size()!=xMaximum.size())
    {
      cout<<"error in the xMin and xMax ranges. Not same size."<<endl;
    }

  //cout<< "xMinimum.size() "<<xMinimum.size() <<" , xMaximum.size() "<<xMaximum.size()<<endl;
  
  // loop over the slices 
  for (size_t i = 0; i < xMinimum.size() ; ++i) { //-1
    
    // Define the name and title for the TH1D histogram
    // Format the histogram name to display 4 decimal places
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(4) << "_"<<i;
    std::string histName = yAxisName+"__"+xAxisName + stream.str() + "_" +type;
    
    TH1D *projY = hist2D->ProjectionY(histName.c_str());
    histVector.push_back(projY); //

    // cout<<endl;
    // cout<<"slice id: "<<i<<endl;
    // cout<<"xMin = "<<xMin<<" , xMax = "<<xMax<<endl;
    // cout<< "FindBin(xMin) = "<<hist2D ->GetXaxis()->FindBin(xMin)<<", using binMin= "<<binMin<< " to be exclusive"<<endl;
    // cout<< "FindBin(xMax) = "<<hist2D ->GetXaxis()->FindBin(xMax)<<", using binMax= "<<binMax<< " to be exclusive"<<endl;

  }// end loop over slices
}// end ProjectAndCloneHistogram


TH2D* Utility::CutXYRangeTH2D(TH2D* hist, double xmin, double xmax, double ymin, double ymax) {
  // Create a new histogram as a copy of the original
  TH2D* histCopy = (TH2D*)hist->Clone();
  histCopy->SetName((std::string(hist->GetName()) + "_cut").c_str());  // Rename to avoid conflict

  // Loop over all bins and set the content to zero if outside the y-range
  int nBinsX = hist->GetNbinsX();
  int nBinsY = hist->GetNbinsY();

  for (int i = 1; i <= nBinsX; i++) {
     double binCenterX = hist->GetXaxis()->GetBinCenter(i);
    for (int j = 1; j <= nBinsY; j++) {
      double binCenterY = hist->GetYaxis()->GetBinCenter(j);
      if (binCenterY < ymin || binCenterY > ymax || binCenterX < xmin || binCenterX>xmax) {
	histCopy->SetBinContent(i, j, 0);
      }
    }
  }

  return histCopy;
}// end CutYRangeTH2D


// Function to check if TF1 goes below zero in a given range
bool Utility::doesFunctionGoBelowZero(TF1* func, Double_t xMin, Double_t xMax, Int_t nSteps = 1000) {
  // Step size for sampling the function
  Double_t step = (xMax - xMin) / nSteps;
    
  // Iterate over the range
  for (Int_t i = 0; i <= nSteps; ++i) {
    Double_t x = xMin + i * step;
    Double_t y = func->Eval(x);  // Evaluate the function at point x
        
    if (y < 0) {
      std::cout<<std::endl<< "WARING WARNING"<<std::endl;
      std::cout << "Function goes below zero at x = " << x << ", f(x) = " << y << std::endl<<std::endl;
      return true;  // Return true if it goes below zero
    }
  }
    
  // If no negative value was found, return false
  //std::cout << "Function does not go below zero in the range [" << xMin << ", " << xMax << "]." << std::endl<<std::endl;
  return false;
}// end doesFunctionGoBelowZero

// calculate the mean on a vector of doubles 
double Utility::CalculateMean(const std::vector<double>& data){
  double sum = 0;
  for (double value : data) {
    sum += value;
  }
  return sum/data.size();
}

// calculate the standard deviation on a vector of doubles
double Utility::CalculatePopulationStDev(const std::vector<double>& data){
  double mean = CalculateMean(data);
  double numerator_sum = 0;
  for (double value : data){
    numerator_sum += pow( (value-mean),2);
  }
  double variance = numerator_sum/data.size();
  return sqrt(variance);
}

// calculate the standard deviation on a vector of doubles
// Using the "sample variance" which divides by (n-1)
double Utility::CalculateStDev(const std::vector<double>& data){
  double mean = CalculateMean(data);
  double numerator_sum = 0;
  for (double value : data){
    numerator_sum += pow( (value-mean),2);
  }
  double variance = numerator_sum/(data.size()-1);
  return sqrt(variance);
}

// calculate the weighted mean on a vector of doubles 
double Utility::CalculateWeightedMean(const std::vector<double>& data, const std::vector<double>& uncert){
  double numerator_sum=0;
  double sum_of_weights=0;
  double value =0;
  double weight = 0;
  if (data.size()!=uncert.size()){
    cout<<"data and uncert vectors not the same size in Calculate Weighted Mean"<<endl;
  }
  for (int i = 0;i<data.size();i++) {
    value = data[i];
    weight = 1/( pow(uncert[i],2));
    numerator_sum += value*weight;
    sum_of_weights += weight;
  }
  return numerator_sum/sum_of_weights;
}

// calculate the weighted stdev on a vector of doubles 
double Utility::CalculateWeightedStDev(const std::vector<double>& data, const std::vector<double>& uncert){
  double numerator_sum=0;
  double sum_of_weights=0;
  double value =0;
  double weight = 0;
  double mean_w = CalculateWeightedMean(data, uncert);
  if (data.size()!=uncert.size()){
    cout<<"data and uncert vectors not the same size in Calculate Weighted Mean"<<endl;
  }
  for (int i = 0;i<data.size();i++) {
    value = data[i];
    weight = 1/( pow(uncert[i],2));
    numerator_sum += weight*pow((value-mean_w),2);
    sum_of_weights += weight;
  }
  double variance = numerator_sum/sum_of_weights;
  return sqrt(variance);
}

TH1D* Utility::sumHistogramsWithPolynomial(TH1D* h1, TH1D* h2, TF1* poly) {
  if (h1->GetNbinsX() != h2->GetNbinsX() || 
      h1->GetXaxis()->GetXmin() != h2->GetXaxis()->GetXmin() ||
      h1->GetXaxis()->GetXmax() != h2->GetXaxis()->GetXmax()) {
    std::cerr << "Histograms must have the same binning and range!" << std::endl;
    return nullptr;
  }

  // Create a new histogram for the sum
  TH1D *h_sum = (TH1D*)h1->Clone("h_sum");
  h_sum->SetTitle("Sum of Histograms and Polynomial");
  h_sum->Reset();


  // Loop over bins and add the content of h1, h2, and the polynomial
  for (int i = 1; i <= h_sum->GetNbinsX(); ++i) {
    double bin_center = h_sum->GetBinCenter(i);
    double content = h1->GetBinContent(i) + h2->GetBinContent(i) + poly->Eval(bin_center);
    h_sum->SetBinContent(i, content);
  }

  return h_sum;
}// end sumHistogramsWithPolynomial

// This function takes a string in the form:
// "(x1,y1),(x2,y2),(x3,y3)...."
// and uses reference params to turn it into two vectors of:
// (x1,x2,x3...) and (y1,y2,y3)
void Utility::parseStringToVectors(const std::string& input,std::vector<double>& xMin, std::vector<double>& xMax) {
  std::stringstream ss(input);
  std::string token;

  // Clear the input vectors in case they contain previous data
  xMin.clear();
  xMax.clear();

  while (std::getline(ss, token, ')')) {
    // Find the opening parenthesis
    std::size_t openParen = token.find('(');
    if (openParen == std::string::npos) continue;

    // Remove parentheses from the token
    std::string cleanToken = token.substr(openParen + 1);
    cleanToken.erase(remove(cleanToken.begin(), cleanToken.end(), ')'), cleanToken.end());

    // Split the token into xMin and xMax values
    std::stringstream pairStream(cleanToken);
    std::string xMinStr, xMaxStr;

    std::getline(pairStream, xMinStr, ',');
    std::getline(pairStream, xMaxStr, ',');

    try {
      // Convert strings to doubles and add to respective vectors
      xMin.push_back(std::stod(xMinStr));
      xMax.push_back(std::stod(xMaxStr));
    } catch (const std::invalid_argument& e) {
      std::cerr << "Invalid argument: " << e.what() << std::endl;
    } catch (const std::out_of_range& e) {
      std::cerr << "Out of range: " << e.what() << std::endl;
    }
  }
}// end parseStringToVectors 


//// Get the title from the histogram and display it on a canvas. 
/// This expects the title to be in the form:
///  y_axis:x_axis {cut1&&cut2&&cut3....}
void Utility::printParsedTitle(const std::string& title,TString outputlocation,const std::string& modifier) {

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
  TCanvas* cuts_canvas = new TCanvas(Form("cuts_canvas_%s",modifier.c_str()), Form("Parsed Histogram Title %s",modifier.c_str()) , 1000, 600);
    
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
    
  // Optionally, save the canvas as an image
  cuts_canvas->SaveAs(Form("%s/global_cuts_%s.pdf", outputlocation.Data(),modifier.c_str()));
}// end PrintParsedTitle

TGraphErrors*  Utility::histogramToGraphErrors(TH1D *hist) {
  int numBins = hist->GetNbinsX();

  double x[numBins];
  double y[numBins];
  double ex[numBins]; // No x errors for simplicity
  double ey[numBins];

  // Fill arrays with histogram data
  for (int i = 0; i < numBins; ++i) {
    x[i] = hist->GetBinCenter(i + 1);
    y[i] = hist->GetBinContent(i + 1);
    ex[i] = 0; // No x errors for simplicity
    ey[i] = hist->GetBinError(i + 1);
  }

  // Create a TGraphErrors
  TGraphErrors *graph = new TGraphErrors(numBins, x, y, ex, ey);

  // // Create a canvas
  // TCanvas *canvas = new TCanvas("canvas", "Graph from Histogram", 800, 600);

  // // Draw the graph
  // graph->Draw("AP");

  // canvas->Update();
  // canvas->Draw();

  graph->SetTitle("");
  return graph;
}// end histogramToGraphErrors 

TH1D* Utility::GetResidualHistogram(TH1D* hist, TF1* fit) {
  // Create a new histogram for residuals
  TH1D* h_residual = new TH1D(Form("%s_residual", hist->GetName()), "", hist->GetNbinsX(), hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());

  // Loop over bins and calculate residuals
  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
    double binError = hist->GetBinError(i);
    double binCenter = hist->GetBinCenter(i);
    double binContent = hist->GetBinContent(i);
    double fitValue = fit->Eval(binCenter);
    double residual = binContent - fitValue;
    h_residual->SetBinContent(i, residual);
    h_residual->SetBinError(i, binError);
     //h_residual->SetBinError(i, sqrt(binContent));
  }

  return h_residual;
}// end GetResidualHistogram


 // Add the bin-by-bin errors in quarature and return a clone of hist1 with the new errors 
TH1D* Utility::AddBinErrorsInQuadrature(TH1D* hist1, TH1D* hist2, TH1D* hist3) {
    // Clone hist1 to create the result histogram
    TH1D* result = (TH1D*)hist1->Clone();
    
    // Loop over all bins
    for (int i = 1; i <= hist1->GetNbinsX(); ++i) {
        // Get errors from each histogram
        double error1 = hist1->GetBinError(i);
        double error2 = hist2->GetBinError(i);
        double error3 = hist3->GetBinError(i);
	// Should cout the bin locatioons to make sure everything is expected
        
        // Calculate the new error as the quadrature sum of the errors
	double newError = std::sqrt(error1 * error1+ error2 * error2 + error3 * error3);
	 //cout<<"old data err: "<< error1<<", hist2 err: "<<error2<<", hist3 err: "<<error3<<", new error: "<<newError<<endl;
        
        // Set the new error in the result histogram
        result->SetBinError(i, newError);
    }
    
    return result;
}// end AddBinErrosInQuadrature


/// Set the errors of hist1 with the errors of hist2
void Utility::SetBinErrors(TH1D* hist1, TH1D* hist2)
{
  for (int i = 1; i <= hist1->GetNbinsX(); ++i){
    // Get the bin error from hist2
    double error =  hist2->GetBinError(i);
    // set the bin error of hist1 to the error of hist2 
    hist1->SetBinError(i, error);
  }
  return;  
}// end SetBinErrors


TGraphErrors* Utility::createGraphFromFit(TH1D* hist, TF1* fit) {
  int numBins = hist->GetNbinsX();

  double x[numBins];
  double y[numBins];
  double ex[numBins]; // No x errors for simplicity
  double ey[numBins];

  // Fill arrays with fit function data
  for (int i = 0; i < numBins; ++i) {
    double binCenter = hist->GetBinCenter(i + 1);
    x[i] = binCenter;
    y[i] = fit->Eval(binCenter);
    ex[i] = 0; // No x errors for simplicity
    ey[i] = 0; // No y errors for simplicity
  }

  // Create and return a TGraphErrors
  TGraphErrors* graph = new TGraphErrors(numBins, x, y, ex, ey);
  graph->SetTitle("");
  return graph;
}// end createGraphFromFit

void Utility::adjustCanvas(TCanvas* canvas,
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
}// end adjustCanvas

void Utility::adjustPad(TPad* pad,
			double leftMargin = 0.15, double rightMargin = 0.05, 
			double bottomMargin = 0.01, double topMargin = 0.01) {
  // cout<<"left "<< leftMargin<<endl;
  // cout<<"right "<< rightMargin<<endl;
  // cout<<"bottom "<< bottomMargin<<endl;
  // cout<<"top "<< bottomMargin<<endl;
  // Set canvas margins
  pad->SetLeftMargin(leftMargin);
  pad->SetRightMargin(rightMargin);
  pad->SetBottomMargin(bottomMargin);
  pad->SetTopMargin(topMargin);
}// end adjustPad

void Utility::adjustPadVirtual(TVirtualPad* pad,
			       double leftMargin = 0.15, double rightMargin = 0.05, 
			       double bottomMargin = 0.01, double topMargin = 0.01) {
  // cout<<"left "<< leftMargin<<endl;
  // cout<<"right "<< rightMargin<<endl;
  // cout<<"bottom "<< bottomMargin<<endl;
  // cout<<"top "<< bottomMargin<<endl;
  // Set canvas margins
  pad->SetLeftMargin(leftMargin);
  pad->SetRightMargin(rightMargin);
  pad->SetBottomMargin(bottomMargin);
  pad->SetTopMargin(topMargin);
}// end adjustPadVirtual

void Utility::AdjustHistPadLabelOffset(TH1D* hist, TPad* pad, double xoffset = 0.02, double yoffset= 0.02, double textsize=0.03){
  double padHeight = pad->GetAbsHNDC();
  // Adjust axis label positions
  hist->GetXaxis()->SetLabelOffset(xoffset); // Adjust X-axis label offset
  hist->GetYaxis()->SetLabelOffset(yoffset); // Adjust Y-axis label offset
  hist->GetXaxis()->SetLabelSize(textsize / padHeight);  // Adjust the X-axis label size relative to the pad height
  hist->GetYaxis()->SetLabelSize(textsize / padHeight);  // Adjust the Y-axis label size relative to the pad height
}// end AdjustHistPadLabelOffset

void Utility::AdjustHistLabelOffset(TH1D* hist, double xoffset = 0.02, double yoffset= 0.02, double textsize=0.05){
  // Adjust axis label positions
  hist->GetXaxis()->SetLabelOffset(xoffset); // Adjust X-axis label offset
  hist->GetYaxis()->SetLabelOffset(yoffset); // Adjust Y-axis label offset
  hist->GetXaxis()->SetLabelSize(textsize);   // Adjust the X-axis label size
  hist->GetYaxis()->SetLabelSize(textsize);   // Adjust the Y-axis label size
}// end AdjustHistLabelOffset


void Utility::AdjustGraphLabelOffset(TGraphErrors* graph, double xoffset = 0.02, double yoffset= 0.02, double textsize=0.05){
  // Adjust axis label positions
  graph->GetXaxis()->SetLabelOffset(xoffset); // Adjust X-axis label offset
  graph->GetYaxis()->SetLabelOffset(yoffset); // Adjust Y-axis label offset
  graph->GetXaxis()->SetLabelSize(textsize);   // Adjust the X-axis label size
  graph->GetYaxis()->SetLabelSize(textsize);   // Adjust the Y-axis label size
}// end AdjustGraphLabelOffset

void Utility::AdjustFitLabelOffset(TF1* fit, double xoffset = 0.02, double yoffset= 0.02, double textsize=0.05){
  // Adjust axis label positions
  fit->GetXaxis()->SetLabelOffset(xoffset); // Adjust X-axis label offset
  fit->GetYaxis()->SetLabelOffset(yoffset); // Adjust Y-axis label offset
  fit->GetXaxis()->SetLabelSize(textsize);   // Adjust the X-axis label size
  fit->GetYaxis()->SetLabelSize(textsize);   // Adjust the Y-axis label size
}// end Adjust Fit LabelOffset


// useful for variable stability studies
std::string Utility::incrementRangeStringStepsize(double min, double max,double plus_or_minus_xmin, double plus_or_minus_xmax,double stepsize, int nDecimals) {
  std::ostringstream result;
  result << std::fixed << std::setprecision(nDecimals);
  bool first_output = true; // To handle comma placement

  double range = max - min;
  double min_range_percent = 50.0;

  int nSteps_xmin = static_cast<int>(std::floor(plus_or_minus_xmin*2 / stepsize));
   
  int nSteps_xmax = static_cast<int>(std::floor(plus_or_minus_xmax*2 / stepsize));
   
  
  for (int i = 0; i <= nSteps_xmin; ++i) {
    double new_min = min - plus_or_minus_xmin + i * stepsize;
    for (int j = 0; j <= nSteps_xmax; ++j) {
      double new_max = max - plus_or_minus_xmax + j * stepsize;
      double new_range = new_max - new_min;
            
      // Ensure new_range meets the minimum percentage of the original range
      if (new_range >= range * min_range_percent * 0.01) {
	if (!first_output) {
	  result << ", ";
	}
	result << "(" << new_min << "," << new_max << ")";
	first_output = false;
      }
    }
  }

  return result.str(); // Return the constructed string
}// end incerement range string stepsize


// useful for fit range stabiliy 
std::string Utility::incrementRangeStringStepsizeNotNested(double min, double max,double stepsize_xmin, double stepsize_xmax, int nSteps, int nDecimals) {
  std::ostringstream result;
  result << std::fixed << std::setprecision(nDecimals);
  bool first_output = true; // To handle comma placement

  std::vector<double> new_min_vector;
  std::vector<double> new_max_vector;
  
  for (int i = 0; i <= nSteps; ++i) {
    double new_min = min + i * stepsize_xmin;
    new_min_vector.push_back(new_min);
  }
  
  for (int j = 0; j <= nSteps; ++j) {
    double new_max = max + j * stepsize_xmax;
    new_max_vector.push_back(new_max);
  }

  for (int i = 0; i <= nSteps; i++){
    if (!first_output) {
      result << ", ";
    }
    result << "(" << new_min_vector[i] << "," << new_max_vector[i] << ")";
    first_output = false;
  }
  

  return result.str(); // Return the constructed string
}// end incerement range not nested 


void Utility::DrawLatexRsfLabels(TLatex latex, double Rsf, double Rsf_err, double ChiSq, double ndf, int nEntries_data, double upperleftX = 0.88, double upperleftY = 0.48, double spacing = 0.01, double textsize = 0.07, int allign =32 ){
  latex.SetNDC(true);
  latex.SetTextSize(textsize); // Increase text size
  latex.SetTextAlign(allign); // 32 is right
  latex.DrawLatex(upperleftX, upperleftY,Form("R_{sf} = %.5f #pm %.5f", Rsf, Rsf_err));
  latex.DrawLatex(upperleftX, (upperleftY-spacing), Form("#chi^{2}/ndf = %.2f/%.0f", ChiSq, ndf));
  latex.DrawLatex(upperleftX, (upperleftY-spacing*2), Form("nEntries: %i", nEntries_data));
}// end DrawLatexRsfLabels



void Utility::RemoveDuplicateSlices(std::vector<double>& xMinSlices, std::vector<double>& xMaxSlices) {
  // Ensure both vectors have the same size
  if (xMinSlices.size() != xMaxSlices.size()) {
    std::cerr << "Error: xMinSlices and xMaxSlices must have the same size!" << std::endl;
    return;
  }

  // Use a set to store the unique ranges
  std::set<std::pair<double, double>> uniqueRanges;

  // New vectors to store non-duplicate ranges
  std::vector<double> new_xMinSlices;
  std::vector<double> new_xMaxSlices;

  for (size_t i = 0; i < xMinSlices.size(); ++i) {
    std::pair<double, double> range(xMinSlices[i], xMaxSlices[i]);

    // Only add the range if it's not a duplicate
    if (uniqueRanges.find(range) == uniqueRanges.end()) {
      uniqueRanges.insert(range);
      new_xMinSlices.push_back(xMinSlices[i]);
      new_xMaxSlices.push_back(xMaxSlices[i]);
    }
    else{
      // std::cout<<"not adding "<<xMinSlices[i]<<", "<<xMaxSlices[i]<<std::endl;
    }
  }

  // Replace the original vectors with the new ones
  xMinSlices = new_xMinSlices;
  xMaxSlices = new_xMaxSlices;
}// end remove dupiilicate slices


// Establish hcal active area excluding N blks from edge, MC DB
std::vector<double> Utility::hcalaa_mc (int exblkN_x=1, int exblkN_y=1) {
  std::vector<double> hcalaa;
  // Positions (mc)
  double hcalposXi_mc = -2.655;    //m, distance from beam center to top of HCal w/75cm offset
  double hcalposXf_mc = 1.155;     //m, distance from beam center to bottom of HCal w/75cm offset
  double hcalposYi_mc = -0.92964;  //m, distance from beam center to opposite-beam side of HCal
  double  hcalposYf_mc = 0.92964;   //m, distance from beam center to beam side of HCal
  double hcalblk_div_h = 0.15494;  //m, horizontal center-to-center dist.
  double hcalblk_div_v = 0.15875;  //m, vertical center-to-center dist.
  double hcalaaXi = hcalposXi_mc + exblkN_x*hcalblk_div_v;
  double hcalaaXf = hcalposXf_mc - exblkN_x*hcalblk_div_v;
  double hcalaaYi = hcalposYi_mc + exblkN_y*hcalblk_div_h;
  double hcalaaYf = hcalposYf_mc - exblkN_y*hcalblk_div_h;
  hcalaa.push_back( hcalaaXi ); 
  hcalaa.push_back( hcalaaXf );
  hcalaa.push_back( hcalaaYi );
  hcalaa.push_back( hcalaaYf );
  return hcalaa;
}// end calculating the hcal active area


// returns the hcal positions as a vector. Xi, Xf, Yi, Yf
std::vector<double> Utility::hcalpos(){
  std::vector<double> hcalpos;
  double hcalposXi_mc = -2.655;    //m, distance from beam center to top of HCal w/75cm offset
  double hcalposXf_mc = 1.155;     //m, distance from beam center to bottom of HCal w/75cm offset
  double hcalposYi_mc = -0.92964;  //m, distance from beam center to opposite-beam side of HCal
  double  hcalposYf_mc = 0.92964;   //m, distance from beam center to beam side of HCal
  hcalpos.push_back( hcalposXi_mc ); 
  hcalpos.push_back( hcalposXf_mc );
  hcalpos.push_back( hcalposYi_mc );
  hcalpos.push_back( hcalposYf_mc );
  return hcalpos;
}

std::vector<double> Utility::hcalfid (double dxsig_p, double dxsig_n, double dysig, std::vector<double> hcalaa, double nsigma_p = 1, double nsigma_n = 1, double nsigma_dy=1) {
  std::vector<double> fid;
  double hcalx_t = hcalaa[0] + nsigma_p*dxsig_p;  // top margin (relevant for proton)
  double hcalx_b = hcalaa[1] - nsigma_n*dxsig_n;  // bottom margin (relevant for neutron)
  double hcaly_r = hcalaa[2] + nsigma_dy*dysig;    // right margin
  double hcaly_l = hcalaa[3] - nsigma_dy*dysig;    // left margin
  fid.push_back( hcalx_t ); 
  fid.push_back( hcalx_b );
  fid.push_back( hcaly_r );
  fid.push_back( hcaly_l );
  return fid;
}// end hcalfid


// Function to create a box of TLine objects
std::vector<TLine*> Utility::CreateBox(const std::vector<double>& coords, int color, int lineWidth) {
  // Check if the coordinates are valid (expecting 4 values: Xi, Xf, Yi, Yf)
  if (coords.size() != 4) {
    std::cerr << "Error: Coordinates vector must have 4 elements (Xi, Xf, Yi, Yf)." << std::endl;
    return {};
  }

  double Xi = coords[0];
  double Xf = coords[1];
  double Yi = coords[2];
  double Yf = coords[3];

  // Create a vector to store the lines
  std::vector<TLine*> lines;

  // Create horizontal lines
  TLine* LineXi = new TLine(Yi, Xi, Yf, Xi); // Horizontal line at Xi
  LineXi->SetLineWidth(lineWidth);
  LineXi->SetLineColor(color);
  lines.push_back(LineXi);

  TLine* LineXf = new TLine(Yi, Xf, Yf, Xf); // Horizontal line at Xf
  LineXf->SetLineWidth(lineWidth);
  LineXf->SetLineColor(color);
  lines.push_back(LineXf);

  // Create vertical lines
  TLine* LineYi = new TLine(Yi, Xi, Yi, Xf); // Vertical line at Yi
  LineYi->SetLineWidth(lineWidth);
  LineYi->SetLineColor(color);
  lines.push_back(LineYi);

  TLine* LineYf = new TLine(Yf, Xi, Yf, Xf); // Vertical line at Yf
  LineYf->SetLineWidth(lineWidth);
  LineYf->SetLineColor(color);
  lines.push_back(LineYf);

  return lines;
}// end create box



TEllipse* Utility::CreateEllipse(double x_center, double y_center, 
                              double radius_x, double radius_y, 
				 double angle = 0, int color = kRed, int LineStyle=1, int FillStyle =0) {
    // Create the ellipse with the specified parameters
    TEllipse* ellipse = new TEllipse(x_center, y_center, radius_x, radius_y, 0, 360, angle);
    
    // Set the ellipse style
    ellipse->SetLineColor(color);  // Set line color
    ellipse->SetFillStyle(FillStyle); // Set fill style (0 means no fill)

    // Optionally adjust other attributes like line width or style
    ellipse->SetLineWidth(2);
    ellipse->SetLineStyle(LineStyle); // 1 is sold line, 2 is dashed line. 
    
    // Return the ellipse so it can be drawn or manipulated externally
    return ellipse;
}// end create ellipse


void Utility::Fit_or_Replace_Pol0(TH1D* hist, std::pair<double, double> range) {
    // Check if the histogram exists
    if (!hist) {
        std::cerr << "Error: Histogram is null!" << std::endl;
        return;
    }

    // Retrieve the list of functions associated with the histogram
    TF1* existingFit = hist->GetFunction("pol0");

    // If there is already a "pol0" fit, remove it
    if (existingFit) {
        std::cout << "Removing existing pol0 fit..." << std::endl;
        hist->GetListOfFunctions()->Remove(existingFit);
        delete existingFit; // Free the memory for the old fit
    }

    // Define the pol0 function and fit the histogram between range.first and range.second
    hist->Fit("pol0", "R Q", "", range.first, range.second);

    
}// end Fit_or_Replace_Pol0



void Utility::Fit_or_Replace_Pol0_xmin(TH1D* hist, double xmin) {
    // Check if the histogram exists
    if (!hist) {
        std::cerr << "Error: Histogram is null!" << std::endl;
        return;
    }

    // Get the upper limit of the X-axis (xmax)
    double xmax = hist->GetXaxis()->GetXmax();

    // Retrieve the list of functions associated with the histogram
    TF1* existingFit = hist->GetFunction("pol0");

    // If there is already a "pol0" fit, remove it
    if (existingFit) {
        std::cout << "Removing existing pol0 fit..." << std::endl;
        hist->GetListOfFunctions()->Remove(existingFit);
        delete existingFit; // Free the memory for the old fit
    }

    // Define the pol0 function and fit the histogram between xmin and xmax
    hist->Fit("pol0", "R Q", "", xmin, xmax);

}// end Fit or Replace Pol0 xmin


// Removes a fit function from a histogram. 
void Utility::RemoveExistingFit(TH1D* hist) {
    // Check if the histogram exists
    if (!hist) {
        std::cerr << "Error: Histogram is null!" << std::endl;
        return;
    }

    // Retrieve the list of functions associated with the histogram
    TList* functions = hist->GetListOfFunctions();

    // Check if there are any functions associated with the histogram
    if (functions && functions->GetSize() > 0) {
      //std::cout << "Removing all fits from the histogram..." << std::endl;

        // Loop through all functions and remove them
        TObject* func = nullptr;
        TIterator* it = functions->MakeIterator();
        while ((func = it->Next())) {
            TF1* fit = dynamic_cast<TF1*>(func);
            if (fit) {
                functions->Remove(fit);
                delete fit; // Free the memory for the old fit
            }
        }

        delete it; // Clean up the iterator
    } else {
      // std::cout << "No fit functions to remove." << std::endl;
    }
}// end RemoveExistingFit


// Remove Zeros from a TGraphErrors 
void Utility::RemoveZeroPoints(TGraphErrors* graph) {
    if (!graph) {
        std::cerr << "Error: Graph is null!" << std::endl;
        return;
    }

    int nPoints = graph->GetN();  // Get the number of points in the graph

    // Create temporary vectors to store non-zero points
    std::vector<double> xVals, yVals, exVals, eyVals;

    for (int i = 0; i < nPoints; ++i) {
        double x, y;
        graph->GetPoint(i, x, y);

        if (y != 0) {  // Only keep points where y is non-zero
            xVals.push_back(x);
            yVals.push_back(y);
            exVals.push_back(graph->GetErrorX(i));
            eyVals.push_back(graph->GetErrorY(i));
        }
    }

    // Clear the existing points in the graph
    graph->Set(0);

    // Add the non-zero points back to the graph
    for (size_t i = 0; i < xVals.size(); ++i) {
        graph->SetPoint(i, xVals[i], yVals[i]);
        graph->SetPointError(i, exVals[i], eyVals[i]);
    }
}// end RemoveZeroPoints
