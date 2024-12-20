#ifndef UTILITY_H
#define UTILITY_H

#include <TH1.h>
#include <TF1.h>

class Utility {
private:

public:
  
  Utility();  // Constructor
  ~Utility();  // Destructor



  /* // Utility function for scaling and shifting, Trying to account for bin ranges*/
  TH1D* ScaleAndShiftHistogram(TH1D *hist, double scale, double shift);

  void FillHistogramFromVector(const std::vector<double>& data, TH1D* hist);
  
  // Add the bin-by-bin errors in quarature and return a clone of hist1 with the new errors 
  TH1D* AddBinErrorsInQuadrature(TH1D* hist1, TH1D* hist2, TH1D* hist3); 

  /// Set the errors of hist1 with the errors of hist2
  void SetBinErrors(TH1D* hist1, TH1D* hist2);

  
  // shift a 2d in Y
  TH2* Shift2DHistogramY(const TH2* hist, double yShift);

  /// get the "rate" over a given range on the histogram
  double GetRateFromHistogram(TH1* hist, double xMin, double xMax);

  // alternative method of calculating the rate using the gaussian+pol0 fit 
  std::pair<double, double> GetRateAndErrorFromFitOffset(TH1* hist, double c, double c_error);
  
  // alternative method of calculating the rate by taking the pol0 offset from the fit that is a combination of a gaussian and pol0
  std::vector<std::pair<double, double>> GetRatesAndErrorsFromTH2D(TH2* hist2D, std::vector<double> c, std::vector<double> c_error);
  
  // get the rates over a TH2D
  std::vector<double> GetRatesFromTH2D(TH2* hist2D, double xMin, double xMax);

  // takes a vector of doubles and divides every entry by a given scalar (double)
  std::vector<double> DivideVectorByScalar(std::vector<double>& vec, double scalar);
  std::vector<double> MultiplyVectorByScalar(std::vector<double>& vec, double scalar);

  std::vector<std::pair<double, double>> DividePairVectorByScalar(const std::vector<std::pair<double, double>>& vec, double scalar);
  std::vector<std::pair<double, double>> MultiplyPairVectorByScalar(const std::vector<std::pair<double, double>>& vec, double scalar);

  std::vector<double> ExtractPairElement( const std::vector<std::pair<double, double>>& vec, bool extractFirst);

  
  // Function to remove elements from a vector based on a list of indices
  std::vector<double> RemoveElementsByIndices(const std::vector<double>& input, const std::vector<int>& indices);
  
  // Utility fucntion to customize a graph 
  void customizeGraph(TGraphErrors *graph, int markerStyle, int markerColor, double markersize);
  void customizeGraphMore(TGraphErrors *graph, int markerStyle, int markerColor, double markerSize, 
			  std::string graphTitle ="", std::string  xAxisLabel="", std::string  yAxisLabel="",
			  double TitleOffsetX = 1.4, double TitleOffsetY = 2, 
			  double LabelOffsetX = 0.01, double LabelOffsetY = 0.01);

  // Function to create a TPolyMarker overlay for recolored points
  TPolyMarker* CreatePointColorOverlay(TGraphErrors* graph, const std::vector<int>& indices, Color_t color);

  // draw a gaussian on top of a background (pol0)
  void DrawGaussianPlusPol0(TH1D* hist, double amp, double mean, double sigma, double p0);

  // make a legend for a stack of histograms with the number of entries
  TLegend* CreateLegendFromStack(THStack* stack, const std::vector<std::string>& labels, const std::vector<std::string>& options, double x1 = 0.7, double y1 = 0.7, double x2 = 0.9, double y2 = 0.9);

  // take a vector of histograms and turns them into a THStack for easy plotting 
  THStack* CreateStackFromVector(const std::vector<TH1*>& histograms, const std::string& stackName = "stack", const std::string& stackTitle = "Histogram Stack");
  
  // ulity fuction to slice up a 2d histogram 
  void SliceAndProjectHistogram_xMinxMax(TH2D* hist2D, const std::vector<double>& xMinimum,const std::vector<double>& xMaximum, std::vector<TH1D*>& histVector, std::string xAxisName, std::string yAxisName, std::string type);
  void SliceAndProjectHistogram_xMinxMax_inclusiveMin(TH2D* hist2D, const std::vector<double>& xMinimum,const std::vector<double>& xMaximum, std::vector<TH1D*>& histVector, std::string xAxisName, std::string yAxisName, std::string type);

  void ProjectAndCloneHistogram(TH2D* hist2D, const std::vector<double>& xMinimum,const std::vector<double>& xMaximum, std::vector<TH1D*>& histVector, std::string xAxisName, std::string yAxisName, std::string type);

  // sets the bin content to zero outside the xmin xmax ymin ymax range
  TH2D* CutXYRangeTH2D(TH2D* hist, double xmin, double xmax, double ymin, double ymax);
  
  bool doesFunctionGoBelowZero(TF1* func, Double_t xMin, Double_t xMax, Int_t nSteps = 1000);

  // statisitcs on vectors of doubles 
  double CalculateMean(const std::vector<double>& data);
  double CalculateStDev(const std::vector<double>& data); // Using the "sample variance" which divides by (n-1)
 
  double CalculatePopulationStDev(const std::vector<double>& data);// Using the "population sample variance" which divides by (n-1)
  double CalculateWeightedMean(const std::vector<double>& data, const std::vector<double>& uncert);
  double CalculateWeightedStDev(const std::vector<double>& data, const std::vector<double>& uncert);

  TH1D* sumHistogramsWithPolynomial(TH1D* h1, TH1D* h2, TF1* poly);

  // converts a string to vetors 
  void parseStringToVectors(const std::string& input, std::vector<double>& xMin, std::vector<double>& xMax);
  // prints the parsed histogram title and makes a canvas 
  void printParsedTitle(const std::string& title,TString outputlocation,const std::string& modifier );
  /// converts a histogram to TGraphErrors
  TGraphErrors*  histogramToGraphErrors(TH1D *hist);
  /// Gets a residual between a histogram and a fit 
  TH1D* GetResidualHistogram(TH1D* hist, TF1* fit);
  // converts a TF1 to a TGraphErrors
  TGraphErrors* createGraphFromFit(TH1D* hist, TF1* fit);
  // adjusts a canvas 
  void adjustCanvas(TCanvas* canvas,
		    double leftMargin = 0.15, double rightMargin = 0.05, 
		    double bottomMargin = 0.15, double topMargin = 0.10);
  void adjustPad(TPad* pad,
		 double leftMargin = 0.15, double rightMargin = 0.05, 
		 double bottomMargin = 0.3, double topMargin = 0.10);
  void adjustPadVirtual(TVirtualPad* pad,
			double leftMargin = 0.15, double rightMargin = 0.05, 
			double bottomMargin = 0.01, double topMargin = 0.01);
  void AdjustHistLabelOffset(TH1D* hist, double xoffset = 0.02, double yoffset= 0.02, double textsize = 0.05);
  void AdjustHistPadLabelOffset(TH1D* hist, TPad* pad, double xoffset = 0.02, double yoffset= 0.02, double textsize=0.03);
  void AdjustGraphLabelOffset(TGraphErrors* graph, double xoffset = 0.02, double yoffset= 0.02, double textsize = 0.05);
  void AdjustFitLabelOffset(TF1* fit, double xoffset = 0.02, double yoffset= 0.02, double textsize=0.05);
  void DrawLatexRsfLabels(TLatex latex, double Rsf, double Rsf_err, double ChiSq, double ndf,  int nEntries_data,double upperleftX = 0.88, double upperleftY = 0.48, double spacing = 0.01, double textsize = 0.07, int allign =32 );
 
  
  
  std::string incrementRangeStringStepsize(double min, double max,double plus_or_minus_xmin, double plus_or_minus_xmax,double stepsize, int nDecimals);

  std::string incrementRangeStringStepsizeNotNested(double min, double max,double stepsize_xmin, double stepsize_xmax, int nSteps, int nDecimals);
  
  // removes duplicates in my range strings 
  void RemoveDuplicateSlices(std::vector<double>& xMinSlices, std::vector<double>& xMaxSlices);

  // calculate the active area boundary for hcal 
  std::vector<double> hcalaa_mc (int exblkN_x=1, int exblkN_y=1);
  
  // calculate the fiducal boundary given the number of sigma for each 
  std::vector<double> hcalfid (double dxsig_p, double dxsig_n, double dysig, std::vector<double> hcalaa, double nsigma_p = 1, double nsigma_n = 1, double nsigma_dy=1);
  // returns the hcal positions as a vector
  std::vector<double> hcalpos();
  
  // Function to create a box of TLine objects
  std::vector<TLine*> CreateBox(const std::vector<double>& coords, int color, int lineWidth);
  // function to creatue an ellipse 
  TEllipse* CreateEllipse(double x_center, double y_center, 
			  double radius_x, double radius_y, 
			  double angle = 0, int color = kRed, int LineStyle=1, int FillStyle=0);

  // takes a TH1D and replaces the fit with a pol0. (Helpful for draw_SpotCutStudy)
  void Fit_or_Replace_Pol0(TH1D* hist,std::pair<double, double> range );
  // takes xmin and uses the full range of the histo to get xmax 
  void Fit_or_Replace_Pol0_xmin(TH1D* hist, double xmin);
  void RemoveExistingFit(TH1D* hist);

  // Remove zeros from a TGraphErrors
  void RemoveZeroPoints(TGraphErrors* graph);
};

#endif
