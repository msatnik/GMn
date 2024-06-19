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

struct JobSummary {
  int jobid = -1; // job number on the batch farm
  int Nthrown = -1; // number of events that passed rejection sampling
  int Ntried = -1; // number of events simc had to make to produce Nthrown
  double genvol = -1.0; //(MeV*sr^2) I should go look up what this is
  double luminosity = -1.0; // (ub^-1) density/time of the beam (?) 
  double ebeam = -1.0; // (GeV) Energy of the thrown electrons 
  double charge = -1.0; // (mC) Charge of the beam (how does this compare to luminosity and energy?)
  int RndmSeed = -1; // seed used in simc generation
  int UsingRS = -1; // if rejection sampling was used in the simuation
  double MaxWtRS = -1.0; // max weight value given to the simulation to conduct the rejection sampling
  int wtGTmaxwt = -1; // number of events that were above that max weight
  double ObsMaxWtRS = -1.0; // the largest value of the weight that was above the max weight
};


void simc_read_in_test(){ //main


  //// ****** Reading in the simc summary file****************************************************
  // Declare file path 
  std::string summaryfilePath("/lustre19/expphy/volatile/halla/sbs/seeds/simc/p390sf_sbs4_sbs30p_simc/simcout/p390sf_sbs4_sbs30p_simc_deep_summary.csv");

  // Open the CSV file
  std::ifstream summaryfile(summaryfilePath);
  std::cout<<"reading in summary file: "<< summaryfilePath<<std::endl;

  // Check if the file is opened successfully
  if (!summaryfile.is_open()) {
    std::cerr << "Error opening file!" << std::endl;
    return 1;
  }

  // Skip the first line (header)
  std::string header;
  std::getline(summaryfile, header);

  std::vector<JobSummary> JobSummaries;

  // Read the file line by line
  std::string line;
  while (std::getline(summaryfile, line)) {
    // Skip empty lines
    if (line.empty())
      continue;

    // Create a string stream from the line
    std::istringstream iss(line);

    // Parse the line using comma as the delimiter
    std::string token;
    std::getline(iss, token, ',');  // Extrac 1st: jobid
    int jobid = std::stoi(token);

    std::getline(iss, token, ',');  // Extract 2nd param: Nthrown
    int Nthrown = std::stoi(token);

    std::getline(iss, token, ',');  // Extract 3nd param: Ntried
    int Ntried = std::stoi(token);

    std::getline(iss, token, ',');  // Extract 4th param: genvol
    double genvol = std::stod(token);

    std::getline(iss, token, ',');  // Extract 5th param: luminosity
    double luminosity = std::stod(token);

    std::getline(iss, token, ',');  // Extract 6th param: ebeam
    double ebeam = std::stod(token);

    std::getline(iss, token, ',');  // Extract 7th param: charge
    double charge = std::stod(token);
    
    std::getline(iss, token, ',');  // Extract 8th param: RndmSeed
    int RndmSeed= std::stoi(token);

    std::getline(iss, token, ',');  // Extract 9th param: UsingRS
    int UsingRS= std::stoi(token);

    std::getline(iss, token, ',');  // Extract 10th param: MaxWtRS
    double MaxWtRS= std::stod(token);

    std::getline(iss, token, ',');  // Extract 11th param: wtGTmaxwt
    int wtGTmaxwt= std::stoi(token);

    std::getline(iss, token, ',');  // Extract 12th param: ObsMaxWtRS
    double ObsMaxWtRS= std::stod(token);

    // Create a JobSummary object and add it to the vector
    JobSummary job;
    job.jobid = jobid;
    job.Nthrown = Nthrown;
    job.Ntried = Ntried;
    job.genvol = genvol;
    job.luminosity = luminosity;
    job.ebeam = ebeam;
    job.charge = charge;
    job.RndmSeed = RndmSeed;
    job.UsingRS = UsingRS;
    job.MaxWtRS = MaxWtRS;
    job.wtGTmaxwt = wtGTmaxwt;
    job.ObsMaxWtRS = ObsMaxWtRS;
    JobSummaries.push_back(job);
  } // end loop over each line of the summary file 

  // Close the file
  summaryfile.close();

  // Print out what was read in from the CSV file
  std::cout<<header<<endl;

  // // Loop over the vector 
  // for (size_t i = 0; i < JobSummaries.size(); ++i) {
  //   // Output values of the current entry
  //   // std::cout << "Job ID: " << JobSummaries[i].jobid <<", Nthrown: " << JobSummaries[i].Nthrown <<  ", Ntried: " << JobSummaries[i].Ntried << std::endl;
  //   std::cout << JobSummaries[i].jobid <<", "<< JobSummaries[i].Nthrown<<", "<< JobSummaries[i].Ntried <<", "<< JobSummaries[i].genvol <<", "<< JobSummaries[i].luminosity <<", "<< JobSummaries[i].ebeam <<", "<< JobSummaries[i].charge <<", "<< JobSummaries[i].RndmSeed <<", "<< JobSummaries[i].UsingRS <<", "<< JobSummaries[i].MaxWtRS <<", "<< JobSummaries[i].wtGTmaxwt <<", "<< JobSummaries[i].ObsMaxWtRS <<std::endl;
  // }

// Loop over the vector and output the values
    for (const auto& job : JobSummaries) {
	std::cout << job.jobid <<", "<< job.Nthrown<<", "<< job.Ntried <<", "<< job.genvol <<", "<< job.luminosity <<", "<< job.ebeam <<", "<< job.charge <<", "<< job.RndmSeed <<", "<< job.UsingRS <<", "<< job.MaxWtRS <<", "<< job.wtGTmaxwt <<", "<<job.ObsMaxWtRS <<std::endl;
    }// end loop over vector of jobs
 

    // loop over and get the sum of Ntried
    double Ntried_sum = 0;
    for (const auto& job : JobSummaries) {
      Ntried_sum += job.Ntried;
    }
    std::cout << "sum of Ntried: "<<Ntried_sum<<std::endl;
    
    // calculate the corrected weight 
    double max_weight =  JobSummaries[0].MaxWtRS; // MaxWtRS is the same for all jobs when runs are submitted together in Jlab-HPC
    double generation_volume = JobSummaries[0].genvol; // genvol is the same for all jobs when runs are submitted together in Jlab-HPC
    double luminosity = JobSummaries[0].luminosity; // luminosity is the same for all jobs when runs are submitted together in Jlab-HPC
    double corrected_weight = max_weight * generation_volume * luminosity / Ntried_sum;
    std::cout<<"corrected weight = max weight * generation volume * luminosity / Ntried_sum" <<std::endl;
    std::cout << corrected_weight <<" = " << max_weight <<" * " << generation_volume << " * "<<luminosity<<" / "<< Ntried_sum<<std::endl;
    

  return;
}//end main
