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
#include <algorithm>


void parser_test(){ //main

  TChain *C = new TChain("T");

  C->Add("/w/halla-scshelf2102/sbs/msatnik/GMn/symlinks/sbs4/sbs4_30p/*11496*");
  // C->Add("/w/halla-scshelf2102/sbs/msatnik/GMn/symlinks/sbs4/sbs4_30p/*11495*");

 // Create output file
    TFile outputFile("../output/parser_test.root", "RECREATE");

    double BBps_e;

    C->SetBranchStatus("bb.ps.e", 1);
    C->SetBranchAddress("bb.ps.e", &BBps_e);
    
    // Create output tree
    TTree *P = new TTree("P", "P");

    double bb_ps_e_out;
    double math_test_out;

     P->Branch("bb_ps_e", &bb_ps_e_out, "bb_ps_e/D");
     P->Branch("math_test", &math_test_out, "math_test/D");

     int counter = 0;
    
// Loop over events
    Long64_t nEntries = C->GetEntries();
    cout<<"nEntries: "<<nEntries<<endl;
    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
        // Get entry
        C->GetEntry(iEntry);

	if (iEntry%50000 == 0) cout<<iEntry<<endl;

	if(BBps_e < 0.1) continue;
	double math_test = BBps_e *2;
	if (math_test > 0.4) continue;


	math_test_out = math_test;
	bb_ps_e_out = BBps_e;

	counter++;
	P->Fill();

    }// end loop over entries




// Write output tree to file
    outputFile.cd();
    P->Write();
    
    // Close files
    outputFile.Close();

    cout<<"counter= "<<counter<<endl;

}//end main
