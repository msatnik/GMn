#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>

void draw_dx_with_different_cuts() {

  gStyle->SetOptFit(11111);
  gStyle->SetCanvasPreferGL(1);
  //gStyle -> SetOptStat(0110);
  gStyle ->SetEndErrorSize(0);
  
  //const char* filename="/lustre24/expphy/volatile/halla/sbs/msatnik/output/sbs8_70p_simc_deep_sept3.root";
  const char* filename="/lustre24/expphy/volatile/halla/sbs/msatnik/output/sbs8_50p_simc_deep_sept3.root";
  //const char* filename="/lustre24/expphy/volatile/halla/sbs/msatnik/output/sbs9_70p_simc_deep_sept12_z.root";
  
  // Open the ROOT file
  TFile *file = TFile::Open(filename);
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }

  // Access the tree (assuming it's named "P")
  TTree *P = (TTree*)file->Get("P");
  if (!P) {
    std::cerr << "Tree P not found in file!" << std::endl;
    return;
  }

 
    
  // Create histograms
  TH1D *orig = new TH1D("orig", "Original Distribution;hcal_dx;Counts", 140, -2, 1.5);
  TH1D *H_arm = new TH1D("H_arm", "H Arm Cuts;hcal_dx;Counts", 140, -2, 1.5);
  TH1D *H_arm_anti_y_exp = new TH1D("H_arm_anti_y_exp", "H Arm Cuts w/ anti y_exp;hcal_dx;Counts", 140, -2, 1.5);
  TH1D *e_arm = new TH1D("e_arm", "e Arm Cuts;hcal_dx;Counts", 140, -2, 1.5);
  TH1D *both = new TH1D("both", "Both Arm Cuts;hcal_dx;Counts", 140, -2, 1.5);


  TH1D *H_arm_anticut = new TH1D("H_arm_anticut", "H Arm Anticut;hcal_dx;Counts", 140, -2, 1.5);
  TH1D *e_arm_anticut = new TH1D("e_arm_anticut", "e Arm Anticut;hcal_dx;Counts", 140, -2, 1.5);
   TH1D *both_anticut = new TH1D("both_anticut", "Both Arm Anticut;hcal_dx;Counts", 140, -2, 1.5);

  // // Create histograms
  // TH1D *orig = new TH1D("orig", "Original Distribution;hcal_dx;Counts", 400, -4, 4);
  // TH1D *H_arm = new TH1D("H_arm", "H Arm Cuts;hcal_dx;Counts", 400, -4, 4);
  // TH1D *e_arm = new TH1D("e_arm", "E Arm Cuts;hcal_dx;Counts", 400, -4, 4);
  // TH1D *both = new TH1D("both", "Both Arm Cuts;hcal_dx;Counts", 400, -4, 4);


    // Setting up cuts 
  std::string EnergyCutString = "bb_ps_e>0.2"; // PS and Shower Energy cut
  std::string TrackQualityCutString = "bb_tr_n==1&&bb_gem_track_nhits>=3 && bb_gem_track_chi2ndf<30 "; 
  std::string TargetVertexCutString = "abs(bb_tr_vz)<0.06";
  std::string W2CutString = "W2>0.66&&W2<1.10"; 
  std::string FidXCutString = "hcal_x_exp>-1.5&&hcal_x_exp<0.63";
  std::string FidYCutString = "hcal_y_exp>-0.4&&hcal_y_exp<0.35";
  std::string dyCutString = "abs(hcal_dy)<0.3";
  std::string e_over_p_CutString = " abs(e_over_p - 0.98) < 0.16";
  std::string HCal_Energy_CutString = "hcal_e>0.04";
  std::string Optics_CutString = "abs(bb_tr_r_x-bb_tr_r_th*0.9)<0.3 && abs(bb_tr_r_y-0.9*bb_tr_r_ph+0.005)<0.1";
   

  // Define cut strings
  std::string H_arm_string ="("+FidXCutString +"&&"+ FidYCutString+"&&"+dyCutString+"&&" +HCal_Energy_CutString+")";

  std::string e_arm_string ="("+TrackQualityCutString+"&&"+TargetVertexCutString+"&&"+EnergyCutString+"&&"+W2CutString+"&&"+e_over_p_CutString+"&&"+Optics_CutString+")";



  
  //"hcal_x_exp>-1.27&& hcal_x_exp<0.63
    
  //hcal_x_exp>-1.27&&hcal_x_exp<0.63&&hcal_y_exp>-0.4&&hcal_y_exp<0.35 && 

    
  // Combine the cuts for the "both" histogram
  std::string both_string = e_arm_string + " && " + H_arm_string;


  // Fill the histograms based on your cuts
  P->Draw(("hcal_dx>>orig"),"corrected_weight");
  P->Draw(("hcal_dx>>H_arm"), Form("corrected_weight*(%s)",H_arm_string.c_str()));
  P->Draw(("hcal_dx>>H_arm_anticut"), Form("corrected_weight*(!%s)",H_arm_string.c_str()));
  P->Draw(("hcal_dx>>e_arm"), Form("corrected_weight*(%s)",e_arm_string.c_str()));
  P->Draw(("hcal_dx>>e_arm_anticut"), Form("corrected_weight*(!%s)",e_arm_string.c_str()));
  P->Draw(("hcal_dx>>both"), Form("corrected_weight*(%s)",both_string.c_str()));
  P->Draw(("hcal_dx>>both_anticut"), Form("corrected_weight*(!(%s))",both_string.c_str()));
  
    
  // Set colors for the histograms
  orig ->SetLineWidth(2);
  orig  ->SetFillColorAlpha(kBlue,0.1);
  e_arm->SetLineColor(kMagenta);
  e_arm->SetLineWidth(2);
  e_arm ->SetFillColorAlpha(kMagenta,0.1);
  e_arm_anticut->SetLineColor(kRed);
  e_arm_anticut->SetLineWidth(2);
  e_arm_anticut ->SetFillColorAlpha(kRed,0.1);
  H_arm->SetLineColor(kGreen+2);
  H_arm->SetLineWidth(2);
  H_arm ->SetFillColorAlpha(kGreen+2,0.1);
  H_arm_anticut->SetLineColor(kRed);
  H_arm_anticut->SetLineWidth(2);
  H_arm_anticut ->SetFillColorAlpha(kRed,0.1);
  both->SetLineColor(kRed);
  both->SetLineWidth(2);
  both ->SetFillColorAlpha(kRed,0.1);
  both_anticut->SetLineColor(kRed+3);
  both_anticut->SetLineWidth(2);
  both_anticut ->SetFillColorAlpha(kRed+3,0.1);

  orig->Scale(1/3.0);

  // Draw histograms on a canvas
  TCanvas *c1 = new TCanvas("c1", "Comparison of Cuts", 800, 600);
  c1->SetGrid();
  orig->Draw("E hist");             // Draw the original histogram
  H_arm->Draw("E same hist");       // Draw H_arm on top
  e_arm->Draw("E same hist");       // Draw e_arm on top
  both->Draw("E same hist");        // Draw both on top


  // Add a legend
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(orig, "Original", "l f");
  legend->AddEntry(e_arm, "e Arm", "l f");
  legend->AddEntry(H_arm, "H Arm", "l ");
  legend->AddEntry(both, "Both", "l f");
  legend->Draw();

  // Save the canvas to a file
  // c1->SaveAs("hcal_dx_comparison.png");


  // Create the second canvas: Divided canvas showing each histogram separately
  TCanvas *c2 = new TCanvas("c2", "Individual Histograms", 800, 800);
  c2->Divide(2, 2);  // 2x2 grid for 4 histograms

  // Draw each histogram in its own pad
  c2->cd(1);   // Go to pad 1
  orig->Draw("E hist");
  TLegend *legend2_1 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend2_1->AddEntry(orig, "Original", "l f");
  legend2_1->Draw();

  c2->cd(2);   // Go to pad 2
  e_arm->Draw("E hist");
  e_arm_anticut->Draw("E hist same");
  TLegend *legend2_2 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend2_2->AddEntry(e_arm, "e Arm", "l f");
  legend2_2->AddEntry(e_arm_anticut, "e Arm Anticut", "l f");
  legend2_2->Draw();

  c2->cd(3);   // Go to pad 3
  H_arm->Draw("E hist");
  H_arm_anticut->Draw("E hist same");
  TLegend *legend2_3 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend2_3->AddEntry(H_arm, "H Arm", "l f");
  legend2_3->AddEntry(H_arm_anticut, "H Arm Anticut", "l f");
  legend2_3->Draw();

  c2->cd(4);   // Go to pad 4
  both->Draw("E hist");
  both_anticut->Draw("E hist same");
  TLegend *legend2_4 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend2_4->AddEntry(both, "Both", "l f");
  legend2_4->AddEntry(both_anticut, "Both Anticut", "l f");
  legend2_4->Draw();

  // Close the file
  // file->Close();
}
