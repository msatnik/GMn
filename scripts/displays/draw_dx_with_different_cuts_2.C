#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>

void draw_dx_with_different_cuts_2() {// main
  
  gStyle->SetOptFit(11111);
  gStyle->SetCanvasPreferGL(1);
  gStyle->SetEndErrorSize(0);

   //const char* filename="/lustre24/expphy/volatile/halla/sbs/msatnik/output/sbs8_70p_simc_deep_sept3.root";
  
  //const char* filename = "/lustre24/expphy/volatile/halla/sbs/msatnik/output/sbs8_50p_simc_deep_sept3.root";

  const char* filename="/lustre24/expphy/volatile/halla/sbs/msatnik/output/sbs9_70p_simc_deep_sept12_z.root";
  
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
  
  // Setting up cuts 
  std::string EnergyCutString = "bb_ps_e>0.2"; // PS and Shower Energy cut
  std::string TrackQualityCutString = "bb_tr_n==1&&bb_gem_track_nhits>=3"; //  coudl also put chi2 cut
  std::string TargetVertexCutString = "abs(bb_tr_vz)<0.06";
  std::string W2CutString = "W2>0.66&&W2<1.10"; 
  std::string FidXCutString = "hcal_x_exp>-2.28&&hcal_x_exp<0.78";
  std::string FidYCutString = "hcal_y_exp>-0.5&&hcal_y_exp<0.5";
  std::string dyCutString = "abs(hcal_dy)<0.3";
  std::string e_over_p_CutString = "abs(e_over_p - 0.98) < 1";// prov says maybe don't use in sim. lose a lot of sym 
  std::string HCal_Energy_CutString = "hcal_e>0.04";
  std::string Optics_CutString = "abs(bb_tr_r_x-bb_tr_r_th*0.9)<0.3 && abs(bb_tr_r_y-0.9*bb_tr_r_ph+0.005)<0.1";

  //std::string W2CutString = "W2>0.66&&W2<1.10"; 
  //std::string FidXCutString = "hcal_x_exp>-1.5&&hcal_x_exp<0.63";
  //std::string FidYCutString = "hcal_y_exp>-0.4&&hcal_y_exp<0.35";
  // std::string dyCutString = "abs(hcal_dy)<0.3";
  
  std::string H_arm_string = "("+FidXCutString +"&&"+ FidYCutString+"&&"+dyCutString+"&&" +HCal_Energy_CutString+")";
  std::string e_arm_string = "("+TrackQualityCutString+"&&"+TargetVertexCutString+"&&"+EnergyCutString+"&&"+W2CutString+"&&"+e_over_p_CutString+"&&"+Optics_CutString+")";
  //std::string H_arm_anti_y_exp_string = "hcal_x_exp>-1.27&& hcal_x_exp<0.63 && abs(hcal_dy) < 0.3 && hcal_e >0.04 && !("+FidYCutString+")";
  std::string both_string = e_arm_string + " && " + H_arm_string;

  // Fill the histograms based on your cuts
  P->Draw("hcal_dx>>orig", "corrected_weight");
  P->Draw(("hcal_dx>>H_arm"), Form("corrected_weight*(%s)", H_arm_string.c_str()));
  P->Draw(("hcal_dx>>H_arm_anticut"), Form("corrected_weight*(!%s)", H_arm_string.c_str()));
  P->Draw(("hcal_dx>>e_arm"), Form("corrected_weight*(%s)", e_arm_string.c_str()));
  P->Draw(("hcal_dx>>e_arm_anticut"), Form("corrected_weight*(!%s)", e_arm_string.c_str()));
  P->Draw(("hcal_dx>>both"), Form("corrected_weight*(%s)", both_string.c_str()));
  P->Draw(("hcal_dx>>both_anticut"), Form("corrected_weight*(!(%s))", both_string.c_str()));
  // P->Draw(("hcal_dx>>H_arm_anti_y_exp"), Form("corrected_weight*(%s)", H_arm_anti_y_exp_string.c_str()));

  // Normalize histograms based on total number of entries
  orig->Scale(1.0 / orig->Integral());
  H_arm->Scale(1.0 / H_arm->Integral());
  H_arm_anticut->Scale(1.0 / H_arm_anticut->Integral());
  e_arm->Scale(1.0 / e_arm->Integral());
  e_arm_anticut->Scale(1.0 / e_arm_anticut->Integral());
  both->Scale(1.0 / both->Integral());
  both_anticut->Scale(1.0 / both_anticut->Integral());
  //H_arm_anti_y_exp->Scale(1.0 / H_arm_anti_y_exp->Integral());

  // Set colors for the histograms
  orig->SetLineWidth(2);
  //orig->SetFillColorAlpha(kBlue, 0.1);
  e_arm->SetLineColor(kMagenta);
  e_arm->SetLineWidth(2);
  //e_arm->SetFillColorAlpha(kMagenta, 0.1);
  e_arm_anticut->SetLineColor(kMagenta);
  e_arm_anticut->SetLineWidth(2);
  //e_arm_anticut->SetFillColorAlpha(kMagenta, 0.1);
  H_arm->SetLineColor(kGreen+2);
  H_arm->SetLineWidth(2);
  //H_arm->SetFillColorAlpha(kGreen+2, 0.1);
  H_arm_anticut->SetLineColor(kGreen+2);
  H_arm_anticut->SetLineWidth(2);
  //H_arm_anticut->SetFillColorAlpha(kGreen+2, 0.1);
  both->SetLineColor(kRed);
  both->SetLineWidth(2);
  // both->SetFillColorAlpha(kRed, 0.1);
  both_anticut->SetLineColor(kRed);
  both_anticut->SetLineWidth(2);
  //both_anticut->SetFillColorAlpha(kRed, 0.1);
  // H_arm_anti_y_exp->SetLineColor(kOrange);
  // H_arm_anti_y_exp->SetLineWidth(2);
  // H_arm_anti_y_exp->SetFillColorAlpha(kOrange, 0.1);

  // Draw histograms on a canvas
  TCanvas *c1 = new TCanvas("c1", "Comparison of Cuts", 800, 600);
  c1->SetGrid();
  both->Draw("E ");
  orig->Draw("E same");             // Draw the original histogram
  H_arm->Draw("E same");
  e_arm->Draw("E same");
  //H_arm_anti_y_exp->Draw("E same hist");

  // Add a legend
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(orig, "Original", "l f");
  legend->AddEntry(e_arm, "e Arm", "l f");
  legend->AddEntry(H_arm, "H Arm", "l");
  legend->AddEntry(both, "Both", "l f");
  //legend->AddEntry(H_arm_anti_y_exp, "Anti-y_exp", "L f");
  legend->Draw();
  c1->Update();

  // Save the canvas to a file
  // c1->SaveAs("hcal_dx_comparison.png");

  // Draw histograms on a canvas
  TCanvas *anticutcanvas = new TCanvas("anticutcanvas", "Comparison of anticuts", 800, 600);
  anticutcanvas->SetGrid();
  both_anticut->Draw("E ");
  orig->Draw("E same");             // Draw the original histogram
  H_arm_anticut->Draw("E same ");
  e_arm_anticut->Draw("E same");
  //H_arm_anti_y_exp->Draw("E same hist");
  anticutcanvas->Update();

   // Add a legend
  TLegend *leged_anticut = new TLegend(0.7, 0.7, 0.9, 0.9);
  leged_anticut->AddEntry(orig, "Original", "l f");
  leged_anticut->AddEntry(e_arm_anticut, "e Arm anticut", "l f");
  leged_anticut->AddEntry(H_arm_anticut, "H Arm anticut", "l");
  leged_anticut->AddEntry(both_anticut, "Both anticut", "l f");
  //legend->AddEntry(H_arm_anti_y_exp, "Anti-y_exp", "L f");
  leged_anticut->Draw();
  anticutcanvas->Update();
  

  // Create ratio plots
  TH1D *ratio_Harm = (TH1D*)H_arm->Clone("ratio_Harm");
  ratio_Harm->Divide(orig);
  
  TH1D *ratio_e_arm = (TH1D*)e_arm->Clone("ratio_e_arm");
  ratio_e_arm->Divide(orig);
  
  TH1D *ratio_both = (TH1D*)both->Clone("ratio_both");
  ratio_both->Divide(orig);

  // Draw the ratio plot
  TCanvas *c2 = new TCanvas("c2", "Ratio Plot", 800, 600);
  ratio_Harm->SetLineColor(kGreen+2);
  ratio_Harm->SetTitle("Ratio to Original;hcal_dx;Ratio");
  ratio_Harm->Draw("E");
  
  ratio_e_arm->SetLineColor(kMagenta);
  ratio_e_arm->Draw("E same");
  
  ratio_both->SetLineColor(kRed);
  ratio_both->Draw("E same");

  // Add legend for the ratio plot
  TLegend *legend2 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend2->AddEntry(ratio_Harm, "H Arm / Orig", "l");
  legend2->AddEntry(ratio_e_arm, "e Arm / Orig", "l");
  legend2->AddEntry(ratio_both, "Both / Orig", "l");
  legend2->Draw();

   // Draw histograms on a canvas
  TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
  c3->SetGrid();
  both->Draw("E");
}//end main
