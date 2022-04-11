#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooUniform.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistFunc.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooRealSumPdf.h"
#include "RooParamHistFunc.h"
#include "RooHistConstraint.h"
#include "RooBinSamplingPdf.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include <iostream>
#include <memory>
using namespace RooFit;
#define _USE_MATH_DEFINES
#include <fstream>
#include <sstream>
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TStyle.h>
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TROOT.h"
#include <TText.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TLegend.h>
using namespace std;

const int eres_n_bin = 53;
const float eres_bin_min = 7;
const float eres_bin_max = 20;
const float eres_bin_width = (eres_bin_max-eres_bin_min)/(eres_n_bin-1);

const int gain_n_bin = 221;
const float gain_bin_min = 10000;
const float gain_bin_max = 65000;
const float gain_bin_width = (gain_bin_max-gain_bin_min)/(gain_n_bin-1);

const int charge_n_bin = 1024;
const float charge_bin_min = 0e-05;
const float charge_bin_max = 200000;

TH2F* charge_spectre = NULL;
TH2F* dead_charge_spectre = NULL;

TH3D *MC_MW8[3] = {NULL, NULL, NULL};
TH3D *MC_MW5[3] = {NULL, NULL, NULL};
TH3D *MC_XW[3]  = {NULL, NULL, NULL};
TH3D *MC_GV[3]  = {NULL, NULL, NULL};

void Load_spectre(){
  TFile *file = new TFile("histo_kolmo/histo_donee/histo_charge_amplitude_energie_611.root", "READ");
  gROOT->cd();
  charge_spectre = (TH2F*)file->Get("histo_pm_charge");
  // charge_spectre->RebinY(4);
  return;
}

void test(){

  for(float gain = gain_bin_min; gain<= gain_bin_max; gain+= gain_bin_width) {
    std::cout << "gain = " << gain  << " bin = " << gain/gain_bin_width - 40 <<'\n';
  }
  // for (float eres = eres_bin_min; eres<= eres_bin_max; eres+= eres_bin_width) {
  //   std::cout << "eres = " << eres << " bin = " << eres/eres_bin_width - 27 <<'\n';
  // }
}

TH1D* spectre_charge_full(int om_number){
  TH1D* spectre_charge = charge_spectre->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
  return spectre_charge;
  // TFile *file = new TFile("Histo_simu_new/Histo_mystere.root", "READ");
  // gROOT->cd();
  // TH1D* spectre_charge = (TH1D*)file->Get("MC");
  // return spectre_charge;
}

void Load_dead_spectre(){
  TFile *file = new TFile("histo_kolmo/histo_donee/histo_charge_amplitude_energie_610.root", "READ");
  gROOT->cd();
  dead_charge_spectre = (TH2F*)file->Get("histo_pm_charge");
  return;
}

TH1D* dead_spectre(int om_number){
  TH1D* dead_charge = dead_charge_spectre->ProjectionY(Form("dead_charge%03d",om_number), om_number+1, om_number+1);
  return dead_charge;
}

void charge_to_ampl() {

  double amplitude, charge;
  int om_number, om;
  double coef;

  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om);
  Result_tree.Branch("coef", &coef);

  TFile *file = new TFile("histo_kolmo/histo_donee/histo_Li_system_611.root", "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("amplitude_tree",1);
  tree->SetBranchAddress("amplitude_tree", &amplitude);
  tree->SetBranchStatus("charge_tree",1);
  tree->SetBranchAddress("charge_tree", &charge);

  for (int i = 0; i < 712; i++) {
    if (i%100 == 0) {
      std::cout << "entry = " << i << '\n';
    }
    TH1D *amplitude_charge = new TH1D ("amplitude_charge", "", 1000, 0, 0.03);
    tree->Project("amplitude_charge", "amplitude_tree/charge_tree", Form("om_number == %i && amplitude_tree > 100", i));
    amplitude_charge->Draw();
    coef = amplitude_charge->GetMean();
    om = i;
    Result_tree.Fill();
    delete amplitude_charge;
  }

  TFile *newfile = new TFile("charge_to_ampl.root", "RECREATE");
  newfile->cd();
  Result_tree.Write();
  newfile->Close();

}

void LoadMC() {
  TFile *histo_file_Tl_MW8 = new TFile("Histo_simu_new/MC_simu_Tl_208_10G_eres_53_gain_221_OC_MW8.root", "READ");
  histo_file_Tl_MW8->cd();
  MC_MW8[0] = (TH3D*)histo_file_Tl_MW8->Get("MC_MW8_Simu_Tl_208_10G_ubc");
  // MC_MW8[0]->RebinZ(4);
  TFile *histo_file_Bi_MW8 = new TFile("Histo_simu_new/MC_simu_Bi_214_10G_eres_53_gain_221_OC_MW8.root", "READ");
  histo_file_Bi_MW8->cd();
  MC_MW8[1] = (TH3D*)histo_file_Bi_MW8->Get("MC_MW8_Simu_Bi_214_10G_ubc");
  // MC_MW8[1]->RebinZ(4);
  TFile *histo_file_K_MW8 = new TFile("Histo_simu_new/MC_simu_K_40_10G_eres_53_gain_221_OC_MW8.root", "READ");
  histo_file_K_MW8->cd();
  MC_MW8[2] = (TH3D*)histo_file_K_MW8->Get("MC_MW8_Simu_K_40_10G_ubc");
  // MC_MW8[2]->RebinZ(4);

  TFile *histo_file_Tl_MW5 = new TFile("Histo_simu_new/MC_simu_Tl_208_10G_eres_53_gain_221_OC_MW5.root", "READ");
  histo_file_Tl_MW5->cd();
  MC_MW5[0] = (TH3D*)histo_file_Tl_MW5->Get("MC__MW5Simu_Tl_208_10G_ubc");
  // MC_MW5[0]->RebinZ(4);
  TFile *histo_file_Bi_MW5 = new TFile("Histo_simu_new/MC_simu_Bi_214_10G_eres_53_gain_221_OC_MW5.root", "READ");
  histo_file_Bi_MW5->cd();
  MC_MW5[1] = (TH3D*)histo_file_Bi_MW5->Get("MC__MW5Simu_Bi_214_10G_ubc");
  // MC_MW5[1]->RebinZ(4);
  TFile *histo_file_K_MW5 = new TFile("Histo_simu_new/MC_simu_K_40_10G_eres_53_gain_221_OC_MW5.root", "READ");
  histo_file_K_MW5->cd();
  MC_MW5[2] = (TH3D*)histo_file_K_MW5->Get("MC__MW5Simu_K_40_10G_ubc");
  // MC_MW5[2]->RebinZ(4);

  TFile *histo_file_TL_XW = new TFile("Histo_simu_new/MC_simu_Tl_208_10G_eres_53_gain_221_OC_XW.root", "READ");
  histo_file_TL_XW->cd();
  MC_XW[0] = (TH3D*)histo_file_TL_XW->Get("MC_XW_Simu_Tl_208_10G_ubc");
  // MC_XW[0]->RebinZ(4);
  TFile *histo_file_Bi_XW = new TFile("Histo_simu_new/MC_simu_Bi_214_10G_eres_53_gain_221_OC_XW.root", "READ");
  histo_file_Bi_XW->cd();
  MC_XW[1] = (TH3D*)histo_file_Bi_XW->Get("MC_XW_Simu_Bi_214_10G_ubc");
  // MC_XW[1]->RebinZ(4);
  TFile *histo_file_K_XW = new TFile("Histo_simu_new/MC_simu_K_40_10G_eres_53_gain_221_OC_XW.root", "READ");
  histo_file_K_XW->cd();
  MC_XW[2] = (TH3D*)histo_file_K_XW->Get("MC_XW_Simu_K_40_10G_ubc");
  // MC_XW[2]->RebinZ(4);

  TFile *histo_file_Tl_GV = new TFile("Histo_simu_new/MC_simu_Tl_208_10G_eres_53_gain_221_OC_GV.root", "READ");
  histo_file_Tl_GV->cd();
  MC_GV[0] = (TH3D*)histo_file_Tl_GV->Get("MC_GV_Simu_Tl_208_10G_ubc");
  // MC_GV[0]->RebinZ(4);
  TFile *histo_file_Bi_GV = new TFile("Histo_simu_new/MC_simu_Bi_214_10G_eres_53_gain_221_OC_GV.root", "READ");
  histo_file_Bi_GV->cd();
  MC_GV[1] = (TH3D*)histo_file_Bi_GV->Get("MC_GV_Simu_Bi_214_10G_ubc");
  // MC_GV[1]->RebinZ(4);
  TFile *histo_file_K_GV = new TFile("Histo_simu_new/MC_simu_K_40_10G_eres_53_gain_221_OC_GV.root", "READ");
  histo_file_K_GV->cd();
  MC_GV[2] = (TH3D*)histo_file_K_GV->Get("MC_GV_Simu_K_40_10G_ubc");
  // MC_GV[2]->RebinZ(4);
}

float dead_time(int om, float gain){
  TH1D* full_spectre = NULL;
  full_spectre = spectre_charge_full(om);
  full_spectre->Scale(1.0/(15*60+2));
  TH1D* deadspectre = NULL;
  deadspectre = dead_spectre(om);
  deadspectre->Scale(1.0/(30*60+1));
  float dead_time = 0;
  for (int i = 0; i < gain*2.2; i++) {
    full_spectre->SetBinContent(i,0);
    deadspectre->SetBinContent(i,0);
  }
  TCanvas* canvas = new TCanvas;
  canvas->SetLogy();
  full_spectre->Draw();
  deadspectre->Draw("same");
  dead_time = 1/(full_spectre->Integral()/deadspectre->Integral());
  return dead_time;
  delete full_spectre;
  delete deadspectre;
  return dead_time;
  delete canvas;
}

TH3D* MC_chooser(int om, int comp){
  if (om <520 && om%13 != 12 && om%13 != 0){
    return MC_MW8[comp];
  }
  else if (om <520 && (om%13 == 12 || om%13 == 0)) {
    return MC_MW5[comp];
  }
  else if (om > 519 && om < 648) {
    return MC_XW[comp];
  }
  else if (om > 647) {
    return MC_GV[comp];
  }

  return NULL;
}

int eres_chooser(int om){
  int eres_chooser = 0;
  if (om <520 && om%13 != 12 && om%13 != 0){
    eres_chooser = 8;
  }
  else if (om <520 && (om%13 == 12 || om%13 == 0)) {
    eres_chooser = 9;
  }
  else if (om > 519 && om < 648) {
    eres_chooser = 12;
  }
  else if (om > 647) {
    eres_chooser = 15;
  }
  return eres_chooser;
}

double* om_gain_fit(int om){
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);

  double mean = 0;
  double sigma = 0;
  double n_evt = 0;
  double nbg = 0;
  double chi = 0;
  double ndf = 0;
  double chin = 0;

  double* tab = new double[7];
  TCanvas* canvas = new TCanvas;
  canvas->SetLogy();

  TH1D* spectre_om = NULL;
  spectre_om = spectre_charge_full(om);
  TF1* f_ComptonEdgePoly = new TF1 ("f_ComptonEdgePoly","[0]/2.0*(1+TMath::Erf(([1]-x)/(TMath::Sqrt(2)*[2])))+[3]*x", 40000, 90000);
  f_ComptonEdgePoly->SetParNames("N_evt","Mean","Sigma","Nbg" );

  if ((om % 13) == 12 )        //om multiple de (13)-1
  {
    f_ComptonEdgePoly->SetParameters(120, 72367, 6893, 3.91e-5);
    f_ComptonEdgePoly->SetRange(6000,100000);
    f_ComptonEdgePoly->Draw("same");
    spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
  }
  else if ((om % 13) == 0)       //om multiple de 13
  {
    f_ComptonEdgePoly->SetParameters(112, 68168, 5604, 1.2e-05);
    f_ComptonEdgePoly->SetRange(50000,100000);
    f_ComptonEdgePoly->Draw("same");
    spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-1.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
  }
  else         //om normaux (8pouces)
  {
    f_ComptonEdgePoly->SetParameters(111, 60978, 3787, 4.19e-05);
    f_ComptonEdgePoly->SetRange(55000,100000);
    f_ComptonEdgePoly->Draw("same");
    spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+7.5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-3.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+7.5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
  }
  chi = (f_ComptonEdgePoly->GetChisquare());
  std::cout << chi << '\n';
  ndf = (f_ComptonEdgePoly->GetNDF());
  chin = (f_ComptonEdgePoly->GetChisquare()/f_ComptonEdgePoly->GetNDF());

  if (chin < 1.5 && mean > 7000) {
    n_evt = f_ComptonEdgePoly->GetParameter(0);
    mean = (f_ComptonEdgePoly->GetParameter(1))/2.6;
    sigma = (f_ComptonEdgePoly->GetParameter(2))/2.6;
    nbg = f_ComptonEdgePoly->GetParameter(3);
  }
  else{
    n_evt = 0;
    mean = 0;
    sigma = 0;
    nbg = 0;
  }
  delete f_ComptonEdgePoly;
  delete canvas;
  tab[0] = mean;
  tab[1] = sigma;
  tab[2] = chin;
  tab[3] = chi;
  tab[4] = ndf;
  tab[5] = n_evt;
  tab[6] = nbg;
  return tab;
}

int* om_chooser(string om){
  int* lim = new int[2];
  if(om == "-MWIT" || om == "-1"){
    lim[0] = 0;
    lim[1] = 260;
  }
  if(om == "-MWFR" || om == "-2"){
    lim[0] = 260;
    lim[1] = 520;
  }
  if(om == "-XW" || om == "-3"){
    lim[0] = 520;
    lim[1] = 648;
  }
  if(om == "-GV" || om == "-4"){
    lim[0] = 648;
    lim[1] = 712;
  }
  if(om == "-FULL"  ){
    lim[0] = 0;
    lim[1] = 712;
  }
  return lim;
}

double events_counter(double &mean, float min, TH1D* spectre_om){
  double frac = spectre_om->Integral((mean)*1024/200000, 1024)/spectre_om->Integral();
  std::cout << "mean = " << (mean)*1024/200000 << '\n';
  while (frac < 0.3) {
    mean-= 250;
    frac = spectre_om->Integral((mean)*1024/200000, 1024)/spectre_om->Integral();
  }
  return mean;
}

double* roofitter(int om, double gain, double eres, TH1D* spectre_om, TH1D* mc0, TH1D* mc1, TH1D* mc2, double *rootab, int bin){
  const float start_x = spectre_om->GetXaxis()->GetBinUpEdge(bin);
  int nbin = 0;
  for (int i = bin; i < spectre_om->GetNbinsX(); i++) {
    if (spectre_om->GetBinContent(i) != 0) {
      nbin++;
    }
  }
  std::cout << "nbin = " << nbin << '\n';

  RooRealVar x("x", "x", 10000, 130000);
  x.setBins(615);
  RooDataHist Tl("Tl", "Tl", x, Import(*mc0));
  RooDataHist Bi("Bi", "Bi", x, Import(*mc1));
  RooDataHist K("K", "K", x, Import(*mc2));

  RooRealVar RooTl("RooTl", "Tl", 0.2, 0.1, 0.4);
  RooRealVar RooBi("RooBi", "Bi", 0.5, 0.3, 0.6);
  RooRealVar RooK("RooK", "K", 0.3, 0.2, 0.4);

  RooHistPdf Tl_pdf ("Tl_pdf", "", x, Tl);
  RooHistPdf Bi_pdf ("Bi_pdf", "", x, Bi);
  RooHistPdf K_pdf ("K_pdf", "", x, K);
  RooBinSamplingPdf Tl_spdf ("Tl_spdf", "", x, Tl_pdf);
  RooBinSamplingPdf Bi_spdf ("Bi_spdf", "", x, Bi_pdf);
  RooBinSamplingPdf K_spdf ("K_spdf", "", x, K_pdf);

  RooDataHist spectre_data("spectre_data", "spectre_data", x, Import(*spectre_om));

  RooParamHistFunc p_ph_Tl("p_ph_Tl","p_ph_Tl",Tl);
  RooParamHistFunc p_ph_Bi("p_ph_Bi","p_ph_Bi",Bi);
  RooParamHistFunc p_ph_K("p_ph_K","p_ph_K",K);

  RooAddPdf sum_simu("sum_simu", "sum_simu",
                      RooArgList(Tl_spdf, Bi_spdf, K_spdf),
                      RooArgList(RooTl, RooBi));



  RooNLLVar RooLogL("LogL", "LogL", sum_simu, spectre_data);
  RooChi2Var RooChi2("Chi2", "Chi2", sum_simu, spectre_data);   // create the variance
  std::cout << "test1 = " << RooChi2.getVal()/(nbin - 2) << '\n';

  // RooMinimizer *miniLog = new RooMinimizer(RooLogL);
  // double ConvLog = miniLog->minimize("Minuit", "Migrad");
  //
  //   std::cout << "test2 = " << RooChi2.getVal()/(nbin - 2) << '\n';

  RooMinimizer *miniChi = new RooMinimizer(RooChi2);


  double ConvChi = miniChi->minimize("Minuit", "Migrad");  //   Create the Chi2/LogL

  double Chi2 = RooChi2.getVal()/(nbin - 2);
  std::cout << "test3 = " << Chi2 << '\n';

  TCanvas* can = new TCanvas;
  can->cd();
  auto frame = x.frame(Title("Fit gain simu"));

  spectre_data.plotOn(frame, DataError(RooAbsData::SumW2), DrawOption("P"));

  // TH1* test3 = test.createHistogram("", RooTl, Binning(100,0,1));
  // test3->Draw();

  sum_simu.plotOn(frame, FillColor(0));
  spectre_data.plotOn(frame, DataError(RooAbsData::SumW2), DrawOption("P"));

  // Plot model components
  sum_simu.plotOn(frame, Components(Tl_spdf), LineColor(kGreen), Name("Tl_curve"));
  sum_simu.plotOn(frame, Components(Bi_spdf), LineColor(kYellow), Name("Bi_curve"));
  sum_simu.plotOn(frame, Components(K_spdf), LineColor(kBlue), Name("K_curve"));

  sum_simu.plotOn(frame, LineColor(kRed), Name("sum_curve"), Range(start_x, 130000));

  double Tl_int = RooTl.getVal();
  double Bi_int = RooBi.getVal();
  double K_int = 1- Tl_int -Bi_int;
  double int_tot = Tl_int + Bi_int + K_int;
  double Tl_int_error = RooTl.getError();
  double Bi_int_error = RooBi.getError();
  double K_int_error = sqrt(pow((Bi_int_error), 2) + pow((Tl_int_error), 2));
  std::cout << "RooTl = " << Tl_int << " +- " << Tl_int_error << '\n';
  std::cout << "RooBi = " << Bi_int << " +- " << Bi_int_error << '\n';
  std::cout << "RooK = " << K_int << " +- " << K_int_error << '\n';
  frame->GetYaxis()->UnZoom();
  frame->GetYaxis()->SetTitle("n events");
  frame->GetXaxis()->SetTitle("charge (u.a)");
  frame ->Draw();
  can->SetLogy();
  // can->SaveAs("test.png");

  rootab[0] = RooLogL.getVal();
  rootab[1] = Tl_int;
  rootab[2] = Bi_int;
  rootab[3] = K_int;
  rootab[5] = Tl_int_error;
  rootab[6] = Bi_int_error;
  rootab[7] = K_int_error;
  rootab[4] = RooChi2.getVal()/(nbin - 2);
  TLatex l;
  l.SetTextFont(40);
  l.DrawLatex(90000, 80, Form("Khi2/NDF = %.2f",rootab[4]));
  can->SaveAs(Form("OM_fit/om_%d/best_fit_om_%d_eres_%.2f_gain_%.0f_bin_%d.png", om, om, eres, gain, bin));

  delete miniChi;
  delete can;
  delete frame;
  // delete miniLog;
  return rootab;
}

double* get_om_eff(string compo, int om, double *om_eff) {
  om_eff[0] = 0;
  om_eff[1] = 0;
  TFile *newfile = new TFile(Form("eff_om/eff_om_%s.root", compo.c_str()), "READ");

  TTree* tree = (TTree*)newfile->Get("new_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("eff_tot",1);
  tree->SetBranchAddress("eff_tot", &om_eff[0]);
  tree->SetBranchStatus("eff_tot_error",1);
  tree->SetBranchAddress("eff_tot_error", &om_eff[1]);

  tree->GetEntry(om);
  newfile->Close();
  return om_eff;
}

void Fit_Gain_Simu(string wall = "", int p = 0; std::vector<int> *om_vector = nullptr, std::vector<int> *param_value = nullptr) {
  LoadMC();
  Load_spectre();
  Load_dead_spectre();
  TH1::SetDefaultSumw2();

  double logL, Chi2, param1, param2, param3, om_counting_Tl, om_counting_Bi, om_counting_K, eff_Tl, eff_Bi, eff_K, eff_error_Tl, eff_error_Bi, eff_error_K;
  double int_cut_mc0, int_cut_mc1, int_cut_mc2, int_full_mc0, int_full_mc1, int_full_mc2;
  float gain, eres = 0;
  int om_number = 0;
  double mean_erf = 0;
  double om_flux_Tl, om_flux_Bi, om_flux_K;
  int lim =0;
  float deadt;
  TFile *newfile = new TFile(Form("histo_fit/fit_test%s.root", wall.c_str()), "RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("logL", &logL);
  Result_tree.Branch("Chi2", &Chi2);
  Result_tree.Branch("param1", &param1);
  Result_tree.Branch("param2", &param2);
  Result_tree.Branch("param3", &param3);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("eres", &eres);
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("om_flux_Tl", &om_flux_Tl);
  Result_tree.Branch("om_flux_Bi", &om_flux_Bi);
  Result_tree.Branch("om_flux_K", &om_flux_K);
  Result_tree.Branch("lim", &lim);
  Result_tree.Branch("deadt", &deadt);
  Result_tree.Branch("om_counting_Tl", &om_counting_Tl);
  Result_tree.Branch("om_counting_Bi", &om_counting_Bi);
  Result_tree.Branch("om_counting_K", &om_counting_K);
  Result_tree.Branch("eff_Tl", &eff_Tl);
  Result_tree.Branch("eff_Bi", &eff_Bi);
  Result_tree.Branch("eff_K", &eff_K);

  std::ifstream  charge("gain_data/Resultat_mean_energie_charge_total.txt");
  float charge_valeur_fit [712];
  memset(charge_valeur_fit, -1, 712*sizeof(float));
  int charge_om_num;
  while (charge >> charge_om_num)
  {
    charge >> charge_valeur_fit[charge_om_num];
  }
  double* rootab = new double[8];
  double* eff = new double[2];

  float min = 0.85;
  float max = 1.15;

  int* borne = new int[2];
  // borne = om_chooser(wall);

  for (int om = 285; om < 286; om++)
  for (int om = borne[0]; om < borne[1]; om = om +1)
  {
    TH3D* MC_Tl_208 = MC_chooser(om, 0);
    TH3D* MC_Bi_214 = MC_chooser(om, 1);
    TH3D* MC_K_40 = MC_chooser(om, 2);
    TH1D* spectre_om = NULL;
    spectre_om = spectre_charge_full(om);

    if (om < 520) {
      if ((spectre_om->GetEntries() < 100) || (spectre_om->GetMean(1) < 1500)){
        std::cout << "om " << om << " deleted"<< '\n';
        delete spectre_om;
        om++;
        spectre_om = spectre_charge_full(om);
        get_om_eff("Tl_208", om, eff);
        eff_Tl = eff[0];
        eff_error_Tl = eff[1];
        get_om_eff("Tl_208", om, eff);
        eff_Bi = eff[0];
        eff_error_Bi = eff[1];
        get_om_eff("Tl_208", om, eff);
        eff_K = eff[0];
        eff_error_K = eff[1];
      }
      if (charge_valeur_fit[om] != -1){
        double *tab = om_gain_fit(om);
        std::cout << tab[0] << '\n';
        if (tab[0] == 0 ) {
          mean_erf = 1/charge_valeur_fit[om];
        }
        else{
          mean_erf = tab[0];
        }
        delete tab;
      }
      else{
        mean_erf = 19000;
        min = 0.7;
        max = 1.2;
      }
    }
    else if (om < 648 && om > 519) {
      mean_erf = 1/charge_valeur_fit[om];
    }
    else{
      mean_erf = 19000;
      min = 0.8;
      max = 1.3;
    }
    std::cout << "mean erf = " << mean_erf<< '\n';
    if (spectre_om->Integral((mean_erf)*1024/200000, 1024)/spectre_om->Integral() < 0.3) {
      mean_erf = events_counter(mean_erf, min, spectre_om);
      std::cout << "mean_erf = " << mean_erf << '\n';
    }
    // if (spectre_om->Integral(mean_erf, 1024) < 400) {
    //   continue;
    // }


    om_number = om;
    int eres = eres_chooser(om);
    int eres_count = (eres-eres_bin_min)/eres_bin_width+1;
    for (int gain_count = 0; gain_count <150; gain_count++) {
      // std::cout << (gain_bin_min + gain_bin_width*(gain_count-1))<< "   et    lim_inf = " << 1/(mean_erf*min) << "   sup  ="  << 1/(mean_erf*max) << '\n';
      if (((gain_bin_min + gain_bin_width*(gain_count-1)) > mean_erf*min) && ((gain_bin_min + gain_bin_width*(gain_count-1))<mean_erf*max)){
        spectre_om = spectre_charge_full(om);

        gain = (gain_bin_min + gain_bin_width*(gain_count-1));
        std::cout << "gain_count = " << gain_count << " and gain = " << gain << '\n';
        // gain_count = 54;
        // eres_count = 5;
        TH1D *mc0 = MC_Tl_208->ProjectionZ("Charge_Tl_208", eres_count, eres_count, gain_count, gain_count);    // first MC histogram
        TH1D *mc1 = MC_Bi_214->ProjectionZ("Charge_Bi_214", eres_count, eres_count, gain_count, gain_count);    // second MC histogram
        TH1D *mc2 = MC_K_40->ProjectionZ("Charge_K_40", eres_count, eres_count, gain_count, gain_count);    // second MC histogram
        int bin = 1.1*gain*1024/200000.0 + 1;

        if (om == 570) {
          delete spectre_om;
          TFile *tfile = new TFile("histo_kolmo/histo_donee/histo_charge_amplitude_energie_611_570.root", "READ");
          gROOT->cd();
          TH2F* charge_spectre_570 = (TH2F*)tfile->Get("histo_pm_charge");
          TH1D* spectre_om = charge_spectre_570->ProjectionY("charge570", 571, 571);
        }

        int_full_mc0 = mc0->Integral();
        int_full_mc1 = mc2->Integral();
        int_full_mc2 = mc2->Integral();

        for (int i = 0; i < bin; i++) {
          mc0->SetBinContent(i, 0);
          mc1->SetBinContent(i, 0);
          mc2->SetBinContent(i, 0);
          spectre_om->SetBinContent(i, 0);
          spectre_om->SetBinError(i, 0);
        }

        for (int i = bin*4/1.1; i < 1024; i++) {
          spectre_om->SetBinContent(i, 0);
          spectre_om->SetBinError(i, 0);
        }

        int_cut_mc0 = mc0->Integral();
        int_cut_mc1 = mc1->Integral();
        int_cut_mc2 = mc2->Integral();

        roofitter(om, gain, eres, spectre_om, mc0, mc1, mc2, rootab, bin);
        logL = rootab[0];
        Chi2 = rootab[4];
        param1 = rootab[1];
        param2 = rootab[2];
        param3 = rootab[3];
        deadt = dead_time(om, bin);

        if (deadt < 1) {
          deadt = 1;
        }
        delete mc0;
        delete mc1;
        delete mc2;
        delete spectre_om;
        mc0 = MC_Tl_208->ProjectionZ("Charge_Tl_208", eres_count, eres_count, gain_count, gain_count);
        mc1 = MC_Bi_214->ProjectionZ("Charge_Bi_214", eres_count, eres_count, gain_count, gain_count);
        mc2 = MC_K_40->ProjectionZ("Charge_K_40", eres_count, eres_count, gain_count, gain_count);
        spectre_om = spectre_charge_full(om);

        mc0->Scale((param1/int_cut_mc0)*spectre_om->Integral());
        om_counting_Tl = mc0->Integral();
        om_flux_Tl = (mc0->Integral()*deadt)/(902*eff_Tl*250000);
        mc1->Scale(spectre_om->Integral()*param2/int_cut_mc1);
        om_counting_Bi = mc1->Integral();
        om_flux_Bi = (mc1->Integral()*deadt)/(902*eff_Bi*250000);
        mc2->Scale(spectre_om->Integral()*param3/int_cut_mc2);
        om_counting_K = mc2->Integral();
        om_flux_K = (mc2->Integral()*deadt)/(902*eff_K*250000);

        Result_tree.Fill();
        delete mc0;
        delete mc1;
        delete mc2;
        delete spectre_om;
      }
    }
  }
  newfile->cd();
  Result_tree.Write();
  newfile->Close();
  std::cout << "good end" << '\n';
}

void Fit_flux() {
  LoadMC();
  Load_spectre();
  Load_dead_spectre();
  TH1::SetDefaultSumw2();

  std::ofstream outFile("gain_data/Mean_value.txt");

  double logL, Chi2, param1, param2, param3, om_counting_Tl, om_counting_Bi, om_counting_K, eff_Tl, eff_Bi, eff_K, eff_error_Tl, eff_error_Bi, eff_error_K, om_counting_Tl_error, om_counting_Bi_error, om_counting_K_error;
  double int_cut_mc0, int_cut_mc1, int_cut_mc2, int_full_mc0, int_full_mc1, int_full_mc2;
  float gain, eres = 0;
  int om_number = 0;
  double mean_erf = 0;
  double om_flux_Tl, om_flux_Bi, om_flux_K;
  int lim =0;
  float deadt;
  double test;

  TFile *gainfile = new TFile("histo_fit/distrib_test_.root", "READ");
  TTree* tree = (TTree*)gainfile->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("gain",1);
  tree->SetBranchAddress("gain", &test);
  double gain_tab[712];
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    gain_tab[i] = test;
  }
  TFile *newfile = new TFile("histo_fit/fit_test3.root", "RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("logL", &logL);
  Result_tree.Branch("Chi2", &Chi2);
  Result_tree.Branch("param1", &param1);
  Result_tree.Branch("param2", &param2);
  Result_tree.Branch("param3", &param3);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("eres", &eres);
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("om_flux_Tl", &om_flux_Tl);
  Result_tree.Branch("om_flux_Bi", &om_flux_Bi);
  Result_tree.Branch("om_flux_K", &om_flux_K);
  Result_tree.Branch("lim", &lim);
  Result_tree.Branch("deadt", &deadt);
  Result_tree.Branch("om_counting_Tl", &om_counting_Tl);
  Result_tree.Branch("om_counting_Bi", &om_counting_Bi);
  Result_tree.Branch("om_counting_K", &om_counting_K);
  Result_tree.Branch("om_counting_Tl_error", &om_counting_Tl_error);
  Result_tree.Branch("om_counting_Bi_error", &om_counting_Bi_error);
  Result_tree.Branch("om_counting_K_error", &om_counting_K_error);
  Result_tree.Branch("eff_Tl", &eff_Tl);
  Result_tree.Branch("eff_Bi", &eff_Bi);
  Result_tree.Branch("eff_K", &eff_K);
  Result_tree.Branch("eff_error_Tl", &eff_error_Tl);
  Result_tree.Branch("eff_error_Bi", &eff_error_Bi);
  Result_tree.Branch("eff_error_K", &eff_error_K);


  double* eff = new double[2];
  double* rootab = new double[8];
  double facteur;
  for (int om = 0; om < 712; om++)
  {
    TH3D* MC_Tl_208 = MC_chooser(om, 0);
    TH3D* MC_Bi_214 = MC_chooser(om, 1);
    TH3D* MC_K_40 = MC_chooser(om, 2);
    TH1D* spectre_om = NULL;
    spectre_om = spectre_charge_full(om);

    get_om_eff("Tl_208", om, eff);
    eff_Tl = eff[0];
    eff_error_Tl = eff[1];

    get_om_eff("Tl_208", om, eff);
    eff_Bi = eff[0];
    eff_error_Bi = eff[1];

    get_om_eff("Tl_208", om, eff);
    eff_K = eff[0];
    eff_error_K = eff[1];

    if ((spectre_om->GetEntries() < 100) || (spectre_om->GetMean(1) < 1500)){
      delete spectre_om;
      gain = 0;
      Chi2 = 0;
      om_flux_Tl = 0;
      om_flux_Bi = 0;
      om_flux_K = 0;
      om_counting_K = 0;
      om_counting_Tl = 0;
      om_counting_Bi = 0;
      deadt = 0;
      eff_Tl = 0;
      eff_Bi = 0;
      eff_K = 0;
      Result_tree.Fill();
      om++;
      spectre_om = spectre_charge_full(om);
      get_om_eff("Tl_208", om, eff);
      eff_Tl = eff[0];
      eff_error_Tl = eff[1];
      get_om_eff("Tl_208", om, eff);
      eff_Bi = eff[0];
      eff_error_Bi = eff[1];
      get_om_eff("Tl_208", om, eff);
      eff_K = eff[0];
      eff_error_K = eff[1];
    }
    om_number = om;
    int eres = eres_chooser(om);
    int eres_count = (eres-eres_bin_min)/eres_bin_width+1;
    spectre_om = spectre_charge_full(om);
    gain = gain_tab[om];
    if (om == 192) {
      gain = 17750;
    }
    if (om == 318) {
      gain = 14500;
    }
    if (om == 493) {
      gain = 29500;
    }
    if (om == 574) {
      gain = 21250;
    }
    if (om == 575) {
      gain = 17500;
    }
    if (om == 593) {
      gain = 24250;
    }
    if (om == 594) {
      gain = 23500;
    }
    if (om == 616) {
      gain = 32000;
    }
    if (om == 662) {
      gain = 12500;
    }
    if (om == 697) {
      gain = 19000;
    }
    if (om == 9) {
      gain = 25000;
    }
    int gain_count =((gain - 10000)/250) + 1;

    TH1D *mc0 = MC_Tl_208->ProjectionZ("Charge_Tl_208", eres_count, eres_count, gain_count, gain_count);    // first MC histogram
    TH1D *mc1 = MC_Bi_214->ProjectionZ("Charge_Bi_214", eres_count, eres_count, gain_count, gain_count);    // second MC histogram
    TH1D *mc2 = MC_K_40->ProjectionZ("Charge_K_40", eres_count, eres_count, gain_count, gain_count);    // second MC histogram
    int bin = 1.1*gain*1024/200000.0/4 + 1;   // rebin 4
    // int bin = 1.1*gain*1024/200000.0 + 1;

    if (om == 570) {
      delete spectre_om;
      TFile *tfile = new TFile("histo_kolmo/histo_donee/histo_charge_amplitude_energie_611_570.root", "READ");
      gROOT->cd();
      TH2F* charge_spectre_570 = (TH2F*)tfile->Get("histo_pm_charge");
      TH1D* spectre_om = charge_spectre_570->ProjectionY("charge570", 571, 571);
    }

    int_full_mc0 = mc0->Integral();
    int_full_mc1 = mc2->Integral();
    int_full_mc2 = mc2->Integral();

    for (int i = 0; i < bin; i++) {
      mc0->SetBinContent(i, 0);
      mc1->SetBinContent(i, 0);
      mc2->SetBinContent(i, 0);
      spectre_om->SetBinContent(i, 0);
      spectre_om->SetBinError(i, 0);
    }

    if (om < 520) {
      for (int i = bin*4/1.1; i < 1024; i++) {
        spectre_om->SetBinContent(i, 0);
        spectre_om->SetBinError(i, 0);
      }
    }

    int_cut_mc0 = mc0->Integral();
    int_cut_mc1 = mc1->Integral();
    int_cut_mc2 = mc2->Integral();
    double sp_om_cut = spectre_om->Integral();
    roofitter(om, gain, eres, spectre_om, mc0, mc1, mc2, rootab, bin);
    // return;
    logL = rootab[0];
    Chi2 = rootab[4];
    param1 = rootab[1];
    param2 = rootab[2];
    param3 = rootab[3];
    double param1error = rootab[5];
    std::cout << "param1 = " << param1error << '\n';
    double param2error = rootab[6];
    double param3error = rootab[7];
    deadt = dead_time(om, bin);
    if (deadt < 1) {
      deadt = 1;
    }
    delete mc0;
    delete mc1;
    delete mc2;
    delete spectre_om;


    mc0 = MC_Tl_208->ProjectionZ("Charge_Tl_208", eres_count, eres_count, gain_count, gain_count);
    mc1 = MC_Bi_214->ProjectionZ("Charge_Bi_214", eres_count, eres_count, gain_count, gain_count);
    mc2 = MC_K_40->ProjectionZ("Charge_K_40", eres_count, eres_count, gain_count, gain_count);
    spectre_om = spectre_charge_full(om);

    for (int i = 0; i < bin; i++) {
      mc0->SetBinContent(i, 0);
      mc1->SetBinContent(i, 0);
      mc2->SetBinContent(i, 0);
      spectre_om->SetBinContent(i, 0);
      spectre_om->SetBinError(i, 0);
    }


    // mc0->Scale(sp_om_cut*(param1/int_cut_mc0));
    // mc1->Scale(sp_om_cut*param2/int_cut_mc1);
    // mc2->Scale(sp_om_cut*param3/int_cut_mc2);

    // for (int i = 0; i < gain*2*1024/200000.; i++) {
    //   mc0->SetBinContent(i, 0);
    // }

    om_counting_Tl = mc0->Integral()*deadt/900*(sp_om_cut*(param1/int_cut_mc0));
    om_counting_Tl_error = mc0->Integral()*deadt/900*(sp_om_cut*(param1error/int_cut_mc0));
    om_counting_Bi = mc1->Integral()*deadt/900*(sp_om_cut*(param2/int_cut_mc1));
    om_counting_Bi_error = mc1->Integral()*deadt/900*(sp_om_cut*(param2error/int_cut_mc1));
    om_counting_K = mc2->Integral()*deadt/900*(sp_om_cut*(param3/int_cut_mc2));
    om_counting_K_error = mc2->Integral()*deadt/900*(sp_om_cut*(param3error/int_cut_mc2));

    outFile << om << "\t" << gain << endl;
    Result_tree.Fill();

    delete mc0;
    delete mc1;
    delete mc2;
    delete spectre_om;
  }
  newfile->cd();
  Result_tree.Write();
  newfile->Close();
  std::cout << "good end" << '\n';
}

void distrib(string name) {
  TFile *eff_file = new TFile(Form("histo_fit/%s.root", name.c_str()), "READ");
  // std::ofstream outFile("Best_Khi2.txt");
  int om_number =0;
  double Chi2, om_counting_Tl, om_counting_Bi, om_counting_K, counting_Tl, counting_Bi, counting_K, effTl, effBi, effK, eff_Tl, eff_Bi, eff_K;
  float gain = 0;
  float eres = 0;
  double om_flux_K = 0;
  double om_flux_Bi = 0;
  double om_flux_Tl = 0;
  double param1 = 0;
  double param2 = 0;
  float deadt = 0;
  TTree* eff_tree = (TTree*)eff_file->Get("Result_tree");
  eff_tree->SetBranchStatus("*",0);
  eff_tree->SetBranchStatus("om_number",1);
  eff_tree->SetBranchAddress("om_number", &om_number);
  eff_tree->SetBranchStatus("Chi2",1);
  eff_tree->SetBranchAddress("Chi2", &Chi2);
  eff_tree->SetBranchStatus("gain",1);
  eff_tree->SetBranchAddress("gain", &gain);
  eff_tree->SetBranchStatus("eres",1);
  eff_tree->SetBranchAddress("eres", &eres);
  eff_tree->SetBranchStatus("om_flux_Tl",1);
  eff_tree->SetBranchAddress("om_flux_Tl", &om_flux_Tl);
  eff_tree->SetBranchStatus("om_flux_Bi",1);
  eff_tree->SetBranchAddress("om_flux_Bi", &om_flux_Bi);
  eff_tree->SetBranchStatus("om_flux_K",1);
  eff_tree->SetBranchAddress("om_flux_K", &om_flux_K);
  eff_tree->SetBranchStatus("param1",1);
  eff_tree->SetBranchAddress("param1", &param1);
  eff_tree->SetBranchStatus("param2",1);
  eff_tree->SetBranchAddress("param2", &param2);
  eff_tree->SetBranchStatus("deadt",1);
  eff_tree->SetBranchAddress("deadt", &deadt);
  eff_tree->SetBranchStatus("om_counting_Tl",1);
  eff_tree->SetBranchAddress("om_counting_Tl", &om_counting_Tl);
  eff_tree->SetBranchStatus("om_counting_Bi",1);
  eff_tree->SetBranchAddress("om_counting_Bi", &om_counting_Bi);
  eff_tree->SetBranchStatus("om_counting_K",1);
  eff_tree->SetBranchAddress("om_counting_K", &om_counting_K);
  eff_tree->SetBranchStatus("eff_Tl",1);
  eff_tree->SetBranchAddress("eff_Tl", &eff_Tl);
  eff_tree->SetBranchStatus("eff_Bi",1);
  eff_tree->SetBranchAddress("eff_Bi", &eff_Bi);
  eff_tree->SetBranchStatus("eff_K",1);
  eff_tree->SetBranchAddress("eff_K", &eff_K);

  double Gain = 0;
  double Tl = 10000;
  double Bi = 10000;
  double K = 10000;
  double par1 = 10000;
  double par2 = 10000;
  int om = 0;
  double test = 10;
  float dt = 100;
  TFile *newfile = new TFile(Form("histo_fit/distrib_%s.root", name.c_str()), "RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om);
  Result_tree.Branch("Chi2", &test);
  Result_tree.Branch("gain", &Gain);
  Result_tree.Branch("om_flux_Tl", &Tl);
  Result_tree.Branch("om_flux_Bi", &Bi);
  Result_tree.Branch("om_flux_K", &K);
  Result_tree.Branch("param1", &par1);
  Result_tree.Branch("param2", &par2);
  Result_tree.Branch("deadt", &dt);
  Result_tree.Branch("om_counting_Tl", &counting_Tl);
  Result_tree.Branch("om_counting_Bi", &counting_Bi);
  Result_tree.Branch("om_counting_K", &counting_K);
  Result_tree.Branch("eff_Tl", &effTl);
  Result_tree.Branch("eff_Bi", &effBi);
  Result_tree.Branch("eff_K", &effK);

  for (int j = 0; j < 712; j++) {
    for (double i = 0; i < eff_tree->GetEntries(); i++) {
      eff_tree->GetEntry(i);
      if (om_number == j) {
        if (Chi2 < test) {
          test = Chi2;
          Gain = gain;
          Tl = om_flux_Tl;
          Bi = om_flux_Bi;
          K = om_flux_K;
          par1 = param1;
          par2 = param2;
          dt = deadt;
          counting_Tl = om_counting_Tl;
          counting_Bi = om_counting_Bi;
          counting_K = om_counting_K;
          effTl = eff_Tl;
          effBi = eff_Bi;
          effK = eff_K;
        }
      }
    }
    om = j;
    Result_tree.Fill();
    Tl = 1000;
    Bi = 1000;
    K = 1000;
    par1 = 1000;
    par2 = 1000;
    test = 1000;
    dt = 100;
  }

  newfile->cd();
  Result_tree.Write();
  newfile->Close();
  std::cout << "distrib done" << '\n';
}

void histo_mystere(){
  TFile *histo_file_Tl = new TFile("Histo_simu_new/MC_simu_Tl_208_10G_eres_53_gain_221_OC_MW8.root", "READ");
  histo_file_Tl->cd();
  TH3D* MC_Tl_208 = (TH3D*)histo_file_Tl->Get("MC_simu_Tl_208_10G_ubc");

  TFile *histo_file_Bi = new TFile("Histo_simu_new/MC_simu_Bi_214_10G_eres_53_gain_221_OC_MW8.root", "READ");
  histo_file_Bi->cd();
  TH3D* MC_Bi_214 = (TH3D*)histo_file_Bi->Get("MC_simu_Bi_214_10G_ubc");

  TFile *histo_file_K = new TFile("Histo_simu_new/MC_simu_K_40_10G_eres_53_gain_221_OC_MW8.root", "READ");
  histo_file_K->cd();
  TH3D* MC_K_40 = (TH3D*)histo_file_K->Get("MC_simu_K_40_10G_ubc");

  TH1D *mc0 = MC_Tl_208->ProjectionZ("Charge_Tl_208", 5, 5, 46, 46);    // first MC histogram
  TH1D *mc1 = MC_Bi_214->ProjectionZ("Charge_Bi_214", 5, 5, 46, 46);    // second MC histogram
  TH1D *mc2 = MC_K_40->ProjectionZ("Charge_K_40", 5, 5, 46, 46);    // second MC histogram

  mc0->Scale(1/mc0->Integral());
  mc1->Scale(1/mc1->Integral());
  mc2->Scale(1/mc2->Integral());

  mc0->Scale(0.2);
  mc1->Scale(0.5);
  mc2->Scale(0.3);

  mc0->Draw();
  mc1->Draw("same");
  mc2->Draw("same");

  TH1D MC("MC", "MC", 1024, 0, 200000);

  MC.Add(mc0);
  MC.Add(mc1);
  MC.Add(mc2);
  MC.Draw("same");
  TFile *file = new TFile("Histo_simu_new/Histo_mystere.root", "RECreate");
  file->cd();
  MC.Write();
  file->Close();

}

void eff_om_prep(string name) {

  double eff_tot =0;
  double eff_tot_error =0;
  int om = 0;
  std::vector<int> *om_id = new std::vector<int>;
  std::vector<double> *energy = new std::vector<double>;
  TFile *file = new TFile(Form("Histo_simu_new/Simu_%s_10G.root", name.c_str()), "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_id_new",1);
  tree->SetBranchAddress("om_id_new", &om_id);
  tree->SetBranchStatus("energy_u",1);
  tree->SetBranchAddress("energy_u", &energy);

  TFile *newfile = new TFile(Form("eff_om/eff_prep_%s.root", name.c_str()), "RECREATE");
  TTree new_tree("new_tree","");
  new_tree.Branch("eff_tot", &eff_tot);
  new_tree.Branch("eff_tot_error", &eff_tot_error);
  new_tree.Branch("om", &om);

  TH1D energy_id("e_id", "e id ", 712, 0, 712);
  TH1D energy_id_cut("e_id_cut", "e id cut ", 712, 0, 712);

  for (int i = 0; i < 1e5; i++) {
    // for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (i % 100000 == 0) {
      std::cout << "entry = " << i << '\n';
    }
    for (size_t k = 0; k < om_id->size(); k++) {
      for (size_t j = 0; j < energy->size(); j++) {
        energy_id.Fill(om_id->at(k)-1);
        om = om_id->at(k)-1;
        eff_tot = (energy->at(j))/5.0e9;
        eff_tot_error = (sqrt(energy->at(j))/5.0e9);
        if (energy->at(j) > 1) {
          energy_id_cut.Fill(om_id->at(k));
        }
        new_tree.Fill();
      }
    }
  }

  newfile->cd();
  new_tree.Write();
  energy_id.Write();
  energy_id_cut.Write();
  newfile->Close();
}

void eff_om(string name) {

  double eff_tot =0;
  double eff_tot_error =0;
  int om = 0;

  TFile *file = new TFile(Form("eff_om/eff_prep_%s.root", name.c_str()), "READ");
  TH1D* eff_om_prep = (TH1D*)file->Get("e_id");

  TH2F eff_tot_histo("eff", "efficacite mur ", 20, 0, 20, 13, 0, 13);

  TFile *newfile = new TFile(Form("eff_om/eff_om_%s.root", name.c_str()), "RECREATE");
  TTree new_tree("new_tree","");
  new_tree.Branch("eff_tot", &eff_tot);
  new_tree.Branch("eff_tot_error", &eff_tot_error);
  new_tree.Branch("om", &om);

  // for (int i = 0; i < eff_om_prep->GetMaximumBin(); i++) {
  for (int i = 0; i < 712; i++) {
    om = i;
    int om_col = (i % 13);
    int om_row = (i / 13);

    eff_tot = (eff_om_prep->GetBinContent(i))/5.0e9;
    eff_tot_error = (sqrt(eff_om_prep->GetBinContent(i))/5.0e9);

    eff_tot_histo.SetBinContent( om_row+1, om_col+1, eff_tot);

    new_tree.Fill();
  }

  file->Close();

  newfile->cd();

  new_tree.Write();
  eff_tot_histo.Write();

  newfile->Close();

}

void new_eff_om(string name) {
  std::vector<int> *om_id = new std::vector<int>;
  std::vector<double> *energy = new std::vector<double>;
  float vertex_position[3];
  double eff_tot =0;
  double eff_tot_error =0;
  int om = 0;
  TFile *file = new TFile(Form("Energie_OM_%s.root", name.c_str()), "READ");

  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_id",1);
  tree->SetBranchAddress("om_id", &om_id);
  tree->SetBranchStatus("energy",1);
  tree->SetBranchAddress("energy", &energy);

  TH1D energy_id("e_id", "e id ", 712, 0, 712);
  TFile *newfile = new TFile(Form("eff_om_%s.root", name.c_str()), "RECREATE");
  TTree new_tree("new_tree","");
  new_tree.Branch("eff_tot", &eff_tot);
  new_tree.Branch("eff_tot_error", &eff_tot_error);
  new_tree.Branch("om", &om);


  for (int i = 0; i < tree->GetEntries(); i++) {
    memset(vertex_position, 0, 3*sizeof(float));
    tree->GetEntry(i);
    if (i % 100000 == 0) {
      std::cout << "entry = " << i << '\n';
    }
    for (size_t k = 0; k < om_id->size(); k++) {
      for (size_t j = 0; j < energy->size(); j++) {
        energy_id.Fill(om_id->at(k)-1);
      }
    }
  }

  for (int i = 0; i < 712; i++) {
    om = i;
    eff_tot = (energy_id.GetBinContent(i))/5.0e8;
    eff_tot_error = (sqrt(energy_id.GetBinContent(i))/5.0e8);


    new_tree.Fill();
  }

  newfile->cd();
  new_tree.Write();
  newfile->Close();
  file->Close();
}

void fusion() {
  double om_flux_Tl, om_flux_Bi, om_flux_K, par1, par2, Tl, Bi, K, p1, p2, test, Gain, Chi2, gain, om_counting_Tl, om_counting_Bi, om_counting_K, counting_Tl, counting_Bi, counting_K, effTl, effBi, effK, eff_Tl, eff_Bi, eff_K;
  int om, om_id;
  float deadt, dt;
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om);
  Result_tree.Branch("Chi2", &test);
  Result_tree.Branch("gain", &Gain);
  Result_tree.Branch("om_flux_Tl", &om_flux_Tl);
  Result_tree.Branch("om_flux_Bi", &om_flux_Bi);
  Result_tree.Branch("om_flux_K", &om_flux_K);
  Result_tree.Branch("param1", &par1);
  Result_tree.Branch("param2", &par2);
  Result_tree.Branch("deadt", &dt);
  Result_tree.Branch("om_counting_Tl", &counting_Tl);
  Result_tree.Branch("om_counting_Bi", &counting_Bi);
  Result_tree.Branch("om_counting_K", &counting_K);
  Result_tree.Branch("eff_Tl", &effTl);
  Result_tree.Branch("eff_Bi", &effBi);
  Result_tree.Branch("eff_K", &effK);

  TFile *file1 = new TFile("histo_fit/distrib_fit_it.root", "READ");
  TTree* tree = (TTree*)file1->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_id);
  tree->SetBranchStatus("om_flux_Tl",1);
  tree->SetBranchAddress("om_flux_Tl", &Tl);
  tree->SetBranchStatus("om_flux_Bi",1);
  tree->SetBranchAddress("om_flux_Bi", &Bi);
  tree->SetBranchStatus("om_flux_K",1);
  tree->SetBranchAddress("om_flux_K", &K);
  tree->SetBranchStatus("param1",1);
  tree->SetBranchAddress("param1", &p1);
  tree->SetBranchStatus("param2",1);
  tree->SetBranchAddress("param2", &p2);
  tree->SetBranchStatus("Chi2",1);
  tree->SetBranchAddress("Chi2", &Chi2);
  tree->SetBranchStatus("gain",1);
  tree->SetBranchAddress("gain", &gain);
  tree->SetBranchStatus("deadt",1);
  tree->SetBranchAddress("deadt", &deadt);
  tree->SetBranchStatus("om_counting_Tl",1);
  tree->SetBranchAddress("om_counting_Tl", &om_counting_Tl);
  tree->SetBranchStatus("om_counting_Bi",1);
  tree->SetBranchAddress("om_counting_Bi", &om_counting_Bi);
  tree->SetBranchStatus("om_counting_K",1);
  tree->SetBranchAddress("om_counting_K", &om_counting_K);
  tree->SetBranchStatus("eff_Tl",1);
  tree->SetBranchAddress("eff_Tl", &eff_Tl);
  tree->SetBranchStatus("eff_Bi",1);
  tree->SetBranchAddress("eff_Bi", &eff_Bi);
  tree->SetBranchStatus("eff_K",1);
  tree->SetBranchAddress("eff_K", &eff_K);

  for (size_t i = 0; i < 260; i++) {
    tree->GetEntry(i);
    dt = deadt;
    om = om_id;
    om_flux_Tl = Tl;
    om_flux_Bi = Bi;
    om_flux_K = K;
    par1 = p1;
    par2 = p2;
    Gain = gain;
    test = Chi2;
    counting_Tl = om_counting_Tl;
    counting_Bi = om_counting_Bi;
    counting_K = om_counting_K;
    effTl = eff_Tl;
    effBi = eff_Bi;
    effK = eff_K;
    Result_tree.Fill();
  }

  deadt = om_id = Tl = Bi = K = p1 = p2 = gain = Chi2 = 0;

  TFile *file2 = new TFile("histo_fit/distrib_fit_fr.root", "READ");
  TTree* tree2 = (TTree*)file2->Get("Result_tree");
  tree2->SetBranchStatus("*",0);
  tree2->SetBranchStatus("om_number",1);
  tree2->SetBranchAddress("om_number", &om_id);
  tree2->SetBranchStatus("om_flux_Tl",1);
  tree2->SetBranchAddress("om_flux_Tl", &Tl);
  tree2->SetBranchStatus("om_flux_Bi",1);
  tree2->SetBranchAddress("om_flux_Bi", &Bi);
  tree2->SetBranchStatus("om_flux_K",1);
  tree2->SetBranchAddress("om_flux_K", &K);
  tree2->SetBranchStatus("param1",1);
  tree2->SetBranchAddress("param1", &p1);
  tree2->SetBranchStatus("param2",1);
  tree2->SetBranchAddress("param2", &p2);
  tree2->SetBranchStatus("Chi2",1);
  tree2->SetBranchAddress("Chi2", &Chi2);
  tree2->SetBranchStatus("gain",1);
  tree2->SetBranchAddress("gain", &gain);
  tree2->SetBranchStatus("deadt",1);
  tree2->SetBranchAddress("deadt", &deadt);
  tree2->SetBranchStatus("om_counting_Tl",1);
  tree2->SetBranchAddress("om_counting_Tl", &om_counting_Tl);
  tree2->SetBranchStatus("om_counting_Bi",1);
  tree2->SetBranchAddress("om_counting_Bi", &om_counting_Bi);
  tree2->SetBranchStatus("om_counting_K",1);
  tree2->SetBranchAddress("om_counting_K", &om_counting_K);
  tree2->SetBranchStatus("eff_Tl",1);
  tree2->SetBranchAddress("eff_Tl", &eff_Tl);
  tree2->SetBranchStatus("eff_Bi",1);
  tree2->SetBranchAddress("eff_Bi", &eff_Bi);
  tree2->SetBranchStatus("eff_K",1);
  tree2->SetBranchAddress("eff_K", &eff_K);

  for (size_t i = 260; i < 520; i++) {
    tree2->GetEntry(i);
    dt = deadt;
    om = om_id;
    om_flux_Tl = Tl;
    om_flux_Bi = Bi;
    om_flux_K = K;
    par1 = p1;
    par2 = p2;
    Gain = gain;
    test = Chi2;
    counting_Tl = om_counting_Tl;
    counting_Bi = om_counting_Bi;
    counting_K = om_counting_K;
    effTl = eff_Tl;
    effBi = eff_Bi;
    effK = eff_K;
    Result_tree.Fill();
  }

  deadt = om_id = Tl = Bi = K = p1 = p2 = gain = Chi2 = 0;

  TFile *file3 = new TFile("histo_fit/distrib_fit_xw.root", "READ");
  TTree* tree3 = (TTree*)file3->Get("Result_tree");
  tree3->SetBranchStatus("*",0);
  tree3->SetBranchStatus("om_number",1);
  tree3->SetBranchAddress("om_number", &om_id);
  tree3->SetBranchStatus("om_flux_Tl",1);
  tree3->SetBranchAddress("om_flux_Tl", &Tl);
  tree3->SetBranchStatus("om_flux_Bi",1);
  tree3->SetBranchAddress("om_flux_Bi", &Bi);
  tree3->SetBranchStatus("om_flux_K",1);
  tree3->SetBranchAddress("om_flux_K", &K);
  tree3->SetBranchStatus("param1",1);
  tree3->SetBranchAddress("param1", &p1);
  tree3->SetBranchStatus("param2",1);
  tree3->SetBranchAddress("param2", &p2);
  tree3->SetBranchStatus("Chi2",1);
  tree3->SetBranchAddress("Chi2", &Chi2);
  tree3->SetBranchStatus("gain",1);
  tree3->SetBranchAddress("gain", &gain);
  tree3->SetBranchStatus("deadt",1);
  tree3->SetBranchAddress("deadt", &deadt);
  tree3->SetBranchStatus("om_counting_Tl",1);
  tree3->SetBranchAddress("om_counting_Tl", &om_counting_Tl);
  tree3->SetBranchStatus("om_counting_Bi",1);
  tree3->SetBranchAddress("om_counting_Bi", &om_counting_Bi);
  tree3->SetBranchStatus("om_counting_K",1);
  tree3->SetBranchAddress("om_counting_K", &om_counting_K);
  tree3->SetBranchStatus("eff_Tl",1);
  tree3->SetBranchAddress("eff_Tl", &eff_Tl);
  tree3->SetBranchStatus("eff_Bi",1);
  tree3->SetBranchAddress("eff_Bi", &eff_Bi);
  tree3->SetBranchStatus("eff_K",1);
  tree3->SetBranchAddress("eff_K", &eff_K);

  for (size_t i = 520; i < 648; i++) {
    tree3->GetEntry(i);
    dt = deadt;
    om = om_id;
    om_flux_Tl = Tl;
    om_flux_Bi = Bi;
    om_flux_K = K;
    par1 = p1;
    par2 = p2;
    Gain = gain;
    test = Chi2;
    counting_Tl = om_counting_Tl;
    counting_Bi = om_counting_Bi;
    counting_K = om_counting_K;
    effTl = eff_Tl;
    effBi = eff_Bi;
    effK = eff_K;
    Result_tree.Fill();

  }

  deadt = om_id = Tl = Bi = K = p1 = p2 = gain = Chi2 = 0;

  TFile *file4 = new TFile("histo_fit/distrib_fit_gv.root", "READ");
  TTree* tree4 = (TTree*)file4->Get("Result_tree");
  tree4->SetBranchStatus("*",0);
  tree4->SetBranchStatus("om_number",1);
  tree4->SetBranchAddress("om_number", &om_id);
  tree4->SetBranchStatus("om_flux_Tl",1);
  tree4->SetBranchAddress("om_flux_Tl", &Tl);
  tree4->SetBranchStatus("om_flux_Bi",1);
  tree4->SetBranchAddress("om_flux_Bi", &Bi);
  tree4->SetBranchStatus("om_flux_K",1);
  tree4->SetBranchAddress("om_flux_K", &K);
  tree4->SetBranchStatus("param1",1);
  tree4->SetBranchAddress("param1", &p1);
  tree4->SetBranchStatus("param2",1);
  tree4->SetBranchAddress("param2", &p2);
  tree4->SetBranchStatus("Chi2",1);
  tree4->SetBranchAddress("Chi2", &Chi2);
  tree4->SetBranchStatus("gain",1);
  tree4->SetBranchAddress("gain", &gain);
  tree4->SetBranchStatus("deadt",1);
  tree4->SetBranchAddress("deadt", &deadt);
  tree4->SetBranchStatus("om_counting_Tl",1);
  tree4->SetBranchAddress("om_counting_Tl", &om_counting_Tl);
  tree4->SetBranchStatus("om_counting_Bi",1);
  tree4->SetBranchAddress("om_counting_Bi", &om_counting_Bi);
  tree4->SetBranchStatus("om_counting_K",1);
  tree4->SetBranchAddress("om_counting_K", &om_counting_K);
  tree4->SetBranchStatus("eff_Tl",1);
  tree4->SetBranchAddress("eff_Tl", &eff_Tl);
  tree4->SetBranchStatus("eff_Bi",1);
  tree4->SetBranchAddress("eff_Bi", &eff_Bi);
  tree4->SetBranchStatus("eff_K",1);
  tree4->SetBranchAddress("eff_K", &eff_K);

  for (size_t i = 648; i < 712; i++) {
    tree4->GetEntry(i);
    dt = deadt;
    om = om_id;
    om_flux_Tl = Tl;
    om_flux_Bi = Bi;
    om_flux_K = K;
    par1 = p1;
    par2 = p2;
    Gain = gain;
    test = Chi2;
    counting_Tl = om_counting_Tl;
    counting_Bi = om_counting_Bi;
    counting_K = om_counting_K;
    effTl = eff_Tl;
    effBi = eff_Bi;
    effK = eff_K;
    Result_tree.Fill();

  }

  TFile *newfile = new TFile("histo_fit/fusion2.root", "RECREATE");
  newfile->cd();
  Result_tree.Write();
  newfile->Close();
  std::cout << "ok " << '\n';

}

void test_param(string input_file) {

  double p1, p2, gain;
  int om_number;
  std::vector<int> *om_id = new std::vector<int>;
  std::vector<int> *param_value = new std::vector<int>;

  TFile *file = new TFile(Form("histo_fit/%s.root", input_file.c_str()), "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om_number);
  tree->SetBranchStatus("param1",1);
  tree->SetBranchAddress("param1", &p1);
  tree->SetBranchStatus("param2",1);
  tree->SetBranchAddress("param2", &p2);
  tree->SetBranchStatus("gain",1);
  tree->SetBranchAddress("gain", &gain);

  for (int i = 0; i < 712; i++) {
    tree->GetEntry(i);
    if (p1 < 0.105) {
      om_id->push_back(om_number);
      param_value->push_back(1);
    }
    else if (p2 > 0.598) {
      om_id->push_back(om_number);
      param_value->push_back(2);
    }
  }
  // Fit_Gain_Simu();


}




// int main(int argc, char const *argv[]){
//   for(int i = 0; i<argc; i++){
//     if (std::string(argv[i]) == "-XW" || std::string(argv[i]) =="-MWFR" || std::string(argv[i]) =="-MWIT" || std::string(argv[i]) =="-GV"){
//       Fit_Gain_Simu(argv[i]);
//       if(std::string(argv[i]) == "--distrib" || std::string(argv[i]) =="-d"){
//         distrib("roofit");
//       }
//     }
//   }
//   return 0;
// }

void Fit_Gain_Res(string wall = "") {
  LoadMC();
  Load_spectre();
  Load_dead_spectre();
  TH1::SetDefaultSumw2();

  double logL, Chi2, param1, param2, param3, om_counting_Tl, om_counting_Bi, om_counting_K, eff_Tl, eff_Bi, eff_K, kolmo, eff_error_Tl, eff_error_Bi, eff_error_K;
  double int_cut_mc0, int_cut_mc1, int_cut_mc2, int_full_mc0, int_full_mc1, int_full_mc2;
  float gain, eres = 0;
  int om_number = 0;
  double mean_erf = 0;
  double om_flux_Tl, om_flux_Bi, om_flux_K;
  int lim =0;
  float deadt;
  double test;

  TFile *gainfile = new TFile("histo_fit/distrib_test_.root", "READ");
  TTree* tree = (TTree*)gainfile->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("gain",1);
  tree->SetBranchAddress("gain", &test);
  double gain_tab[712];
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    gain_tab[i] = test;
  }

  TFile *newfile = new TFile(Form("histo_fit/fit_eres_1%s.root", wall.c_str()), "RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("logL", &logL);
  Result_tree.Branch("Chi2", &Chi2);
  Result_tree.Branch("kolmo", &kolmo);
  Result_tree.Branch("param1", &param1);
  Result_tree.Branch("param2", &param2);
  Result_tree.Branch("param3", &param3);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("eres", &eres);
  Result_tree.Branch("om_number", &om_number);
  Result_tree.Branch("om_flux_Tl", &om_flux_Tl);
  Result_tree.Branch("om_flux_Bi", &om_flux_Bi);
  Result_tree.Branch("om_flux_K", &om_flux_K);
  Result_tree.Branch("lim", &lim);
  Result_tree.Branch("deadt", &deadt);
  Result_tree.Branch("om_counting_Tl", &om_counting_Tl);
  Result_tree.Branch("om_counting_Bi", &om_counting_Bi);
  Result_tree.Branch("om_counting_K", &om_counting_K);
  Result_tree.Branch("eff_Tl", &eff_Tl);
  Result_tree.Branch("eff_Bi", &eff_Bi);
  Result_tree.Branch("eff_K", &eff_K);
  Result_tree.Branch("eres", &eres);
  TH2D gain_eres_Khi("gain_eres_Khi","gain_eres_Khi", 221, 10000,65000, 53, 7, 20);
  TH2D gain_eres_kolmo("gain_eres_kolmo","gain_eres_kolmo", 221, 10000,65000, 53, 7, 20);

  std::ifstream  charge("gain_data/Resultat_mean_energie_charge_total.txt");
  float charge_valeur_fit [712];
  memset(charge_valeur_fit, -1, 712*sizeof(float));
  int charge_om_num;
  while (charge >> charge_om_num)
  {
    charge >> charge_valeur_fit[charge_om_num];
  }
  double* rootab = new double[6];
  double* eff = new double[2];

  float min = 0.85;
  float max = 1.15;

  int* borne = new int[2];
  borne = om_chooser(wall);

  for (int om = 1; om < 2; om++)
  // for (int om = borne[0]; om < borne[1]; om = om +1)
  {
    TH3D* MC_Tl_208 = MC_chooser(om, 0);
    TH3D* MC_Bi_214 = MC_chooser(om, 1);
    TH3D* MC_K_40 = MC_chooser(om, 2);
    TH1D* spectre_om = NULL;
    spectre_om = spectre_charge_full(om);

    if (om < 520) {
      if ((spectre_om->GetEntries() < 100) || (spectre_om->GetMean(1) < 1500)){
        delete spectre_om;
        om++;
        spectre_om = spectre_charge_full(om);
      }
      if (charge_valeur_fit[om] != -1){
        double *tab = om_gain_fit(om);
        if (tab[0] == 0 ) {
          mean_erf = charge_valeur_fit[om];
        }
        else{
          mean_erf = 1/tab[0];
        }
        delete tab;
      }
    }
    else if (om < 648 && om > 519) {
      mean_erf = charge_valeur_fit[om];
    }
    else{
      mean_erf = 1.0/19000;
      min = 0.8;
      max = 1.5;
    }
    om_number = om;

    int eres_count = (eres-eres_bin_min)/eres_bin_width+1;
    gain = gain_tab[om];
    int gain_count =((gain - 10000)/250) + 1;
    // for (int gain_count = 0; gain_count <150; gain_count++) {
      for (int eres_count = 0; eres_count <53; eres_count++) {
        std::cout << (gain_bin_min + gain_bin_width*(gain_count-1))<< "   et    lim_inf = " << 1/(mean_erf*min) << "   sup  ="  << 1/(mean_erf*max) << '\n';
        if ((1.0/(gain_bin_min + gain_bin_width*(gain_count-1)) > mean_erf*min) && (1.0/(gain_bin_min + gain_bin_width*(gain_count-1))<mean_erf*max)){
          spectre_om = spectre_charge_full(om);
          eres = eres_bin_min + eres_bin_width*(eres_count-1);
          gain = (gain_bin_min + gain_bin_width*(gain_count-1));
          std::cout << "gain_count = " << gain_count << " and gain = " << gain << '\n';

          TH1D *mc0 = MC_Tl_208->ProjectionZ("Charge_Tl_208", eres_count, eres_count, gain_count, gain_count);    // first MC histogram
          TH1D *mc1 = MC_Bi_214->ProjectionZ("Charge_Bi_214", eres_count, eres_count, gain_count, gain_count);    // second MC histogram
          TH1D *mc2 = MC_K_40->ProjectionZ("Charge_K_40", eres_count, eres_count, gain_count, gain_count);    // second MC histogram
          int bin = gain*(1024/200000.0);

          if (om == 570) {
            delete spectre_om;
            TFile *tfile = new TFile("histo_kolmo/histo_donee/histo_charge_amplitude_energie_611_570.root", "READ");
            gROOT->cd();
            TH2F* charge_spectre_570 = (TH2F*)tfile->Get("histo_pm_charge");
            TH1D* spectre_om = charge_spectre_570->ProjectionY("charge570", 571, 571);
          }

          int_full_mc0 = mc0->Integral();
          int_full_mc1 = mc2->Integral();
          int_full_mc2 = mc2->Integral();

          for (int i = 0; i < bin; i++) {
          mc0->SetBinContent(i, 0);
          mc0->SetBinError(i, 0);
          mc1->SetBinContent(i, 0);
          mc2->SetBinError(i, 0);
          mc2->SetBinContent(i, 0);
          mc2->SetBinError(i, 0);
          spectre_om->SetBinContent(i, 0);
          spectre_om->SetBinError(i, 0);
          }

          int_cut_mc0 = mc0->Integral();
          int_cut_mc1 = mc1->Integral();
          int_cut_mc2 = mc2->Integral();

          roofitter(om, gain, eres, spectre_om, mc0, mc1, mc2, rootab, bin);
          logL = rootab[0];
          Chi2 = rootab[4];
          param1 = rootab[1];
          param2 = rootab[2];
          param3 = rootab[3];
          kolmo = rootab[5];
          deadt = dead_time(om, bin);
          gain_eres_Khi.SetBinContent(gain_count, eres_count, Chi2);
          gain_eres_kolmo.SetBinContent(gain_count, eres_count, kolmo);
          if (deadt < 1) {
            deadt = 1;
          }
          //
          // mc0 = MC_Tl_208->ProjectionZ("Charge_Tl_208", eres_count, eres_count, gain_count, gain_count);
          // mc1 = MC_Bi_214->ProjectionZ("Charge_Bi_214", eres_count, eres_count, gain_count, gain_count);
          // mc2 = MC_K_40->ProjectionZ("Charge_K_40", eres_count, eres_count, gain_count, gain_count);
          // spectre_om = spectre_charge_full(om);
          //
          // for (int i = 0; i < bin; i++) {
          //   mc0->SetBinContent(i, 0);
          //   mc1->SetBinContent(i, 0);
          //   mc2->SetBinContent(i, 0);
          //   spectre_om->SetBinContent(i, 0);
          //   spectre_om->SetBinError(i, 0);
          // }
          //
          // mc0->Scale((param1/int_cut_mc0)*spectre_om->Integral());
          // mc1->Scale(spectre_om->Integral()*param2/int_cut_mc1);
          // mc2->Scale(spectre_om->Integral()*param3/int_cut_mc2);
          //
          // TH1D *sumfit = new TH1D(*mc0);
          // sumfit->Add(mc1);
          // sumfit->Add(mc2);
          //
          // sumfit->Draw();
          // spectre_om->Draw("same");
          //
          // double kolmo2 = spectre_om->KolmogorovTest(sumfit);
          // std::cout << kolmo2 << '\n';

          delete mc0;
          delete mc1;
          delete mc2;
          delete spectre_om;
          mc0 = MC_Tl_208->ProjectionZ("Charge_Tl_208", eres_count, eres_count, gain_count, gain_count);
          mc1 = MC_Bi_214->ProjectionZ("Charge_Bi_214", eres_count, eres_count, gain_count, gain_count);
          mc2 = MC_K_40->ProjectionZ("Charge_K_40", eres_count, eres_count, gain_count, gain_count);
          spectre_om = spectre_charge_full(om);
          get_om_eff("Tl_208", om, eff);
          eff_Tl = eff[0];
          eff_error_Tl = eff[1];
          get_om_eff("Tl_208", om, eff);
          eff_Bi = eff[0];
          eff_error_Bi = eff[1];
          get_om_eff("Tl_208", om, eff);
          eff_K = eff[0];
          eff_error_K = eff[1];

          mc0->Scale((param1/int_cut_mc0)*spectre_om->Integral());
          om_counting_Tl = mc0->Integral();
          om_flux_Tl = (mc0->Integral()*deadt)/(902*eff_Tl*250000);
          mc1->Scale(spectre_om->Integral()*param2/int_cut_mc1);
          om_counting_Bi = mc1->Integral();
          om_flux_Bi = (mc1->Integral()*deadt)/(902*eff_Bi*250000);
          mc2->Scale(spectre_om->Integral()*param3/int_cut_mc2);
          om_counting_K = mc2->Integral();
          om_flux_K = (mc2->Integral()*deadt)/(902*eff_K*250000);

          Result_tree.Fill();
          delete mc0;
          delete mc1;
          delete mc2;
          delete spectre_om;
        }
      }
    // }
  }
  newfile->cd();
  Result_tree.Write();
  gain_eres_Khi.Write();
  gain_eres_kolmo.Write();
  newfile->Close();
  std::cout << "good end" << '\n';
}
