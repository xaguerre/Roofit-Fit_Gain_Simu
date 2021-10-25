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
#include "RooProdPdf.h"
#include "RooPlot.h"
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

TH1D* spectre_charge_it(int om_number){
  TFile *file = new TFile("histo_kolmo/histo_donee/histo_charge_amplitude_422.root", "READ");
  gROOT->cd();
  TH2F* charge = (TH2F*)file->Get("charge");

  TH1D* spectre_charge = charge->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
  file->Close();
  return spectre_charge;
}

TH1D* spectre_charge_XW(int om_number){
  TFile *file = new TFile("histo_kolmo/histo_donee/histo_charge_amplitude_energie_449.root", "READ");
  gROOT->cd();
  TH2F* charge = (TH2F*)file->Get("histo_pm_charge");

  TH1D* spectre_charge = charge->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
  file->Close();
  return spectre_charge;
}

TH1D* spectre_charge_fr(int om_number){
  TFile *file = new TFile("histo_kolmo/histo_donee/histo_charge_amplitude_energie_435.root", "READ");
  gROOT->cd();
  TH2F* charge = (TH2F*)file->Get("histo_pm_charge");

  TH1D* spectre_charge = charge->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
  file->Close();
  return spectre_charge;
}

TH1D* spectre_chooser(int om){
  TH1D* spectre_chooser = NULL;
  if (om <260){
    spectre_chooser = spectre_charge_it(om);
  }
  else if (om <520) {
    spectre_chooser = spectre_charge_fr(om);
  }
  else if (om >519) {
    spectre_chooser = spectre_charge_XW(om);
  }
  return spectre_chooser;
}

TH3D *MC_MW8[3] = {NULL, NULL, NULL};
TH3D *MC_MW5[3] = {NULL, NULL, NULL};
TH3D *MC_XW[3]  = {NULL, NULL, NULL};
TH3D *MC_GV[3]  = {NULL, NULL, NULL};

void LoadMC() {
  TFile *histo_file_Tl_MW8 = new TFile("Histo_simu_new/MC_simu_Tl_208_recOC_new_eres_53_gain_221_OC_MW8.root", "READ");
  histo_file_Tl_MW8->cd();
  MC_MW8[0] = (TH3D*)histo_file_Tl_MW8->Get("MC_simu_Tl_208_recOC_new_ubc");
  TFile *histo_file_Bi_MW8 = new TFile("Histo_simu_new/MC_simu_Bi_214_recOC_new_eres_53_gain_221_OC_MW8.root", "READ");
  histo_file_Bi_MW8->cd();
  MC_MW8[1] = (TH3D*)histo_file_Bi_MW8->Get("MC_simu_Bi_214_recOC_new_ubc");
  TFile *histo_file_K_MW8 = new TFile("Histo_simu_new/MC_simu_K_40_recOC_new_eres_53_gain_221_OC_MW8.root", "READ");
  histo_file_K_MW8->cd();
  MC_MW8[2] = (TH3D*)histo_file_K_MW8->Get("MC_simu_K_40_recOC_new_ubc");

  TFile *histo_file_Tl_MW5 = new TFile("Histo_simu_new/MC_simu_Tl_208_recOC_new_eres_53_gain_221_OC_MW5.root", "READ");
  histo_file_Tl_MW5->cd();
  MC_MW5[0] = (TH3D*)histo_file_Tl_MW5->Get("MC_simu_Tl_208_recOC_new_ubc");
  TFile *histo_file_Bi_MW5 = new TFile("Histo_simu_new/MC_simu_Bi_214_recOC_new_eres_53_gain_221_OC_MW5.root", "READ");
  histo_file_Bi_MW5->cd();
  MC_MW5[1] = (TH3D*)histo_file_Bi_MW5->Get("MC_simu_Bi_214_recOC_new_ubc");
  TFile *histo_file_K_MW5 = new TFile("Histo_simu_new/MC_simu_K_40_recOC_new_eres_53_gain_221_OC_MW5.root", "READ");
  histo_file_K_MW5->cd();
  MC_MW5[2] = (TH3D*)histo_file_K_MW5->Get("MC_simu_K_40_recOC_new_ubc");

  TFile *histo_file_TL_XW = new TFile("Histo_simu_new/MC_simu_Tl_208_recOC_new_eres_53_gain_221_OC_XW.root", "READ");
  histo_file_TL_XW->cd();
  MC_XW[0] = (TH3D*)histo_file_TL_XW->Get("MC_simu_Tl_208_recOC_new_ubc");
  TFile *histo_file_Bi_XW = new TFile("Histo_simu_new/MC_simu_Bi_214_recOC_new_eres_53_gain_221_OC_XW.root", "READ");
  histo_file_Bi_XW->cd();
  MC_XW[1] = (TH3D*)histo_file_Bi_XW->Get("MC_simu_Bi_214_recOC_new_ubc");
  TFile *histo_file_K_XW = new TFile("Histo_simu_new/MC_simu_K_40_recOC_new_eres_53_gain_221_OC_XW.root", "READ");
  histo_file_K_XW->cd();
  MC_XW[2] = (TH3D*)histo_file_K_XW->Get("MC_simu_K_40_recOC_new_ubc");

  TFile *histo_file_Tl_GV = new TFile("Histo_simu_new/MC_simu_Tl_208_recOC_new_eres_53_gain_221_OC_GV.root", "READ");
  histo_file_Tl_GV->cd();
  MC_GV[0] = (TH3D*)histo_file_Tl_GV->Get("MC_simu_Tl_208_recOC_new_ubc");
  TFile *histo_file_Bi_GV = new TFile("Histo_simu_new/MC_simu_Bi_214_recOC_new_eres_53_gain_221_OC_GV.root", "READ");
  histo_file_Bi_GV->cd();
  MC_GV[1] = (TH3D*)histo_file_Bi_GV->Get("MC_simu_Bi_214_recOC_new_ubc");
  TFile *histo_file_K_GV = new TFile("Histo_simu_new/MC_simu_K_40_recOC_new_eres_53_gain_221_OC_GV.root", "READ");
  histo_file_K_GV->cd();
  MC_GV[2] = (TH3D*)histo_file_K_GV->Get("MC_simu_K_40_recOC_new_ubc");
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
  if (om <260){
    spectre_om = spectre_charge_it(om);
  }
  else if (om <520) {
    spectre_om = spectre_charge_fr(om);
  }
  else if (om >519) {
    spectre_om = spectre_charge_XW(om);
  }

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

  if (chin < 1.5) {
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
  if(om == "-MW"  ){
    lim[0] = 0;
    lim[1] = 520;
  }
  if(om == "-XW"  ){
    lim[0] = 520;
    lim[1] = 648;
  }
  if(om == "-GV"  ){
    lim[0] = 648;
    lim[1] = 712;
  }
  if(om == "-FULL"  ){
    lim[0] = 0;
    lim[1] = 712;
  }
  return lim;
}

double* roofitter(int om, double gain, double eres, TH1D* spectre_om, TH1D* mc0, TH1D* mc1, TH1D* mc2, double *rootab, int bin){
  const float start_x = spectre_om->GetXaxis()->GetBinUpEdge(bin);

  mc0->Scale(1/mc0->Integral());
  mc1->Scale(1/mc1->Integral());
  mc2->Scale(1/mc2->Integral());

  RooRealVar x("x", "x", 0, 100000);
  RooDataHist Tl("Tl", "Tl", x, Import(*mc0));
  RooDataHist Bi("Bi", "Bi", x, Import(*mc1));
  RooDataHist K("K", "K", x, Import(*mc2));

  RooHistPdf Tl_pdf ("Tl_pdf", "", x, Tl);
  RooHistPdf Bi_pdf ("Bi_pdf", "", x, Bi);
  RooHistPdf K_pdf ("K_pdf", "", x, K);

  RooDataHist spectre_data("spectre_data", "spectre_data", x, Import(*spectre_om));

  RooParamHistFunc p_ph_Tl("p_ph_Tl","p_ph_Tl",Tl);
  RooParamHistFunc p_ph_Bi("p_ph_Bi","p_ph_Bi",Bi);
  RooParamHistFunc p_ph_K("p_ph_K","p_ph_K",K);

  RooRealVar RooTl("RooTl", "Tl", 0.2, 0.1, 0.4);
  RooRealVar RooBi("RooBi", "Bi", 0.5, 0.3, 0.6);
  RooRealVar RooK("RooK", "K", 0.3, 0.2, 0.4);

  RooRealSumPdf sum_simu("sum_simu", "sum_simu",
  RooArgList(p_ph_Tl, p_ph_Bi, p_ph_K),
  RooArgList(RooTl,RooBi,RooK),
  false);

  RooFitResult* result1 = sum_simu.fitTo(spectre_data, PrintLevel(-1), SumW2Error(true), Range(start_x,80000), Save(), IntegrateBins(-1), Minimizer("Minuit"));
  sum_simu.setStringAttribute("fitrange", nullptr);
  TCanvas* can = new TCanvas;
  can->cd();
  auto frame = x.frame(Title("Fit gain simu"));
  // Plot data to enable automatic determination of model0 normalisation:
  sum_simu.plotOn(frame, FillColor(0), VisualizeError(*result1));

  spectre_data.plotOn(frame);
  // Plot data again to show it on top of model0 error bands:
  spectre_data.plotOn(frame, Name("spectre_hist"));
  sum_simu.plotOn(frame, LineColor(kRed), Name("sum_curve"));
  rootab[4] = frame->chiSquare();

  // Plot model components
  sum_simu.plotOn(frame, Components(p_ph_Tl), LineColor(kGreen), Name("Tl_curve"));
  sum_simu.plotOn(frame, Components(p_ph_Bi), LineColor(kYellow), Name("Bi_curve"));
  sum_simu.plotOn(frame, Components(p_ph_K), LineColor(kBlue), Name("K_curve"));
  // sum_simu.paramOn(frame);

  RooCurve * Tl_curve = frame->getCurve("Tl_curve");
  double Tl_int = Tl_curve->average(start_x,65000);
  RooCurve * Bi_curve = frame->getCurve("Bi_curve");
  double Bi_int = Bi_curve->average(start_x,65000);
  RooCurve * K_curve = frame->getCurve("K_curve");
  double K_int = K_curve->average(start_x,65000);
  RooCurve * sum_curve = frame->getCurve("sum_curve");
  double int_data = sum_curve->average(start_x,65000);


  // TLatex latex;
  // latex.DrawLatex(.2,.9,"K_{S}");
  frame->GetYaxis()->UnZoom();
  frame->GetYaxis()->SetTitle("n events");
  frame->GetXaxis()->SetTitle("charge (u.a)");
  frame ->Draw();
  can->SetLogy();
  TLatex l;
  l.SetTextFont(40);
  l.DrawLatex(60000, 5000, Form("Khi2/NDF = %.3f",rootab[4]));
  // l->DrawLatex(60000, 3000, Form("loglikelihood = %f",rootab[0]));
  // TLine* xline = new TLine(58000, 2500, 100000, 2500);
  // xline->SetLineWidth(2);
  // xline->Draw();
  // TLine* yline = new TLine(58000, 2500, 58000, 9500);
  // yline->SetLineWidth(2);
  // yline->Draw();

  can->SaveAs(Form("OM_fit/om_%d/best_fit_om_%d_eres_%.2f_gain_%.0f_bin_%d.png", om, om, eres, gain, bin));
  // result1->Print("v");
  // result1->Print();
  rootab[0] = result1->minNll();
  rootab[1] = Tl_int/int_data;
  rootab[2] = Bi_int/int_data;
  rootab[3] = K_int/int_data;

  delete sum_curve;
  delete K_curve;
  delete Bi_curve;
  delete Tl_curve;
  delete can;
  delete frame;
  delete result1;
  return rootab;
}

double get_om_eff(string compo, int om) {
  double om_eff = 0;
  TFile *newfile = new TFile(Form("eff_om_%s.root", compo.c_str()), "READ");

  TTree* tree = (TTree*)newfile->Get("new_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("eff_tot",1);
  tree->SetBranchAddress("eff_tot", &om_eff);

  tree->GetEntry(om);
  newfile->Close();
  return om_eff;
}

void Fit_Gain_Simu(string wall) {
  LoadMC();
  TH1::SetDefaultSumw2();

  double logL, Chi2, param1, param2, param3 = 0;
  float gain, eres = 0;
  int om_number = 0;
  double mean_erf = 0;
  double om_flux_Tl, om_flux_Bi, om_flux_K;
  int lim =0;

  TFile *newfile = new TFile("histo_fit/histo_gv2.root", "RECREATE");
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
  Result_tree.Branch("om_flux_K", &om_flux_K);
  Result_tree.Branch("lim", &lim);

  std::ifstream  charge("gain_data/Resultat_mean_energie_charge_total.txt");
  float charge_valeur_fit [712];
  memset(charge_valeur_fit, -1, 712*sizeof(float));
  int charge_om_num;
  while (charge >> charge_om_num)
  {
    charge >> charge_valeur_fit[charge_om_num];
  }
  double* rootab = new double[5];

  // TFile *histo = new TFile("Histo_simu_new/Histo_mystere.root", "READ");
  // histo->cd();
  // TH1D* MC = (TH1D*)histo->Get("MC");

  float min = 0.85;
  float max = 1.15;

  // int* borne = om_chooser(wall);

  for (int om = 648; om <712; om++)
  // for (int om = borne[0]; om < borne[1]; om = om +1)
  {

    TH3D* MC_Tl_208 = MC_chooser(om, 0);
    TH3D* MC_Bi_214 = MC_chooser(om, 1);
    TH3D* MC_K_40 = MC_chooser(om, 2);

    TH1D* spectre_om = NULL;
    spectre_om = spectre_chooser(om);

    if (om < 520) {
      if ((spectre_om->GetEntries() < 100) || (spectre_om->GetMean(1) < 1500)){
        delete spectre_om;
        om++;
        spectre_om = spectre_chooser(om);
      }
      if (charge_valeur_fit[om] != -1){
        double *tab = om_gain_fit(om);
        if (tab[0] == 0) {
          mean_erf = charge_valeur_fit[om];
        }
        else{
          mean_erf = 1/tab[0];
        }
        std::cout << "mean erf = " << mean_erf << '\n';
        delete tab;
      }
    }
    else if (om < 648 && om > 519) {
      mean_erf = charge_valeur_fit[om];
    }
    else{
      mean_erf = 1.0/22000;
      min = 0.8;
      max = 1.8;
    }
    om_number = om;
    int eres = eres_chooser(om);
    int eres_count = (eres-eres_bin_min)/eres_bin_width+1;
    for (int gain_count = 0; gain_count <150; gain_count++) {
      std::cout << (gain_bin_min + gain_bin_width*(gain_count-1))<< "   et    lim_inf = " << 1/(mean_erf*min) << "   sup  ="  << 1/(mean_erf*max) << '\n';
      if ((1.0/(gain_bin_min + gain_bin_width*(gain_count-1)) > mean_erf*min) && (1.0/(gain_bin_min + gain_bin_width*(gain_count-1))<mean_erf*max)){

        gain = (gain_bin_min + gain_bin_width*(gain_count-1));
        std::cout << "gain_count = " << gain_count << " and gain = " << gain << '\n';

        TH1D *mc0 = MC_Tl_208->ProjectionZ("Charge_Tl_208", eres_count, eres_count, gain_count, gain_count);    // first MC histogram
        TH1D *mc1 = MC_Bi_214->ProjectionZ("Charge_Bi_214", eres_count, eres_count, gain_count, gain_count);    // second MC histogram
        TH1D *mc2 = MC_K_40->ProjectionZ("Charge_K_40", eres_count, eres_count, gain_count, gain_count);    // second MC histogram
        int bin = gain*1024/200000.0;

        for (int i = 0; i < bin; i++) {
          mc0->SetBinContent(i, 0);
          mc1->SetBinContent(i, 0);
          mc2->SetBinContent(i, 0);
          spectre_om->SetBinContent(i, 0);
        }

        roofitter(om, gain, eres, spectre_om, mc0, mc1, mc2, rootab, bin);
        logL = rootab[0];
        Chi2 = rootab[4];
        param1 = rootab[1];
        param2 = rootab[2];
        param3 = rootab[3];

        mc0->Scale(param1*spectre_om->Integral());
        om_flux_Tl = (mc0->Integral()*(10/9.0))/(1800*get_om_eff("Tl", om)*250000);
        std::cout << "om flux Tl = " << om_flux_Tl << '\n';
        mc1->Scale(param2*spectre_om->Integral());
        om_flux_Bi = (mc1->Integral()*(10/9.0))/(1800*get_om_eff("Bi", om)*250000);
        std::cout << "om flux Bi = " << om_flux_Bi << '\n';
        mc2->Scale(param3*spectre_om->Integral());
        om_flux_K = (mc2->Integral()*(10/9.0))/(1800*get_om_eff("K", om)*250000);
        std::cout << "om flux K = " << om_flux_K << '\n';


        Result_tree.Fill();

        delete mc0;
        delete mc1;
        delete mc2;
      }
    }
    delete spectre_om;

  }
  newfile->cd();
  Result_tree.Write();
  newfile->Close();
  std::cout << "good end" << '\n';
}

void distrib(string name) {
  TFile *eff_file = new TFile(Form("histo_fit/histo_%s.root", name.c_str()), "READ");
  // std::ofstream outFile("Best_Khi2.txt");
  int om_number =0;
  double Chi2 = 0;
  float gain = 0;
  float eres = 0;
  double om_flux_K = 0;
  double om_flux_Bi = 0;
  double om_flux_Tl = 0;
  double param1 = 0;
  double param2 = 0;
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
  double Gain = 0;
  double Tl = 10000;
  double Bi = 10000;
  double K = 10000;
  double par1 = 10000;
  double par2 = 10000;
  int om = 0;
  double test = 10;
  TFile *newfile = new TFile(Form("histo_fit/histo_%s_distrib.root", name.c_str()), "RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om);
  Result_tree.Branch("Chi2", &test);
  Result_tree.Branch("gain", &Gain);
  Result_tree.Branch("om_flux_Tl", &Tl);
  Result_tree.Branch("om_flux_Bi", &Bi);
  Result_tree.Branch("om_flux_K", &K);
  Result_tree.Branch("param1", &par1);
  Result_tree.Branch("param2", &par2);

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
        }
      }
    }
    if (om+1 == j ) {
      om = j;
      Result_tree.Fill();
    }
    Tl = 1000;
    Bi = 1000;
    K = 1000;
    par1 = 1000;
    par2 = 1000;
    test = 1000;
  }

  newfile->cd();
  Result_tree.Write();
  newfile->Close();
  std::cout << "distrib done" << '\n';
}

void histo_mystere(){

  TFile *histo_file_Tl = new TFile("Histo_simu_new/MC_simu_Tl_208_recOC_new_eres_53_gain_221_OC_MW8.root", "READ");
  histo_file_Tl->cd();
  TH3D* MC_Tl_208 = (TH3D*)histo_file_Tl->Get("MC_simu_Tl_208_recOC_new_ubc");

  TFile *histo_file_Bi = new TFile("Histo_simu_new/MC_simu_Bi_214_recOC_new_eres_53_gain_221_OC_MW8.root", "READ");
  histo_file_Bi->cd();
  TH3D* MC_Bi_214 = (TH3D*)histo_file_Bi->Get("MC_simu_Bi_214_recOC_new_ubc");

  TFile *histo_file_K = new TFile("Histo_simu_new/MC_simu_K_40_recOC_new_eres_53_gain_221_OC_MW8.root", "READ");
  histo_file_K->cd();
  TH3D* MC_K_40 = (TH3D*)histo_file_K->Get("MC_simu_K_40_recOC_new_ubc");

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

void name(/* arguments */) {
 TH3* MC = new TH3D("MC", "MC", 1024, 0, 200000,1,0,1,100,0,1);
 TRandom3 rando;
 MC->FillRandom("",1000);
 MC->Draw("colz");
}


void eff_om_prep(string name) {

  std::vector<int> *om_id = new std::vector<int>;
  std::vector<double> *energy = new std::vector<double>;
  float vertex_position[3];


  TFile *file = new TFile(Form("Histo_simu_new/Simu_%s_1G_OC.root", name.c_str()), "READ");

  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_id_new",1);
  tree->SetBranchAddress("om_id_new", &om_id);
  tree->SetBranchStatus("energy_u",1);
  tree->SetBranchAddress("energy_u", &energy);
  tree->SetBranchStatus("vertex_position",1);
  tree->SetBranchAddress("vertex_position", &vertex_position);

  TH1D energy_id("e_id", "e id ", 712, 0, 712);
  TH1D energy_id_cut("e_id_cut", "e id cut ", 712, 0, 712);

  for (int i = 0; i < tree->GetEntries(); i++) {
    memset(vertex_position, 0, 3*sizeof(float));
    tree->GetEntry(i);
    if (i % 100000 == 0) {
      std::cout << "entry = " << i << '\n';
    }
    for (size_t k = 0; k < om_id->size(); k++) {
      for (size_t j = 0; j < energy->size(); j++) {
        energy_id.Fill(om_id->at(k)-1);
        if (energy->at(j) > 1) {
          energy_id_cut.Fill(om_id->at(k)-1);
        }
      }
    }
  }

  TFile *newfile = new TFile(Form("eff_om/eff_prep_%s.root", name.c_str()), "RECREATE");
  newfile->cd();
  energy_id.Write();
  energy_id_cut.Write();
  newfile->Close();
}

void eff_om(string name) {

  double eff_tot =0;
  double eff_cut = 0;
  double eff_tot_error =0;
  double eff_cut_error = 0;
  int om = 0;

  TFile *file = new TFile(Form("eff_om/eff_prep_%s.root", name.c_str()), "READ");
  TH1D* eff_om_prep = (TH1D*)file->Get("e_id");
  TH1D* eff_om_prep_cut = (TH1D*)file->Get("e_id_cut");

  TH2F eff_tot_histo("eff", "efficacite mur ", 20, 0, 20, 13, 0, 13);
  TH2F eff_cut_histo("eff_cut", "efficacite mur cut", 20, 0, 20, 13, 0, 13);

  TFile *newfile = new TFile(Form("eff_om/eff_om_%s.root", name.c_str()), "RECREATE");
  TTree new_tree("new_tree","");
  new_tree.Branch("eff_tot", &eff_tot);
  new_tree.Branch("eff_tot_error", &eff_tot_error);
  new_tree.Branch("eff_cut", &eff_cut);
  new_tree.Branch("eff_cut_error", &eff_cut_error);
  new_tree.Branch("om", &om);

  // for (int i = 0; i < eff_om_prep->GetMaximumBin(); i++) {
  for (int i = 0; i < 712; i++) {
    om = i;
    int om_col = (i % 13);
    int om_row = (i / 13);

    eff_tot = (eff_om_prep->GetBinContent(i))/1.0e5;
    eff_tot_error = (sqrt(eff_om_prep->GetBinContent(i))/1.0e5);
    eff_cut = (eff_om_prep_cut->GetBinContent(i))/1.0e5;
    eff_cut_error = (sqrt(eff_om_prep_cut->GetBinContent(i))/1.0e5);

    eff_tot_histo.SetBinContent( om_row+1, om_col+1, eff_tot);
    eff_cut_histo.SetBinContent( om_row+1, om_col+1, eff_cut);

    new_tree.Fill();
  }

  file->Close();

  newfile->cd();

  new_tree.Write();
  eff_tot_histo.Write();
  eff_cut_histo.Write();

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
    eff_tot = (energy_id.GetBinContent(i))/5.0e5;
    eff_tot_error = (sqrt(energy_id.GetBinContent(i))/5.0e5);


    new_tree.Fill();
  }

  newfile->cd();
  new_tree.Write();
  newfile->Close();
  file->Close();
}

void fusion() {


  double om_flux_Tl, om_flux_Bi, om_flux_K, par1, par2, Tl, Bi, K, p1, p2, test, Gain, Chi2, gain;
  int om, om_id;


  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om);
  Result_tree.Branch("Chi2", &test);
  Result_tree.Branch("gain", &Gain);
  Result_tree.Branch("om_flux_Tl", &om_flux_Tl);
  Result_tree.Branch("om_flux_Bi", &om_flux_Bi);
  Result_tree.Branch("om_flux_K", &om_flux_K);
  Result_tree.Branch("param1", &par1);
  Result_tree.Branch("param2", &par2);

  TFile *file1 = new TFile("histo_fit/histo_it2_distrib.root", "READ");
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

  for (size_t i = 0; i < 260; i++) {
    tree->GetEntry(i);

    om = om_id;
    om_flux_Tl = Tl;
    om_flux_Bi = Bi;
    om_flux_K = K;
    par1 = p1;
    par2 = p2;
    Gain = gain;
    test = Chi2;
    Result_tree.Fill();

  }

  TFile *file2 = new TFile("histo_fit/histo_fr2_distrib.root", "READ");
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

  for (size_t i = 260; i < 520; i++) {
    tree2->GetEntry(i);

    om = om_id;
    om_flux_Tl = Tl;
    om_flux_Bi = Bi;
    om_flux_K = K;
    par1 = p1;
    par2 = p2;
    Gain = gain;
    test = Chi2;
    Result_tree.Fill();

  }

  TFile *file3 = new TFile("histo_fit/histo_xw2_distrib.root", "READ");
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

  for (size_t i = 520; i < 648; i++) {
    tree3->GetEntry(i);

    om = om_id;
    om_flux_Tl = Tl;
    om_flux_Bi = Bi;
    om_flux_K = K;
    par1 = p1;
    par2 = p2;
    Gain = gain;
    test = Chi2;
    Result_tree.Fill();

  }

  TFile *file4 = new TFile("histo_fit/histo_gv2_distrib.root", "READ");
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

  for (size_t i = 648; i < 712; i++) {
    tree4->GetEntry(i);

    om = om_id;
    om_flux_Tl = Tl;
    om_flux_Bi = Bi;
    om_flux_K = K;
    par1 = p1;
    par2 = p2;
    Gain = gain;
    test = Chi2;
    Result_tree.Fill();

  }

  TFile *newfile = new TFile("histo_fit/test2.root", "RECREATE");
  newfile->cd();
  Result_tree.Write();
  newfile->Close();
  std::cout << "ok " << '\n';

}

int main(int argc, char const *argv[]){
  for(int i = 0; i<argc; i++){
    if (std::string(argv[i]) == "-XW" || std::string(argv[i]) =="-MW" || std::string(argv[i]) =="-GV"){
      Fit_Gain_Simu(argv[i]);
      if(std::string(argv[i]) == "--distrib" || std::string(argv[i]) =="-d"){
        distrib("roofit");
      }
    }
    else{
      Fit_Gain_Simu("FULL");
    }
  }


  return 0;
}
