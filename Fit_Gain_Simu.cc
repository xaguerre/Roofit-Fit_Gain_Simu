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

TH1D* spectre_charge_it(int om_number )
{
  TFile *file = new TFile("histo_kolmo/histo_donee/histo_charge_amplitude_422.root", "READ");
  gROOT->cd();
  TH2F* charge = (TH2F*)file->Get("charge");

  TH1D* spectre_charge = charge->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
  file->Close();
  return spectre_charge;
}

TH1D* spectre_charge_XW(int om_number )
{
  TFile *file = new TFile("histo_kolmo/histo_donee/histo_charge_amplitude_energie_449.root", "READ");
  gROOT->cd();
  TH2F* charge = (TH2F*)file->Get("histo_pm_charge");

  TH1D* spectre_charge = charge->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
  file->Close();
  return spectre_charge;
}

TH1D* spectre_charge_fr(int om_number )
{
  TFile *file = new TFile("histo_kolmo/histo_donee/histo_charge_amplitude_energie_435.root", "READ");
  gROOT->cd();
  TH2F* charge = (TH2F*)file->Get("histo_pm_charge");

  TH1D* spectre_charge = charge->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
  file->Close();
  return spectre_charge;
}

TH1D* spectre_chooser(int om)
{
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

TH3D* MC(int value, string name) {
  TH3D* MC = NULL;
  if (value == 1) {
    TFile *histo_file_Tl = new TFile(Form("Histo_simu_new/MC_simu_Tl_208_recOC_new_eres_53_gain_221_OC_%s.root", name.c_str()), "READ");
    histo_file_Tl->cd();
    MC = (TH3D*)histo_file_Tl->Get("MC_simu_Tl_208_recOC_new_ubc");
  }
  if (value == 2) {
    TFile *histo_file_Bi = new TFile("Histo_simu_new/MC_simu_Bi_214_recOC_new_eres_53_gain_221_OC_MW8.root", "READ");
    histo_file_Bi->cd();
    MC = (TH3D*)histo_file_Bi->Get("MC_simu_Bi_214_recOC_new_ubc");
  }
  if (value == 3) {
    TFile *histo_file_K = new TFile("Histo_simu_new/MC_simu_K_40_recOC_new_eres_53_gain_221_OC_MW8.root", "READ");
    histo_file_K->cd();
    MC = (TH3D*)histo_file_K->Get("MC_simu_K_40_recOC_new_ubc");
  }
  return MC;
}

TH3D* MC_chooser(int om, int comp){
  TH3D* MC_chooser = NULL;
  if (om <520 && om%13 != 12 && om%13 != 0){
    MC_chooser = MC(comp, "MW8");
  }
  else if (om <520 && (om%13 == 12 || om%13 == 0)) {
    MC_chooser = MC(comp, "MW5");
  }
  else if (om > 519 && om < 648) {
    MC_chooser = MC(comp, "XW");
  }
  else if (om > 647) {
    MC_chooser = MC(comp, "GV");
  }
  return MC_chooser;
}

double* om_gain_fit(int om)
{
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

double* roofitter(int om, double gain, double eres, TH1D* spectre_om, TH1D* mc0, TH1D* mc1, TH1D* mc2, double *rootab)
{
  const float start_x = spectre_om->GetXaxis()->GetBinUpEdge(105);

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

  RooRealVar RooTl("RooTl", "Tl", 0.2, 0, 1);
  RooRealVar RooBi("RooBi", "Bi", 0.5, 0, 1);
  RooRealVar RooK("RooK", "K", 0.3, 0, 1);

  RooRealSumPdf sum_simu("sum_simu", "sum_simu",
  RooArgList(p_ph_Tl, p_ph_Bi, p_ph_K),
  RooArgList(RooTl,RooBi,RooK),
  false);

  RooFitResult* result1 = sum_simu.fitTo(spectre_data, PrintLevel(-1), SumW2Error(true), Range(start_x,65000), Save(), IntegrateBins(-1), Minimizer("Minuit"));
  sum_simu.setStringAttribute("fitrange", nullptr);
  TCanvas* can = new TCanvas("can", "", 1500, 600);
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
  sum_simu.plotOn(frame, Components(p_ph_Bi), LineColor(kAzure), Name("Bi_curve"));
  sum_simu.plotOn(frame, Components(p_ph_K), LineColor(kBlack), Name("K_curve"));
  sum_simu.paramOn(frame);

  RooCurve * Tl_curve = frame->getCurve("Tl_curve");
  double Tl_int = Tl_curve->average(start_x,65000);
  RooCurve * Bi_curve = frame->getCurve("Bi_curve");
  double Bi_int = Bi_curve->average(start_x,65000);
  RooCurve * K_curve = frame->getCurve("K_curve");
  double K_int = K_curve->average(start_x,65000);
  RooCurve * sum_curve = frame->getCurve("sum_curve");
  double int_data = sum_curve->average(start_x,65000);

  frame->GetYaxis()->UnZoom();
  frame ->Draw();
  can->SetLogy();
  can->SaveAs(Form("Best_fit/best_fit_om_%d_eres_%.2f_gain_%.0f.png", om, eres, gain));
  // result1->Print("v");
  // result1->Print();
  rootab[0] = result1->minNll();
  rootab[1] = Tl_int/int_data;
  rootab[2] = Bi_int/int_data;
  rootab[3] = K_int/int_data;

  std::cout << "Tl  =" << rootab[1]<< '\n';
  std::cout << "Bi  =" << rootab[2]<< '\n';
  std::cout << "K  =" << rootab[3]<< '\n';
  std::cout << "Tot  =" << rootab[3]+rootab[1]+rootab[2]<< '\n';
  //
  // delete can;
  // delete frame;
  // delete result1;
  return rootab;
}

double get_om_eff(string compo, int om) {
  double om_eff = 0;
  TFile *newfile = new TFile(Form("eff_om/eff_om_%s.root", compo.c_str()), "READ");

  TTree* tree = (TTree*)newfile->Get("new_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("eff_tot",1);
  tree->SetBranchAddress("eff_tot", &om_eff);

  tree->GetEntry(om);
  newfile->Close();
  return om_eff;
}

void Fit_Gain_Simu() {
  TH1::SetDefaultSumw2();

  double logL, Chi2, param1, param2, param3 = 0;
  float gain, eres = 0;
  int om_number = 0;
  double mean_erf, sigma_erf = 0;
  double* tab = new double[7];
  double om_flux_Tl, om_flux_Bi, om_flux_K, integrale_tot_Tl, integrale_droite_Tl, integrale_tot_Bi, integrale_droite_Bi, integrale_tot_K, integrale_droite_K = 0;

  TFile *newfile = new TFile("histo_fit/histo_roofit.root", "RECREATE");
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

  std::ifstream  charge("/home/aguerre/Bureau/Thèse/Fit_Gain_Simu/gain_data/Resultat_mean_energie_charge_total.txt");
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

  for (int om = 0; om < 13; om++)
  {
    int lim = 140;
    TH3D* MC_Tl_208 = MC_chooser(om, 1);
    TH3D* MC_Bi_214 = MC_chooser(om, 2);
    TH3D* MC_K_40 = MC_chooser(om, 3);

    om_number = om;
    TH1D* spectre_om = NULL;
    spectre_om = spectre_chooser(om);

    if ((spectre_om->GetEntries() < 100) || (spectre_om->GetMean(1) < 15000)){
      delete spectre_om;
      om++;
      spectre_om = spectre_chooser(om);
    }
    if (charge_valeur_fit[om] != -1){
      if (om < 648){
        tab = om_gain_fit(om);
        if (tab[0] == 0) {
          mean_erf = charge_valeur_fit[om];
          sigma_erf = 0;
        }
        else{
          sigma_erf = 1/tab[1];
          mean_erf = 1/tab[0];
        }
      }
    }

    for (int bin =1; bin < lim; bin++) {
      spectre_om->SetBinContent(bin, 0);
    }

    for (int eres_count = 13; eres_count < 14; eres_count++) {
      for (int gain_count = 0; gain_count <150; gain_count++) {
        std::cout << (gain_bin_min + gain_bin_width*(gain_count-1))<< "   et    lim_inf = " << 1/charge_valeur_fit[om]*0.9 << "   sup  ="  << 1/charge_valeur_fit[om]*1.1 << '\n';
        if ((1/(gain_bin_min + gain_bin_width*(gain_count-1)) > mean_erf*0.85) && (1/(gain_bin_min + gain_bin_width*(gain_count-1))<mean_erf*1.15)){
          gain = (gain_bin_min + gain_bin_width*(gain_count-1));
          eres = eres_bin_min + eres_bin_width*(eres_count-1);
          std::cout << "gain_count = " << gain_count << " and gain = " << gain << '\n';
          TH1D *mc0 = MC_Tl_208->ProjectionZ("Charge_Tl_208", eres_count, eres_count, gain_count, gain_count);    // first MC histogram
          TH1D *mc1 = MC_Bi_214->ProjectionZ("Charge_Bi_214", eres_count, eres_count, gain_count, gain_count);    // second MC histogram
          TH1D *mc2 = MC_K_40->ProjectionZ("Charge_K_40", eres_count, eres_count, gain_count, gain_count);    // second MC histogram

          integrale_tot_Tl = mc0->Integral(0, 1024);
          integrale_tot_Bi = mc1->Integral(0, 1024);
          integrale_tot_K = mc2->Integral(0, 1024);

          for (size_t i = 0; i < lim; i++) {
            mc0->SetBinContent(i, 0);
            mc1->SetBinContent(i, 0);
            mc2->SetBinContent(i, 0);
          }
          integrale_droite_Tl = mc0->Integral(lim+1, 1024);
          integrale_droite_Bi = mc1->Integral(lim+1, 1024);
          integrale_droite_K = mc2->Integral(lim+1, 1024);

          roofitter(om, gain, eres, spectre_om, mc0, mc1, mc2, rootab);
          logL = rootab[0];
          Chi2 = rootab[4];
          param1 = rootab[1];
          param2 = rootab[2];
          param3 = rootab[3];

          mc0->Scale(param1*spectre_om->Integral());
          om_flux_Tl = mc0->Integral()/1800*get_om_eff("Tl_208", om)/655;
          std::cout << "om flux Tl = " << om_flux_Tl << '\n';
          mc1->Scale(param2*spectre_om->Integral());
          om_flux_Bi = mc1->Integral()/1800*get_om_eff("Bi_214", om)/655;
          std::cout << "om flux Bi = " << om_flux_Bi << '\n';
          mc2->Scale(param3*spectre_om->Integral());
          om_flux_K = mc2->Integral()/1800*get_om_eff("K_40", om)/655;
          std::cout << "om flux K = " << om_flux_K << '\n';

          Result_tree.Fill();

          delete mc0;
          delete mc1;
          delete mc2;
        }
      }
    }
    delete spectre_om;
    delete MC_Tl_208;
    delete MC_Bi_214;
    delete MC_K_40;
  }
  newfile->cd();
  Result_tree.Write();
  newfile->Close();
}

void distrib(string name) {
  TFile *eff_file = new TFile(Form("histo_fit/histo_%s.root", name.c_str()), "READ");
  std::ofstream outFile("Best_Khi2.txt");
  int om_number =0;
  double Chi2 = 0;
  float gain = 0;
  float eres = 0;
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
  double n_entry = 0;

  double test = 10;
  if (name.compare("MW8") == 0){
    TH1D* dist = new TH1D ("distribution Chi2","distribution Chi2 OM MW 8p",100, 0, 2);
    dist->GetXaxis()->SetTitle("Khi2");
    for (int j = 1; j < 520; j++) {
      if (((j%13)!=0) && ((j%13)!=12)){
        for (double i = 0; i < eff_tree->GetEntries(); i++) {
          eff_tree->GetEntry(i);
          if (om_number == j) {
            if (Chi2 < test) {
              test = Chi2;
              n_entry = i;
            }
          }
        }
      }
      dist->Fill(test);
      test =10;
      if (((j%13)!=0) && ((j%13)!=12)){
        eff_tree->GetEntry(n_entry);
        outFile << om_number << "\t" << gain << "\t" << eres << endl;
      }
    }
    dist->Draw();
  }
  else if (name.compare("MW5") == 0){
    TH1D* dist = new TH1D ("distribution Chi2","distribution Chi2 OM 5p MW",100, 0, 2);
    dist->GetXaxis()->SetTitle("Khi2");
    for (int j = 0; j < 520; j++) {
      if (((j%13)==0) || ((j%13)==12)){
        for (double i = 0; i < eff_tree->GetEntries(); i++) {
          eff_tree->GetEntry(i);
          if (om_number == j) {
            if (Chi2 < test) {
              test = Chi2;
              n_entry = i;
            }
          }
        }
      }
      dist->Fill(test);
      cout << j << endl;
      test = 10;
      if (((j%13)==0) || ((j%13)==12)){
        eff_tree->GetEntry(n_entry);
        outFile << om_number << "\t" << gain << "\t" << eres << endl;
      }
    }
    dist->Draw();
  }
  else if (name.compare("XW") == 0){
    TH1D* dist = new TH1D ("distribution Chi2","distribution Chi2 OM XW",100, 0, 2);
    dist->GetXaxis()->SetTitle("Khi2");
    for (int j = 520; j < 648; j++) {
      for (double i = 0; i < eff_tree->GetEntries(); i++) {
        eff_tree->GetEntry(i);
        if (om_number == j) {
          if (Chi2 < test) {
            test = Chi2;
            n_entry = i;
          }
        }
      }
      dist->Fill(test);
      test =10;
      eff_tree->GetEntry(n_entry);
      outFile << om_number << "\t" << gain << "\t" << eres << endl;
    }
    dist->Draw();
  }
  else if (name.compare("GV") == 0){
    TH1D* dist = new TH1D ("distribution Chi2","distribution Chi2 OM GV",100, 0, 11);
    dist->GetXaxis()->SetTitle("Khi2");
    for (int j = 648; j < 712; j++) {
      for (double i = 0; i < eff_tree->GetEntries(); i++) {
        eff_tree->GetEntry(i);
        if (om_number == j) {
          if (Chi2 < test) {
            test = Chi2;
            n_entry = i;
          }
        }
      }
      dist->Fill(test);
      test =10;
      eff_tree->GetEntry(n_entry);
      outFile << om_number << "\t" << gain << "\t" << eres << endl;
    }
    dist->Draw();
  }
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

  TCanvas* can = new TCanvas("can", "", 1500, 600);

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

    if (vertex_position[1] > 4999) {
      for (int k = 0; k < om_id->size(); k++) {
        for (int j = 0; j < energy->size(); j++) {
          energy_id.Fill(om_id->at(k)-1);
          if (energy->at(j) > 1) {
            energy_id_cut.Fill(om_id->at(k)-1);
          }
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

  double lim =0;
  double eff_tot =0;
  double eff_cut = 0;
  double eff_tot_error =0;
  double eff_cut_error = 0;
  int om = 0;
  std::vector<int> *om_id = new std::vector<int>;
  std::vector<double> *energy = new std::vector<double>;

  TFile *file = new TFile(Form("eff_om/eff_prep_%s.root", name.c_str()), "READ");
  TH1D* eff_om_prep = (TH1D*)file->Get("e_id");
  TH1D* eff_om_prep_cut = (TH1D*)file->Get("e_id_cut");

  TH2F eff_tot_histo("eff", "efficacite mur ", 20, 0, 20, 13, 0, 13);
  TH2F eff_cut_histo("eff_cut", "efficacite mur cut", 20, 0, 20, 13, 0, 13);

  TFile *newfile = new TFile(Form("eff_om/eff_om_%s.root", name.c_str()), "RECREATE");
  TTree new_tree("new_tree","");
  new_tree.Branch("eff_tot", &eff_tot);
  new_tree.Branch("eff_cut", &eff_cut);
  new_tree.Branch("om", &om);

  // for (int i = 0; i < eff_om_prep->GetMaximumBin(); i++) {
  for (int i = 0; i < 712; i++) {
    om = i;
    int om_col = (i % 13);
    int om_row = (i / 13);

    eff_tot = (eff_om_prep->GetBinContent(i))/1.0e6;
    eff_tot_error = (sqrt(eff_om_prep->GetBinContent(i))/1.0e6);
    eff_cut = (eff_om_prep_cut->GetBinContent(i))/1.0e6;
    eff_cut_error = (sqrt(eff_om_prep_cut->GetBinContent(i))/1.0e6);

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

int main(int argc, char const *argv[])
{
  Fit_Gain_Simu();
  return 0;
}
