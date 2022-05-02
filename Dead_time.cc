#include "Dead_time.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include <iostream>
#include <memory>
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
#include <TSystem.h>
using namespace std;

void dead_time(string run){

  double hic = 0;
  int om = 0;
  uint16_t lt_counter = 0;
  uint32_t lt_time = 0;
  ULong64_t tdc = 0;
  double deadt = 0;
  int ch_lt = 0;
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("om_number", &om);
  Result_tree.Branch("deadt", &deadt);
  Result_tree.Branch("hic", &hic);

  TFile *file = new TFile(Form("histo_kolmo/histo_donee/%s.root", run.c_str()), "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_number",1);
  tree->SetBranchAddress("om_number", &om);
  tree->SetBranchStatus("lt_trigger_counter",1);
  tree->SetBranchAddress("lt_trigger_counter", &lt_counter);
  tree->SetBranchStatus("lt_trigger_time",1);
  tree->SetBranchAddress("lt_trigger_time", &lt_time);
  tree->SetBranchStatus("tdc",1);
  tree->SetBranchAddress("tdc", &tdc);
  tree->SetBranchStatus("ch_lt",1);
  tree->SetBranchAddress("ch_lt", &ch_lt);
  double compteur_time[717];
  memset(compteur_time, 0, 717*sizeof(double));
  double compteur_count[717];
  memset(compteur_count, 0, 717*sizeof(double));
  double compteur_chlt[717];
  memset(compteur_chlt, 0, 717*sizeof(double));
  double run_time = 0;
  double cut_tdc = 0;
  double tdc0;

  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    if (i == 0) {
      tdc0 = tdc;
    }
    if (i == tree->GetEntries() - 1) {
      run_time = tdc - tdc0;
    }
    if (i%10000 == 0 || i%9999 == 0) {
      // std::cout << "compteur_time = " << compteur_time[om] << " and lt_time = " << lt_time<< '\n';
    }
    compteur_time[om] += (double)lt_time;
    compteur_count[om] += lt_counter;

    if (ch_lt > 0 ) {
      compteur_chlt[om]+= 1;
      // std::cout << "lt _counter = " << ch_lt << '\n';
    }
    if (tdc-cut_tdc > 3e8) {
      run_time+=(cut_tdc-tdc);
    }
    cut_tdc = tdc;
  }

  for (int i = 0; i < 260; i++) {
    std::cout <<"om : " << i << " ch_lt = " << compteur_chlt[i] << " and run time = " << run_time << " and count = " << compteur_count[i] << " and time = " << compteur_time[i] << '\n';
    deadt = (compteur_chlt[i]/run_time)/(compteur_count[i]/compteur_time[i]);
    hic = compteur_chlt[i]/compteur_count[i];
    if (compteur_count[i] == 0) {
      deadt = 0;
      hic =0;
    }
    std::cout << "deadt = " << deadt << " and hic = " << hic << '\n';
    Result_tree.Fill();
  }

  TFile* newfile = new TFile(Form("deadt/deadt_cut_%s.root", run.c_str()),"RECREATE");
  newfile->cd();
  Result_tree.Write();
  newfile->Close();
}
