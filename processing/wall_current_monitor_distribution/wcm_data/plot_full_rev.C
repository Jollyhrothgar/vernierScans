#include<fstream>
#include<iostream>
#include "TFile.h"
#include "TGraph.h"
void plot_full_rev(){
  TFile* f = new TFile("rev.root","RECREATE");
  TGraph* g = new TGraph();
  g->SetName("full_revolution");
  g->SetTitle("Full Revolution;time;intensity");
  g->SetLineWidth(2);
  g->SetLineColor(kRed);
  std::ifstream in_file("blueWCM_359711_16444_full_revolution.dat");
  double time_ns = 0.;
  double intensity = 0.;
  while(in_file >> time_ns >> intensity){
    g->SetPoint(g->GetN(),time_ns,intensity);
  }
  g->Write();
  f->Write();
  delete g;
  f->Close();
  delete f;
}
