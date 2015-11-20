#include <string>
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
int Run_ParamExplore(
  const std::string param_1,
  const std::string param_2,
  const std::string param_name
) {
  std::string dist_name = "zdc_zvertex_sim";
  TFile* f1 = new TFile(param_1.c_str(),"READ");
  TFile* f2 = new TFile(param_2.c_str(),"READ");

  TH1F* h1 = (TH1F*) f1->Get(dist_name.c_str());
  TH1F* h2 = (TH1F*) f2->Get(dist_name.c_str());

  h1->SetLineColor(kGreen-2);
  h2->SetLineColor(kPink-2);

  std::string out_file = param_name + ".pdf";

  TCanvas * c = new TCanvas("can","compare",1200,800);
  TPad* p1 = new TPad("p1","left pad",0,0,0.67,1.);
  TPad* p2 = new TPad("p2","right top",0.67,0.5,1.,1.);
  TPad* p3 = new TPad("p3","right bottom",0.67,0,1.,0.5);

  p2->cd();
  h1->Draw();
  p3->cd();
  h2->Draw();
  p1->cd();
  h1->Draw();
  h2->Draw("same");
  c->cd();
  p1->Draw();
  p2->Draw();
  p3->Draw();

  c->SaveAs(out_file.c_str());
 //delete h2;
 //delete h1;
 //delete c;
 //delete f1;
 //delete f2;
  return 0;
}
  
