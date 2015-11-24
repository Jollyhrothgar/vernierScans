#include <string>
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
int make_pretty_graphs(
  const std::string run_number,
  const std::string out_dir
){
  std::cout << run_number << std::endl;
  TFile* f = new TFile("wcm_dists.root","READ");
  bname = "bwcm_zprofile_"+run_number+"_space"; 
  yname = "ywcm_zprofile_"+run_number+"_space";

  TGraph* blue = (TGraph*) f->Get(bname.c_str());
  TGraph* yell = (TGraph*) f->Get(yname.c_str());

  blue->SetLineWidth(2);
  yell->SetLineWidth(2);

  blue->SetLineColor(kBlue-2);
  yell->SetLineColor(kYellow-2);
  
  std::string out_file = out_dir + run_number + "_wcm_zprofile.pdf";
  
  TCanvas * c = new TCanvas("can","compare",1200,800);
  TPad* p1 = new TPad("p1","left pad",0,0,0.67,1.);
  TPad* p2 = new TPad("p2","right top",0.67,0.5,1.,1.);
  TPad* p3 = new TPad("p3","right bottom",0.67,0,1.,0.5);

  p2->cd();
  blue->Draw("AL");
  p3->cd();
  yell->Draw("AL");
  p1->cd();
  blue->Draw("AL");
  yell->Draw("L");
  c->cd();
  p1->Draw();
  p2->Draw();
  p3->Draw();
  c->SaveAs(out_file.c_str());
  return 0;
}
