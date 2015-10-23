#include "../../../FileManagement.h"
#include <iostream>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TH1F.h"

int Run_HourglassCompare(
  std::string canvas_name  = "config_and_vertex_compare",
  std::string norm_canvas_name = "norm_config_and_vertex_compare",
  std::string compare_file = "test_compare.root" 
){

  TFile* f = new TFile(compare_file.c_str());
  TCanvas* c1 = (TCanvas*)f->Get(canvas_name.c_str());
  TCanvas* c2 = (TCanvas*)f->Get(norm_canvas_name.c_str());
  c2->Draw();
  c1->Draw();

  // These must be binned the same!!
  TH1F* sim = (TH1F*)f->Get("zdc_zvertex_sim");
  TH1F* dat = (TH1F*)f->Get("zdc_zvertex_data");

  // Get Chi2Test
  std::cout << dat->Chi2Test(sim,"UW P") << std::endl;
  
  // Generate average residual
  double average_res = 0;
  for(int i = 0; i < dat->GetNbinsX(); i++) {
    average_res +=  pow(dat->GetBinContent(i) - sim->GetBinContent(i),2.0);
  }
  std::cout << "Average residual: " << average_res/(dat->GetNbinsX()) << std::endl;

  return 0;
}
