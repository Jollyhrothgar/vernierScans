#include "../../../FileManagement.h"

int Run_HourglassCompare(
  std::string canvas_name  = "config_and_vertex_compare",
  std::string norm_canvas_name = "norm_config_and_vertex_compare",
  std::string compare_file = "test_compare.root" 
){

  TFile* f = new TFile(compare_file.c_str());
  TCanvas* c1 = (TCanvas*)f->Get(canvas_name.c_str());
  TCanvas* c2 = (TCanvas*)f->Get(norm_canvas_name.c_str());
  c1->Draw();
  c2->Draw();

  return 0;
}
