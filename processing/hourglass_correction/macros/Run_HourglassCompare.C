#include "../../../FileManagement.h"

int Run_HourglassCompare(
  int run_index = 0,
  std::string sim_file =  "/direct/phenix+spin2/beaumim/vernierScans/plots/359711_h0.1_v0_beta85_xing8e-05_HourglassSimulation.root",
  std::string data_file = "/direct/phenix+spin2/beaumim/vernierScans/plots/359711_HourglassData.root" ,
  std::string sim_hist_name = "zvertex_simulation" ,
  std::string data_hist_name = "zdc_zvtx_step_0"  ,
  const std::string lib = "libVernierHourglass.so"
){
  TFile* data = new TFile(data_file.c_str(),"READ");
  TFile* simu = new TFile(sim_file.c_str(),"READ");

  TH1F* data_hist = (TH1F*) data->Get(data_hist_name.c_str());
  TH1F* sim_hist  = (TH1F*) simu->Get(sim_hist_name.c_str());

  TCanvas* compare = new TCanvas("compare","Compare",800,600);
  data_hist->Draw();
  sim_hist->SetLineColor(kRed);
  sim_hist->SetLineWidth(2);
  sim_hist->Draw("same");

  return 0;
}
