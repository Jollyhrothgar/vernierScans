#include "../../FileManagement.h"

// "359711" index 0
// "360879" index 1
// "362492" index 2
// "364636" index 3
// "365866" index 4
// "366605" index 5
// "367138" index 6
Double_t par_b[7][9];
Double_t par_y[7][9];

Double_t g1b_range[7][2] = 
{
  {-180.,-90.},  // index 0
  {-180.,-100.}, // index 1
  {-150.,-90.},  // index 2
  {-180.,-90.},
  {-180.,-90.},
  {-180.,-90.},
  {-180.,-90.}
};

Double_t g2b_range[7][2] = 
{
  {-40.,30.},  // index 0
  {-40.,30.},  // index 1
  {-30.,30.},  // index 2
  {-40.,30.},
  {-40.,30.},
  {-40.,30.},
  {-40.,30.}
};

Double_t g3b_range[7][2] = 
{
  {120.,170.}, // index 0
  {110.,175.}, // index 1
  {80.,180.},  // index 2
  {120.,170.},
  {121.,170.},
  {120.,170.},
  {120.,170.}
};

Double_t g1y_range[7][2] = 
{
  {-180,-90.},   // index 0
  {-180.,-120.}, // index 1
  {-160,-90.},   // index 2
  {-180,-90.},
  {-180,-90.},
  {-180,-90.},
  {-180,-90.}
};

Double_t g2y_range[7][2] = 
{
  {-40.,-30.}, // index 0
  {-30.,30.},  // index 1
  {-30.,-30.}, // index 2
  {-40.,-30.},
  {-40.,-30.},
  {-40.,-30.},
  {-40.,-30.}
};

Double_t g3y_range[7][2] = 
{
  {120.,170.}, // index 0
  {110.,170.}, // index 1
  {90.,160.},  // index 2
  {120.,170.},
  {120.,170.},
  {120.,170.},
  {120.,170.}
};

void fit_blue_yellow(int index = 0){
  const int run_index = index;
  std::string in_tfile_name = "wcm_dists.root";
  std::string out_tfile_name = "fits.root";
  
  TFile* f = new TFile(in_tfile_name.c_str(),"READ");
  std::string g_blue_name = "bwcm_zprofile_"+run_number[run_index]+"_max_centered";
  std::string g_yell_name = "ywcm_zprofile_"+run_number[run_index]+"_max_centered";
  TGraph* g_blue = (TGraph*)f->Get(g_blue_name.c_str());
  TGraph* g_yell = (TGraph*)f->Get(g_yell_name.c_str());

  TF1* g1_b = new TF1("g1_b","gaus", g1b_range[run_index][0], g1b_range[run_index][1]);
  TF1* g2_b = new TF1("g2_b","gaus", g2b_range[run_index][0], g2b_range[run_index][1]);
  TF1* g3_b = new TF1("g3_b","gaus", g3b_range[run_index][0], g3b_range[run_index][1]);
  
  TF1* g1_y = new TF1("g1_y","gaus", g1y_range[run_index][0], g1y_range[run_index][1]);
  TF1* g2_y = new TF1("g2_y","gaus", g2y_range[run_index][0], g2y_range[run_index][1]);
  TF1* g3_y = new TF1("g3_y","gaus", g3y_range[run_index][0], g3y_range[run_index][1]);

  TF1* simple_gaus_blue = new TF1("simple_gaus_blue","gaus",-1500,1500);
  TF1* simple_gaus_yell = new TF1("simple_gaus_yell","gaus",-1500,1500);
  g_blue->Fit(simple_gaus_blue);
  g_yell->Fit(simple_gaus_yell);

  g1_b->SetLineColor(kBlue);
  g2_b->SetLineColor(kGreen);
  g3_b->SetLineColor(kOrange);
  
  g1_y->SetLineColor(kBlue);
  g2_y->SetLineColor(kGreen);
  g3_y->SetLineColor(kOrange);

  TF1* gb_tot = new TF1("blue_zprofile","gaus(0)+gaus(3)+gaus(6)",g1b_range[run_index][0],g3b_range[run_index][1]);
  TF1* gy_tot = new TF1("yellow_zprofile","gaus(0)+gaus(3)+gaus(6)",g1y_range[run_index][0],g3y_range[run_index][1]);

  gb_tot->SetLineColor(kRed);

  g_blue->Fit(g1_b,"R");
  g_blue->Fit(g2_b,"R+");
  g_blue->Fit(g3_b,"R+");

  g_yell->Fit(g1_y,"R");
  g_yell->Fit(g2_y,"R+");
  g_yell->Fit(g3_y,"R+");

  g1_b->GetParameters(&par_b[run_index][0]);
  g2_b->GetParameters(&par_b[run_index][3]);
  g2_b->GetParameters(&par_b[run_index][6]);
  
  g1_y->GetParameters(&par_y[run_index][0]);
  g2_y->GetParameters(&par_y[run_index][3]);
  g2_y->GetParameters(&par_y[run_index][6]);

  gb_tot->SetParameters(par_b[run_index]);
  gy_tot->SetParameters(par_y[run_index]);

  g_blue->Fit(gb_tot,"R");
  g_yell->Fit(gy_tot,"R");

  std::string out_file_name = "zprofile_fit_run_"+run_number[run_index]+".root";
  TFile* f_out = new TFile(out_file_name.c_str(),"RECREATE");
  f_out->cd();
  g_blue->Write();
  g1_b->Write();
  g2_b->Write();
  g3_b->Write();
  gb_tot->Write();
  g_yell->Write();
  g1_y->Write();
  g2_y->Write();
  g3_y->Write();
  gy_tot->Write();
  simple_gaus_blue->Write();
  simple_gaus_yell->Write();
  f_out->Write();
}
