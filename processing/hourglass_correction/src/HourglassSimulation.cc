#include "HourglassSimulation.h"
#include "HourglassConfiguration.h"

#include <cstdlib>
#include <fstream>
#include <cmath>
#include <ctime>
#include <math.h>
#include <chrono> 
#include <map>
#include <vector> 
#include <sstream>
#include <iostream>
#include <sstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <algorithm>

#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraph.h"

#include "TError.h"

#define PI 3.14159265

// Local function definition to avoid dealing with ROOT's immense stupidity
std::chrono::milliseconds GetTime() {
  std::chrono::milliseconds the_time;
  the_time = std::chrono::duration_cast < std::chrono::milliseconds >(std::chrono::system_clock::now().time_since_epoch());
  return the_time;
}

HourglassSimulation::HourglassSimulation() {
  std::stringstream ss;
  ss << "HourglassSimulation_" << std::hex << this;
  this_name = ss.str();
  std::cout << "Instantiating " << this_name << std::endl;
  override_save_file_ = false;
  how_many_things = 0; 
  amaresh_model_run = false;
  new_model_run = false;
  new_fit_model_run = false;
  simple_gaus_model_run = false;
  first_init = false;
  zdc_zvertex_sim = NULL;
  zdc_zvertex_dat = NULL;
  z_bunch_profile_blue_ = NULL;
  z_bunch_profile_yell_ = NULL;
  zvertex_comparison_canvas = NULL;
  simulation_config_canvas = NULL;
  config_and_vertex_compare = NULL;
  amaresh_compare = NULL;
  amaresh_z_blue = NULL;
  amaresh_z_yell = NULL;
  gaus_compare = NULL;
  gaus_z_blue = NULL;
  gaus_z_yell = NULL;
  z_lookup_diff_blue = NULL;
  z_lookup_diff_yell = NULL;
  new_model_z_blue[0] = NULL;
  new_model_z_yell[0] = NULL;
  new_model_z_blue[1] = NULL;
  new_model_z_yell[1] = NULL;
  new_model_z_blue[2] = NULL;
  new_model_z_yell[2] = NULL;
  how_many_runs = 0;
  save_all = false;
  save_directory_ = "./";
  save_file_stub_ = "file";
  number_of_iterations = 0;
  luminosity_normalization_ = 0;
}

HourglassSimulation::~HourglassSimulation() {
 // delete stuff properly here
}

int HourglassSimulation::InitSpacetime() {
  // Defining spatial corrdinates, to be used in arrays
  z_low = -300.0; // physical length of three beam buckets (one bunch)
  z_high = 300.0;
  z_range = z_high - z_low;
  x_low = -0.3;
  x_high = 0.3;
  x_range = x_high - x_low;
  y_low = -0.3;
  y_high = 0.3;
  y_range = y_high - y_low;
  t_low =  -10.e-8; // Bunches take 106 nanoseconds to cross-eachother 
  t_high =  10.e-8; // This value is in nanoseconds already
  t_range = t_high - t_low;
  vel = 2.99792e+10; // speed of light (cm/s)
  N_bin_t = 81; 
  N_bin_z = 599; // ZDC timing resolution is 100 picoseconds, so we set z-vertex resoltuion to 3 cm.
  N_bin_x = 59;
  N_bin_y = 59;
  binsizeZ = z_range/N_bin_z;
  binsizeX = x_range/N_bin_x;
  binsizeY = y_range/N_bin_y;
  binsizeT = t_range/N_bin_t;

  z_position_.resize(N_bin_z,0.0); 
  x_position_.resize(N_bin_x,0.0); 
  y_position_.resize(N_bin_y,0.0); 
  t_position_.resize(N_bin_t,0.0);

  for(int t_ct = 0; t_ct < N_bin_t; t_ct++) { 
    t_position_[t_ct] = (t_low + static_cast<double>(t_ct)*binsizeT + binsizeT/2.0);
  }
  for(int z_ct = 0; z_ct < N_bin_z; z_ct++) { 
    z_position_[z_ct] = (z_low + static_cast<double>(z_ct)*binsizeZ + binsizeZ/2.0); 
  }
  for(int x_ct = 0; x_ct < N_bin_x; x_ct++) {
    x_position_[x_ct] = (x_low + static_cast<double>(x_ct)*binsizeX + binsizeX/2.0); 
  }
  for(int y_ct = 0; y_ct < N_bin_y; y_ct++) {
    y_position_[y_ct] = (y_low + static_cast<double>(y_ct)*binsizeY + binsizeY/2.0); 
  }

  std::cout << "Spacetime Discreetization: " << std::endl 
    << "X Resolution: " << binsizeX << std::endl
    << "Y Resolution: " << binsizeY << std::endl
    << "Z Resolution: " << binsizeZ << std::endl
    << "T Resolution: " << binsizeT << std::endl
    << "T Res (space): " << binsizeT*vel << std::endl;
  return 0;
}

int HourglassSimulation::ShowConfig() {
  std::cout << "HOURGLASS SIMULATION CONFIGURATION" << std::endl;
  std::cout << std::left << std::setw(30) << "RUN_NUMBER" << std::left << std::setw(30) << run_number_ << std::endl;
  for(auto i = config_.par_.begin(); i != config_.par_.end(); ++i) {
    auto name = i->first;
    auto val  = i->second;
    std::cout << std::left << std::setw(30) <<  name << std::left << std::setw(30) << val << std::endl;
  }
  return 0;
}

int HourglassSimulation::InitConfig() {
  run_number_ = config_.GetPar("RUN_NUMBER");
  count_norm  = std::stoi(config_.GetPar("ZDC_COUNTS"));
  xoff        = std::stod(config_.GetPar("X_OFFSET"));
  yoff        = std::stod(config_.GetPar("Y_OFFSET"));
  sigma_x     = std::stod(config_.GetPar("HORIZONTAL_BEAM_WIDTH"));
  sigma_y     = std::stod(config_.GetPar("VERTICAL_BEAM_WIDTH")  );
  N_blue      = std::stod(config_.GetPar("AVG_NUMBER_IONS_BLUE_BEAM"));
  N_yell      = std::stod(config_.GetPar("AVG_NUMBER_IONS_YELLOW_BEAM"));
  z_vtx_off   = std::stod(config_.GetPar("BBC_ZDC_Z_VERTEX_OFFSET"));
  beta_star   = std::stod(config_.GetPar("BETA_STAR" ));
  angle_xz    = std::stod(config_.GetPar("CROSSING_ANGLE_XZ"));
  angle_yz    = std::stod(config_.GetPar("CROSSING_ANGLE_YZ"));
  n_bunch     = std::stod(config_.GetPar("FILLED_BUNCHES"));
  freq        = std::stod(config_.GetPar("BUNCH_CROSSING_FREQUENCY"));
  scale       = std::stod(config_.GetPar("Z_PROFILE_SCALE_VALUE"));
  multi_coll  = std::stod(config_.GetPar("MULTIPLE_COLLISION_RATE"));
  max_coll    = std::stod(config_.GetPar("MAX_COLLISIONS"));
  MAX_COLL    = max_coll + 1;
  sigma_zl    = std::stod(config_.GetPar("Z_BUNCH_WIDTH_LEFT_GAUSSIAN")  );
  sigma_zr    = std::stod(config_.GetPar("Z_BUNCH_WIDTH_RIGHT_GAUSIAN")  );
  sigma_zc    = std::stod(config_.GetPar("Z_BUNCH_WIDTH_CENTRAL_GAUSIAN"));
  mu_zl       = std::stod(config_.GetPar("Z_BUNCH_WIDTH_LEFT_OFFSET")    );
  mu_zr       = std::stod(config_.GetPar("Z_BUNCH_WIDTH_RIGHT_OFFSET")   );
  sc_sigma_zl = sigma_zl*scale; // applying constants in lumi calculation
  sc_sigma_zr = sigma_zr*scale; // applying constants in lumi calculation
  sc_sigma_zc = sigma_zc*scale; // applying constants in lumi calculation
  sc_mu_zl    = mu_zl   *scale; // applying constants in lumi calculation
  sc_mu_zr    = mu_zr   *scale; // applying constants in lumi calculation
  zdc_compare_histo_name_ = config_.GetPar("ZDC_VERTEX_DISTRIBUTION_NAME");

  // Swap x and y variables for vertical scan
  // Config file will still have the correct values, but the simulation only
  // does the scan in one dimension.
  if(yoff > 0){
    double temp_sig_x = sigma_x;
    double temp_sig_y = sigma_y;
    double temp_x_off = xoff;
    double temp_y_off = yoff;
    double angle_xz_temp = angle_xz;
    double angle_yz_temp = angle_yz;

    sigma_x = temp_sig_y;
    sigma_y = temp_sig_x;
    xoff = temp_y_off;
    yoff = temp_x_off;
    angle_xz = angle_yz_temp;
    angle_yz = angle_xz_temp;
  }
  return 0;
}

int HourglassSimulation::InitFromConfig(const std::string& config_file) {
  HourglassConfiguration config;
  // in case for whatever god-forsaken reason you use a config file without all
  // the config parameters...
  config.SetDefaultValues(); 
  config.LoadConfigFile(config_file);
  config_ = config;
  InitConfig();
  return 0;
}

int HourglassSimulation::InitDefaultConfig() {
  HourglassConfiguration config;
  config.SetDefaultValues();
  config_ = config;
  InitConfig();
  return 0;
}

int HourglassSimulation::Compare() {
  auto pars = config_.GetAllPar();
  for(auto i = pars.begin(); i != pars.end(); ++i) {
    std::string par = i->first;
    std::string val = i->second;
    std::string line = par + " " + val;
    config_text->AddText(line.c_str());
  }
  config_and_vertex_compare->cd(1);
  config_text->Draw();
  config_and_vertex_compare->cd(2);
  zdc_zvertex_dat->DrawCopy();
  zdc_zvertex_sim->DrawCopy("same");

  simulation_config_canvas->cd();
  config_text->Draw();

  zvertex_comparison_canvas->cd();
  zdc_zvertex_dat->DrawCopy();
  zdc_zvertex_sim->DrawCopy("same");

  zvertex_comparison_canvas->Update();
  config_and_vertex_compare->Update();
  simulation_config_canvas->Update();
  // Get Chi2Test
  // WW - both histograms are weighted
  // P  - Prints: chi2, ndf, p_value, igood
  // CHI2/NDF - Returns Chi2/NDF instead of p_value
  chi2_test = zdc_zvertex_dat->Chi2Test(zdc_zvertex_sim,"WW P CHI2/NDF");

  // Generate average residual
  double average_res = 0;
  for(int i = 0; i < zdc_zvertex_dat->GetNbinsX(); i++) {
    average_res +=  pow(zdc_zvertex_dat->GetBinContent(i) - zdc_zvertex_sim->GetBinContent(i),2.0);
  }
  squares_residual = average_res/(zdc_zvertex_dat->GetNbinsX());
  return 0;
}

// Hardcode Configuration here, used for last-resort debugging
int HourglassSimulation::ManualInit() {
  // DEFAULT Configuration - max x offset for 359711
  run_number_ = "359711";
  xoff = -0.1; 
  yoff = 0.0;  
  count_norm = 592;
  multi_coll = 0.001; 
  MAX_COLL = 5 + 1;
  n_bunch = 107;
  freq = 78213.0;
  scale = 1.5; // What is this for..? wp
  N_blue = 120.029e9; // blue
  N_yell = 88.167e9;  // yellow
  sigma_xstar = 0.0245674/sqrt(2.0); // From run 12 scan (no weighting on beam width)
  sigma_ystar = 0.0238342/sqrt(2.0); 
  sc_sigma_zl = 35.15*scale; // need to use WCM profile directly.
  sc_sigma_zc = 27.65*scale; 
  sc_sigma_zr = 55.95*scale;  
  sc_mu_zl = -70.2*scale ; 
  sc_mu_zr =  56.7*scale; 
  z_vtx_off = 9.38; // 
  beta_star = 85; // 85
  angle_xz = -0.08e-3; // wp was -0.08e-3
  return 0;
}

int HourglassSimulation::Factorial(int n){ 
  if (n > 1) return n * Factorial(n - 1);
  else return 1;
}

double HourglassSimulation::SmearZVertex(double rand_prob_res, double orig_z){
  double temp_z; 
  const double sigma_res = 15.0;
  double res_prob;
  double sum_res_prob[101];
  double add_prob = 0.0;
  double smeared_z = -999.9;

  //making gaussian for smearing resolution
  for(int i=0; i<=100; i++){
    temp_z = -50.0 + 1.0*i;
    res_prob = exp(-0.5*(pow((temp_z)/(sigma_res), 2.0)));
    add_prob += res_prob;
    sum_res_prob[i] = add_prob;
  }

  //getting smeared position from gaussian
  rand_prob_res = rand_prob_res*sum_res_prob[100];//normalizing rand no 0<=r<=1 to 0<=r<=tot_prob
  bool condres = false;
  int l = 0;
  while((condres == false) && (l<=99)){
    if(rand_prob_res <= sum_res_prob[0]){
      smeared_z = -50.0;
      condres = true;
    } else if((rand_prob_res>sum_res_prob[l]) && (rand_prob_res<=sum_res_prob[l + 1])){
      smeared_z = -50.0 + 1.0*(l+1);
      condres = true;
    } else  {
      condres = false;
    }
    l++;
  }
  smeared_z = smeared_z + orig_z;
  return smeared_z;
}

int HourglassSimulation::RunRootFinder(int model_opt ,const std::string& compare_file) {
  // these are for 200GeV Run15 pp running
  // COLLISIONS/BUNCH CROSSING:  0.435, 0.402, 0.267, 0.126, 0.027, 0.001
  // BEAM OFFSET (mm)         :  0.0  , 0.01 , 0.025, 0.040, 0.060, 0.090
  
  TGraph* g_mc_guess = new TGraph();
  g_mc_guess->SetName("g_mc_guess" );
  g_mc_guess->SetTitle("Multiple Collisions Guess;Beam Offset;Multiple Collision Rate Per Bunch Crossing");
  g_mc_guess->SetPoint(0, 0.0  , 0.435);
  g_mc_guess->SetPoint(1, 0.01 , 0.402);
  g_mc_guess->SetPoint(2, 0.025, 0.267);
  g_mc_guess->SetPoint(3, 0.040, 0.126);
  g_mc_guess->SetPoint(4, 0.060, 0.027);
  g_mc_guess->SetPoint(5, 0.090, 0.001);
  g_mc_guess->SetPoint(6, 0.1  , 0.001);

  bool not_converged = true;
  double prev_residual = squares_residual;
  double prev_chi2     = chi2_test;
  double mc_center = 0.;
  
  //Define ranges - these must contain the "true values"
  if(fabs (xoff) > fabs(yoff)) {
    mc_center = g_mc_guess->Eval(fabs(xoff));
  } else {
    mc_center = g_mc_guess->Eval(fabs(yoff));
  }
  

  std::cout << "Note, overriding intitial Multiple Collision Rate Guess before running: " << mc_center << std::endl;
  multi_coll = mc_center;
  
  Run(model_opt); // Run once to get initial values.
  Compare();
  SaveFigures();

  double mc_rate_min = mc_center - 0.5*mc_center;
  double mc_rate_max = mc_center + 0.5*mc_center;
  double mc_rate_mid = (mc_rate_min+mc_rate_max)/2.0;
  double mc_rate_step = (mc_rate_max - mc_rate_mid)*0.5;
  double beta_max = beta_star+(beta_star*0.1);
  double beta_mid = beta_star;
  double beta_step = (beta_max - beta_mid)*0.5;
  double angle_max = 0.0025;
  double angle_mid = 0.0;
  double angle_step = (angle_max - angle_mid)*0.5;
  
  std::vector<double> v_mc;
  std::vector<double> v_beta;
  std::vector<double> v_angle;

  v_mc   .push_back(mc_rate_mid - mc_rate_step);
  v_mc   .push_back(mc_rate_mid + mc_rate_step);
  v_beta .push_back(beta_mid    - beta_step );
  v_beta .push_back(beta_mid                );
  v_beta .push_back(beta_mid    + beta_step );
  v_angle.push_back(angle_mid   - angle_step);
  v_angle.push_back(angle_mid               );
  v_angle.push_back(angle_mid   + angle_step);

  double best_beta = 0.;
  double best_ang = 0.;
  double best_residual = 0.;
  double best_chi2 = 0.;
  double best_mc = 0.;

  double prev_beta  = beta_star  ;
  double prev_ang   = angle_xz   ;
  double prev_mc    = multi_coll;

  // double curr_sig_x = 0.;
  double curr_beta  = 0.;
  double curr_ang   = 0.;
  double curr_mc    = 0.;

  std::map<double,HourglassConfiguration> least_squares_map;
  std::map<double,HourglassConfiguration> chi2_map;

  least_squares_map[squares_residual] = config_;
  chi2_map[chi2_test] = config_;

  number_of_iterations++;
  while(not_converged) {
    bool first_test = true;
    for(auto b = v_beta.begin(); b != v_beta.end(); ++b) {
      for(auto a = v_angle.begin(); a != v_angle.end(); ++a) {
        for(auto m = v_mc.begin(); m != v_mc.end(); ++m) {
          // Here, we set the variables directly
          config_.ModifyConfigParameter("BETA_STAR",std::to_string(*b));
          if(fabs(yoff) > 0) {
            config_.ModifyConfigParameter("CROSSING_ANGLE_YZ",std::to_string(*a));
          } else {
            config_.ModifyConfigParameter("CROSSING_ANGLE_XZ",std::to_string(*a));
          }
          config_.ModifyConfigParameter("MULTIPLE_COLLISION_RATE",std::to_string(*m));
          InitConfig(); // Ensure that the config file, and the internal varialbes are set to the same thing.

          // Run a new model after setting the member variables
          ResetModel();
          Run(model_opt);
          Compare();
          if(save_all) {
            SaveFigures();
          }
          if(first_test) { 
            best_beta = beta_star;
            best_ang = angle_xz;
            best_residual = squares_residual;
            best_chi2 = chi2_test;
            best_mc = multi_coll;
            first_test = false;
          } else {
            if(squares_residual < best_residual) {
              best_beta     = beta_star;
              best_ang      = angle_xz;
              best_residual = squares_residual;
              best_chi2     = chi2_test;
              best_mc       = multi_coll;
            }
          }
        }// mc loop 
      }//xing angle loop
    }//beta star loop
    beta_step = 0.5*beta_step;
    angle_step = 0.5*angle_step;
    mc_rate_step = 0.5*mc_rate_step;
    v_beta .clear();
    v_angle.clear();
    v_mc   .clear();

    // Redefine the search area
    v_mc   .push_back(best_mc - mc_rate_step);
    v_mc   .push_back(best_mc + mc_rate_step);
    v_beta .push_back(best_beta  - beta_step);
    v_beta .push_back(best_beta);
    v_beta .push_back(best_beta  + beta_step);
    v_angle.push_back(best_ang   - angle_step);
    v_angle.push_back(best_ang);
    v_angle.push_back(best_ang   + angle_step);

    std::string save_file_stub_orig = save_file_stub_;
    std::stringstream ss;
    ss <<  save_file_stub_ << "_iteration_" <<  std::setfill('0') << std::setw(3) << std::to_string(number_of_iterations);
    save_file_stub_ = ss.str();
    SaveFigures();
    save_file_stub_ = save_file_stub_orig;
    ss.str("");

    // Convergence Testing 
    curr_beta  = beta_star  ;
    curr_ang   = angle_xz   ;
    curr_mc    = multi_coll ;

    double beta_diff  = fabs(curr_beta  - prev_beta )/curr_beta ;
    double ang_diff   = fabs(curr_ang   - prev_ang  )/curr_ang  ;
    double mc_diff    = fabs(curr_mc - prev_mc      )/curr_mc   ;

    prev_beta  = curr_beta  ;
    prev_ang   = curr_ang   ;
    prev_mc    = curr_mc    ;

    std::cout << "beta_diff : " << beta_diff     << std::endl;
    std::cout << "ang_diff  : " << ang_diff      << std::endl;
    std::cout << "mc_diff   : " << mc_diff       << std::endl;
    std::cout << "Chi2      : " << chi2_test     << std::endl;
    std::cout << "residual  : " << best_residual << std::endl;

    double convergence_res  = fabs(best_residual - prev_residual );
    double convergence_chi2 = fabs(best_chi2     - prev_chi2     );
    std::cout << "Iteration: " << number_of_iterations 
      << ", Chi2 Convergence: " << convergence_chi2 
      << ", Residual Convergence: " << convergence_res << std::endl;
    std::cout << "Convergence Variables: chi2_test: " << chi2_test << ", squares_residual: " << squares_residual << std::endl;
    std::cout << "Search Variables (beta_star, angle_xz, multi_coll): " << beta_star << ", " << angle_xz << ", " << multi_coll << std::endl;
    if(number_of_iterations == 10 ) { // this value is determined experimentally by observing behavior of chi2 vs iteraiton.
      not_converged = false;
    } else { 
      prev_residual = best_residual;
      prev_chi2 = best_chi2;
    }
    number_of_iterations++;
  }// while loop
  std::cout << "Best Parameters: " << std::endl
    << "best_beta     : " << best_beta     << std::endl
    << "best_ang      : " << best_ang      << std::endl
    << "best_residual : " << best_residual << std::endl
    << "best_chi2     : " << best_chi2     << std::endl
    << "best_mc       : " << best_mc       << std::endl;
  return 0;
}


int HourglassSimulation::Run(int model_opt = 0) {
  time_tracker[0] = GetTime().count();
  std::cout << "phase " << 0 << std::endl;
  ShowConfig();
  GenerateNewModelFromFit();
  time_tracker[1] = GetTime().count();
  std::cout << "phase " << 1 << std::endl;
  GenerateZVertexProfile();
  time_tracker[2] = GetTime().count();
  // FINISHED - show some timing statistics
  auto first_itr = time_tracker.begin();
  auto second_itr = time_tracker.begin();
  ++second_itr;
  std::cout << "TIME Analysis: " << std::endl;
  while(first_itr != time_tracker.end() && second_itr != time_tracker.end()) {
    int phase = first_itr->first; 
    double time = second_itr->second - first_itr->second;
    std::cout << "Phase " << phase  << " time: " << time << " ms (" << (double)time/1000. << " s)" << std::endl;
    ++first_itr;
    ++second_itr;
  }
  std::cout << "In the modeling loop, we performed: " << how_many_things << " iterations." << std::endl;
  how_many_runs++;
  return 0;
}

int HourglassSimulation::LoadZProfile(const std::string& blue_f_name, const::std::string& yell_f_name, const std::string& fit_file_name) {
  std::ifstream z_bunch_profile_blue(blue_f_name.c_str());
  std::ifstream z_bunch_profile_yell(yell_f_name.c_str());
  double z_pos = 0;
  double density = 0;

  while(z_bunch_profile_blue >> z_pos >> density) {
    z_profile_blue_[z_pos] = density;
    z_bunch_profile_blue_->SetPoint(z_bunch_profile_blue_->GetN(), z_pos, density);
  }
  z_pos = 0;
  density = 0;
  while(z_bunch_profile_yell >> z_pos >> density) {
    z_profile_yell_[z_pos] = density;
    z_bunch_profile_yell_->SetPoint(z_bunch_profile_yell_->GetN(), z_pos, density);
  }

  TFile* f = new TFile(fit_file_name.c_str(),"READ");
  if(!f) {
    std::cout << "was not able to open the file containing wcm fits! You will not be able to run ::GenerateNewModelFromFit!" << std::endl;
  }
  f_z_profile_blue_ = (TF1*)f->Get("blue_zprofile");
  f_z_profile_yell_ = (TF1*)f->Get("yellow_zprofile");
  f_simple_gaus_zprofile_blue_ = (TF1*)f->Get("simple_gaus_blue");
  f_simple_gaus_zprofile_yell_ = (TF1*)f->Get("simple_gaus_yell");
  f_simple_gaus_zprofile_blue_->SetRange(-6000,6000);
  f_simple_gaus_zprofile_yell_->SetRange(-6000,6000);
  f_z_profile_blue_->SetRange(-6000,6000); // to account for travelling around the PHENIX IR
  f_z_profile_yell_->SetRange(-6000,6000);
  save_registry_.push_back(f_z_profile_blue_);
  save_registry_.push_back(f_z_profile_yell_);
  save_registry_.push_back(z_bunch_profile_blue_);
  save_registry_.push_back(z_bunch_profile_yell_);
  return 0;
}

double HourglassSimulation::BetaSqueeze(double beam_width, double z ){
  // Apply Beta Squeeze
  double sigma = beam_width*sqrt(1+pow(z/beta_star,2.0));
  return sigma;
}

void HourglassSimulation::RotateBlueCoordinates(double& c, double& z, double angle) {
  double c_orig = c;
  double z_orig = z;
  c = c_orig*cos(angle/2.)-z_orig*sin(angle/2.);
  z = z_orig*cos(angle/2.)+c_orig*sin(angle/2.);
}

void HourglassSimulation::RotateYellCoordinates(double& c, double& z, double angle) {
  double c_orig = c;
  double z_orig = z;
  c = c_orig*cos(angle/2.)+z_orig*sin(angle/2.);
  z = z_orig*cos(angle/2.)-c_orig*sin(angle/2.);
}

void HourglassSimulation::NormalizeDensity() {
  density_normalization_ = 0.;
  double z_temp = 0.;
  double volume = binsizeX*binsizeY*binsizeZ;
  for(auto x_i = x_position_.begin(); x_i != x_position_.end(); ++x_i){
    for(auto z_i = z_position_.begin(); z_i != z_position_.end(); ++z_i){
      for(auto y_i = y_position_.begin(); y_i != y_position_.end(); ++y_i){
        density_normalization_ += 
          GetGaussianDensity(*x_i, sigma_x)
          *GetGaussianDensity(*y_i,sigma_y)
          *LookupZDensityBlue(*z_i,z_temp)
          *volume;
      }
    }
  }
  // Now get the luminosity normalization, while we're at it.
  luminosity_normalization_ = 
    2.0*
    (0.5/density_normalization_)
    *N_blue*N_yell*freq*n_bunch;
}

double HourglassSimulation::GetGaussianDensity(double x, double sigma){
  return exp(-0.5*pow(x/sigma,2.0));
}

// Assumes that we have our wcm density binned at a resolution of
// 1.5 cm
double HourglassSimulation::LookupZDensityBlue(double z, double& found_z) {
  auto first_element = z_profile_blue_.begin();
  auto last_element  = z_profile_blue_.end();
  --last_element; 
  if( z < first_element->first || z > last_element->first) {
    return 0.;
  }
  auto found = z_profile_blue_.lower_bound(z);
  if(found == last_element && fabs(found->first - z) > 3.) { 
    std::cout << "Lookup failed, there is no density to use!" << std::endl
      << "Lookup Z: " << z << std::endl;
  }
  found_z = found->first; 
  return found->second;
}

double HourglassSimulation::LookupZDensityYell(double z, double& found_z) {
  auto first_element = z_profile_yell_.begin();
  auto last_element  = z_profile_yell_.end();
  --last_element; 
  if( z < first_element->first || z > last_element->first) {
    return 0.;
  }
  auto found = z_profile_yell_.lower_bound(z);
  if(found == last_element && fabs(found->first - z) > 3.) { 
    std::cout << "Lookup failed, there is no density to use!" << std::endl
      << "Lookup Z: " << z << std::endl;
  }
  found_z = found->first; 
  return found->second;
}

int HourglassSimulation::CreateCumulativePoissonDistribution() {
  double tot_prob = 0.0; 
  //================ creating an array with cumulative Poisson Disribution ======================
  for(int no_count = 0; no_count < MAX_COLL; ++no_count) { 
    double p = ((exp(-multi_coll))*(pow(multi_coll, static_cast<double>(no_count))))/Factorial(no_count);
    if(p > 0.0) {
      tot_prob += p;
    } else {
      std::cout << "Poisson prob. is negative. Check code." <<  std::endl;
      return 0;
    }
    if(tot_prob > 1.0){
      std::cout << "check poisson prob calculation. it is not normalised." << std::endl;
      std::cout << "prob is: " << tot_prob << std::endl;
    } else {
      poisson_dist_[no_count] = tot_prob;
    }
    //std::cout << poisson_dist_[no_count] << std::endl;//debug
  }
  std::cout << "done accumulating Poisson distbn." <<std::endl;
  return 0;
}

// Amaresh's last known "working" model. Critically, note how normalization was
// done.
void HourglassSimulation::GenerateAmareshModel() {
  std::cout << "Using Amaresh Beam Profile Model" << std::endl;
  double sum_T;
  double add_T;
  double sigma_xstar = sigma_x/sqrt(2.);
  double sigma_ystar = sigma_y/sqrt(2.);
  double half_angle_xz = angle_xz/2.0;
  double cos_half_angle_xz = cos(angle_xz/2.0);
  double spacetime_volume = binsizeX*binsizeY*binsizeZ*vel*binsizeT;

  sum_T = 0.0;
  luminosity_tot_ = 0.0;

  // Calculate Luminosity
  for(int ct=0; ct<N_bin_t; ct++) {
    double t = t_position_[ct];
    add_T = 0.0;
    z_norm_[ct] = 0.0;
    for(int cz=0; cz<N_bin_z; cz++) {
      double z = z_position_[cz];
      //corrected for angle_xz dependence, 2015
      double sigma_xz = sigma_xstar*sqrt(1 + pow(cos_half_angle_xz*z/beta_star, 2.0));

      //corrected for angle_xz dependence, 2015
      double sigma_yz = sigma_ystar*sqrt(1 + pow(cos_half_angle_xz*z/beta_star, 2.0));

      // Z-profile density for first bunch
      double density_z_blue_left = 
        1.522*exp(-0.5*(pow((z*cos_half_angle_xz-sc_mu_zl-vel*t)/sc_sigma_zl,2)))
        *(1.0/sc_sigma_zl);//corrected for rotaion in x-z plane, 2015
      double density_z_blue_center = 
        2.157*exp(-0.5*(pow((z*cos_half_angle_xz-vel*t)/sc_sigma_zc,2)))
        *(1.0/sc_sigma_zc);//corrected for rotaion in x-z plane, 2015
      double density_z_blue_right =
        1.999*exp(-0.5*(pow((z*cos_half_angle_xz-sc_mu_zr-vel*t)/sc_sigma_zr,2)))
        *(1.0/sc_sigma_zr);

      // Z-profile density for second bunch
      double density_z_yell_left =
        1.522*exp(-0.5*(pow((z*cos_half_angle_xz+sc_mu_zl+vel*t)/sc_sigma_zl,2)))
        *(1.0/sc_sigma_zl);
      double density_z_yell_center = 
        2.157*exp(-0.5*(pow((z*cos_half_angle_xz+vel*t)/sc_sigma_zc,2)))
        *(1.0/sc_sigma_zc);
      double density_z_yell_right = 
        1.999*exp(-0.5*(pow((z*cos_half_angle_xz+sc_mu_zr+vel*t)/sc_sigma_zr,2)))
        *(1.0/sc_sigma_zr); 

      double density_blue_z = density_z_blue_left + density_z_blue_center + density_z_blue_right;
      double density_yell_z = density_z_yell_left + density_z_yell_center + density_z_yell_right;

      if( ! amaresh_model_run ){
        if(ct == 40){ // set for t == 0
          amaresh_z_blue->SetPoint( amaresh_z_blue->GetN(), z, density_blue_z);
          amaresh_z_yell->SetPoint( amaresh_z_yell->GetN(), z, density_yell_z);
        }

        if(!new_model_run && !amaresh_model_run ) {
          if(ct == 35){
            new_model_z_blue[0]->SetPoint( new_model_z_blue[0]->GetN(), z, density_blue_z);
            new_model_z_yell[0]->SetPoint( new_model_z_yell[0]->GetN(), z, density_yell_z);
          }
          if(ct == 40){
            new_model_z_blue[1]->SetPoint( new_model_z_blue[1]->GetN(), z, density_blue_z);
            new_model_z_yell[1]->SetPoint( new_model_z_yell[1]->GetN(), z, density_yell_z);
          }
          if(ct == 45){
            new_model_z_blue[2]->SetPoint( new_model_z_blue[2]->GetN(), z, density_blue_z);
            new_model_z_yell[2]->SetPoint( new_model_z_yell[2]->GetN(), z, density_yell_z);
          }
        }
      }

      luminosity_normalization_ = ((n_bunch*freq*N_blue*N_yell)/pow(2*PI, 1.5))/pow(sigma_xz*sigma_yz, 2.0);
      for(int cx=0; cx<N_bin_x; cx++) {
        double x = x_position_[cx];
        double density_x1 = exp(-0.5*pow((x*cos_half_angle_xz-xoff+half_angle_xz*z)/sigma_xz, 2.0)); // only one bunch is offset
        double density_x2 = exp(-0.5*pow((x*cos_half_angle_xz     -half_angle_xz*z)/sigma_xz, 2.0));
        for(int cy=0; cy<N_bin_y; cy++) { 
          double y = y_position_[cy];
          double density_y1 = exp(-0.5*pow((y-yoff)/sigma_yz,2.0)); // offset bunch
          double density_y2 = exp(-0.5*pow((y )/sigma_yz,2.0)); // fixed bunch

          double d_luminosity_ = (luminosity_normalization_ * spacetime_volume
              *(density_y1 * density_x1 * density_blue_z) // bunch 1
              *(density_y2 * density_x2 * density_yell_z) // bunch 2
              ); // this is the summed value
          if(d_luminosity_ >= 0.0) {
            luminosity_tot_ += d_luminosity_; // this is literally just the value of the integral
            add_T += d_luminosity_; // this is the same thing as luminosity_tot_?
          } else {
            std::cout << "Luminosity integral generated negative value, there is an error in the code." << std::endl;
            break;
          }
          how_many_things++;
        }//end loop on y
      }//end loop on x
      gaussian_dist_[ct][cz] = luminosity_tot_; // storing prob for 2-D z-t grid summed over x,y
    }//end loop on z 
    z_norm_[ct] = add_T;
    sum_T += add_T;
    t_dist_[ct] = sum_T;//storing prob in t with z-prob summed over
  }//end loop on t
  // should match our actual luminosity when parameters are configured
  // correctly.
  std::cout << "Luminosity = " << luminosity_tot_ << std::endl;
  std::cout << "done accumulating Gaussian distbns." << std::endl;

  if( !amaresh_model_run) {
    amaresh_compare = new TCanvas("amaresh_compare","Comparing Real WCM Dist to Model",800,1200);
    amaresh_compare->Divide(1,2);
    amaresh_compare->cd(1);
    amaresh_z_blue->SetLineWidth(2);
    amaresh_z_yell->SetLineWidth(2);
    amaresh_z_blue->SetLineColor(kBlue+2);
    amaresh_z_yell->SetLineColor(kOrange+2);
    amaresh_z_blue->Draw("AL");
    amaresh_z_yell->Draw("L");
    amaresh_compare->cd(2);
    z_bunch_profile_blue_->Draw("AL");
    z_bunch_profile_yell_->Draw("L");
    z_bunch_profile_blue_->GetXaxis()->SetRangeUser(-300.,300.);
    amaresh_compare->Update();
    save_registry_.push_back(amaresh_compare);
  }
  amaresh_model_run = true;
}

// This is broken - because the normalization of the gaussian must be
// transformed according to the beta-squeezing.
void HourglassSimulation::GenerateSimpleGausModel() {
  std::cout << "Using Simple Gaussian Beam Profile Model" << std::endl;
  double sum_T;
  double add_T;
  double sigma_xstar = sigma_x/sqrt(2.);
  double sigma_ystar = sigma_y/sqrt(2.);
  double half_angle_xz = angle_xz/2.0;
  double cos_half_angle_xz = cos(angle_xz/2.0);
  double spacetime_volume = binsizeX*binsizeY*binsizeZ*vel*binsizeT;

  sum_T = 0.0;
  luminosity_tot_ = 0.0;

  // Calculate Luminosity
  for(int ct=0; ct<N_bin_t; ct++) {
    double t = t_position_[ct];
    add_T = 0.0;
    z_norm_[ct] = 0.0;
    for(int cz=0; cz<N_bin_z; cz++) {
      double z = z_position_[cz];
      //corrected for angle_xz dependence, 2015
      double sigma_xz = sigma_xstar*sqrt(1 + pow(cos_half_angle_xz*z/beta_star, 2.0));

      //corrected for angle_xz dependence, 2015
      double sigma_yz = sigma_ystar*sqrt(1 + pow(cos_half_angle_xz*z/beta_star, 2.0));

      // Z-profile density for first bunch
      double density_z_blue_center = GetGaussianDensity((z*cos_half_angle_xz-vel*t),138.);

      // Z-profile density for second bunch
      double density_z_yell_center = GetGaussianDensity((z*cos_half_angle_xz+vel*t),138.);

      double density_blue_z = density_z_blue_center;
      double density_yell_z = density_z_yell_center;

      if( ! simple_gaus_model_run ){
        if(ct == 40){ // set for t == 0
          gaus_z_blue->SetPoint( gaus_z_blue->GetN(), z, density_blue_z);
          gaus_z_yell->SetPoint( gaus_z_yell->GetN(), z, density_yell_z);
        }

        if(!new_model_run && !simple_gaus_model_run ) {
          if(ct == 35){
            new_model_z_blue[0]->SetPoint( new_model_z_blue[0]->GetN(), z, density_blue_z);
            new_model_z_yell[0]->SetPoint( new_model_z_yell[0]->GetN(), z, density_yell_z);
          }
          if(ct == 40){
            new_model_z_blue[1]->SetPoint( new_model_z_blue[1]->GetN(), z, density_blue_z);
            new_model_z_yell[1]->SetPoint( new_model_z_yell[1]->GetN(), z, density_yell_z);
          }
          if(ct == 45){
            new_model_z_blue[2]->SetPoint( new_model_z_blue[2]->GetN(), z, density_blue_z);
            new_model_z_yell[2]->SetPoint( new_model_z_yell[2]->GetN(), z, density_yell_z);
          }
        }
      }

      luminosity_normalization_ = ((n_bunch*freq*N_blue*N_yell)/pow(2*PI, 1.5))/pow(sigma_xz*sigma_yz, 2.0);
      for(int cx=0; cx<N_bin_x; cx++) {
        double x = x_position_[cx];
        double density_x1 = exp(-0.5*pow((x*cos_half_angle_xz-xoff+half_angle_xz*z)/sigma_xz, 2.0)); // only one bunch is offset
        double density_x2 = exp(-0.5*pow((x*cos_half_angle_xz     -half_angle_xz*z)/sigma_xz, 2.0));
        for(int cy=0; cy<N_bin_y; cy++) { 
          double y = y_position_[cy];
          double density_y1 = exp(-0.5*pow((y-yoff)/sigma_yz,2.0)); // offset bunch
          double density_y2 = exp(-0.5*pow((y )/sigma_yz,2.0)); // fixed bunch

          double d_luminosity_ = (luminosity_normalization_ * spacetime_volume
              *(density_y1 * density_x1 * density_blue_z) // bunch 1
              *(density_y2 * density_x2 * density_yell_z) // bunch 2
              ); // this is the summed value
          if(d_luminosity_ >= 0.0) {
            luminosity_tot_ += d_luminosity_; // this is literally just the value of the integral
            add_T += d_luminosity_; // this is the same thing as luminosity_tot_?
          } else {
            std::cout << "Luminosity integral generated negative value, there is an error in the code." << std::endl;
            break;
          }
          how_many_things++;
        }//end loop on y
      }//end loop on x
      gaussian_dist_[ct][cz] = luminosity_tot_; // storing prob for 2-D z-t grid summed over x,y
    }//end loop on z 
    z_norm_[ct] = add_T;
    sum_T += add_T;
    t_dist_[ct] = sum_T;//storing prob in t with z-prob summed over
  }//end loop on t
  // should match our actual luminosity when parameters are configured
  // correctly.
  std::cout << "Luminosity = " << luminosity_tot_ << std::endl;
  std::cout << "done accumulating Gaussian distbns." << std::endl;

  if( !simple_gaus_model_run) {
    gaus_compare = new TCanvas("gaus_compare","Comparing Real WCM Dist to Gaus Model",800,1200);
    gaus_compare->Divide(1,2);
    gaus_compare->cd(1);
    gaus_z_blue->SetLineWidth(2);
    gaus_z_yell->SetLineWidth(2);
    gaus_z_blue->SetLineColor(kBlue+2);
    gaus_z_yell->SetLineColor(kOrange+2);
    gaus_z_blue->Draw("AL");
    gaus_z_yell->Draw("L");
    gaus_compare->cd(2);
    z_bunch_profile_blue_->Draw("AL");
    z_bunch_profile_yell_->Draw("L");
    z_bunch_profile_blue_->GetXaxis()->SetRangeUser(-300.,300.);
    gaus_compare->Update();
    save_registry_.push_back(gaus_compare);
  }
  simple_gaus_model_run = true;
}

// This is broken - the normalization of the densities must be transormed
// according to the beta-star squeezing. I have not currently implemented a
// means of doing this.
void HourglassSimulation::GenerateNewModel() {
  std::cout << "Using WCM Data model" << std::endl;
  double sum_T;
  double add_T;
  double interaction_time = 0.;
  double half_angle_xz = angle_xz/2.0;
  double cos_half_angle_xz = cos(angle_xz/2.0);
  double spacetime_volume = binsizeX*binsizeY*binsizeZ*1.0e-9*vel*binsizeT;

  sum_T = 0.0;
  luminosity_tot_ = 0.0;
          
  std::cout << "LUMINOSITY_NORMALIZATION: " << luminosity_normalization_ << std::endl;

  // Calculate luminsosity 
  for(int ct=0; ct<N_bin_t; ct++) {
    double t = t_position_[ct];
    add_T = 0.0;
    z_norm_[ct] = 0.0;
    for(int cz=0; cz<N_bin_z; cz++) {
      double z = z_position_[cz];
      double lookup_z_blue = -9e9;
      double lookup_z_yell = -9e9;
      // double density_blue_z = GetGaussianDensity(z - vel*t,80.0); // super basic model for z
      // double density_yell_z = GetGaussianDensity(z + vel*t,80.0);
      double z_blue = z*cos_half_angle_xz-vel*t;
      double z_yell = z*cos_half_angle_xz+vel*t;
      double density_blue_z = LookupZDensityBlue(z_blue, lookup_z_blue);
      double density_yell_z = LookupZDensityYell(z_yell, lookup_z_yell);
      if(!new_model_run) {
        if( lookup_z_blue > -9e9 ) {
          z_lookup_diff_blue->Fill(fabs(z_blue-lookup_z_blue));
        }
        if( lookup_z_yell > -9e9) {
          z_lookup_diff_yell->Fill(fabs(z_yell-lookup_z_yell));
        }
      }

      // Watch the bunchs collide!
      // if(cz == 300) std::cout << "\r" <<  z-vel*t << ", " << z+vel*t << std::endl;
      //     if(cz == 300) {
      //       if(fabs(density_blue_z*density_yell_z) > 0.) {
      //         interaction_time += binsizeT;
      //       }

      if(!new_model_run) {
        if(ct == 35){
          new_model_z_blue[0]->SetPoint( new_model_z_blue[0]->GetN(), z, density_blue_z);
          new_model_z_yell[0]->SetPoint( new_model_z_yell[0]->GetN(), z, density_yell_z);
        }
        if(ct == 40){
          new_model_z_blue[1]->SetPoint( new_model_z_blue[1]->GetN(), z, density_blue_z);
          new_model_z_yell[1]->SetPoint( new_model_z_yell[1]->GetN(), z, density_yell_z);
        }
        if(ct == 45){
          new_model_z_blue[2]->SetPoint( new_model_z_blue[2]->GetN(), z, density_blue_z);
          new_model_z_yell[2]->SetPoint( new_model_z_yell[2]->GetN(), z, density_yell_z);
        }
      }
      for(int cx=0; cx<N_bin_x; cx++) {
        double x = x_position_[cx];
        double x_blue = x*cos_half_angle_xz - xoff + half_angle_xz*z;
        double x_yell = x*cos_half_angle_xz        - half_angle_xz*z;
        double density_blue_x = GetGaussianDensity(x_blue,BetaSqueeze(sigma_x,z*cos_half_angle_xz));
        double density_yell_x = GetGaussianDensity(x_yell,BetaSqueeze(sigma_x,z*cos_half_angle_xz));
        for(int cy=0; cy<N_bin_y; cy++) {    
          double y = y_position_[cy];
          double y_blue = y;
          double y_yell = y;

          double density_blue_y = GetGaussianDensity(y_blue,BetaSqueeze(sigma_y,z*cos_half_angle_xz));
          double density_yell_y = GetGaussianDensity(y_yell,BetaSqueeze(sigma_y,z*cos_half_angle_xz));

          if(fabs(luminosity_normalization_) < 0.1) {std::cout << "LUMINOSITY IS NOT PROPERLY NORMALIZED!!" << std::endl;}
          double d_luminosity_ = ( 
              luminosity_normalization_
              *spacetime_volume
              *(density_blue_x * density_blue_y * density_blue_z) // blue
              *(density_yell_x * density_yell_y * density_yell_z) // yellow
              ); // this is the summed value
          if(d_luminosity_ >= 0.0) {
            luminosity_tot_ += d_luminosity_; // this is literally just the value of the integral
            add_T += d_luminosity_; // this is the same thing as luminosity_tot_?
          } else {
            std::cout << "Luminosity integral generated negative value, there is an error in the code." <<  std::endl;
            break;
          }
          how_many_things++;
        }//end loop on y
      }//end loop on x
      gaussian_dist_[ct][cz] = luminosity_tot_; // storing prob for 2-D z-t grid summed over x,y
    }//end loop on z  
    z_norm_[ct] = add_T;
    sum_T += add_T;
    t_dist_[ct] = sum_T;//storing prob in t with z-prob summed over
  }//end loop on t

  // should match our actual luminosity when parameters are configured
  // correctly.
  std::cout << "Luminosity = " << luminosity_tot_ << std::endl;
  std::cout << "done accumulating Gaussian distbns." << std::endl;
  std::cout << "Interaction time: " << interaction_time << std::endl;
  new_model_run = true;
}

// This is currently the only method aside from Amaresh's method which handles
// the normalization of the densities under beta-star squeezing appropriately.
void HourglassSimulation::GenerateNewModelFromFit() {
  std::cout << "Using fitted WCM Data Model" << std::endl;
  double sum_T;
  double add_T;
  // Factor of sqrt(2) is coming from the fact that we are looking at the RMS
  // value of simga_x and sigma_y (which in reality change as a function of z)
  double sigma_xstar = sigma_x/sqrt(2.); 
  double sigma_ystar = sigma_y/sqrt(2.);
  double spacetime_volume = binsizeX*binsizeY*binsizeZ*vel*binsizeT;
  double luminosity_constant_ = spacetime_volume*(n_bunch*freq*N_blue*N_yell);

  double xing_angle = angle_xz;
  double local_sigma_xstar = sigma_xstar;
  double local_sigma_ystar = sigma_ystar;
  double beam_offset = xoff;
  double half_angle = xing_angle/2.0;
  double cos_half_angle = cos(xing_angle/2.0);
  if(fabs(yoff) > 0.) {
    xing_angle = angle_yz;
    local_sigma_xstar = sigma_ystar;
    local_sigma_ystar = sigma_xstar;
    beam_offset = yoff;
    xing_angle = angle_yz;
    half_angle = xing_angle/2.0;
    cos_half_angle = cos(xing_angle/2.0);
  }

  sum_T = 0.0;
  luminosity_tot_ = 0.0;

  // Calculate Luminosity
  for(int ct=0; ct<N_bin_t; ct++) {
    double t = t_position_[ct];
    add_T = 0.0;
    z_norm_[ct] = 0.0;
    for(int cz=0; cz<N_bin_z; cz++) {
      double z = z_position_[cz];
      //corrected for xing_angle dependence, 2015
      double sigma_xz = local_sigma_xstar*sqrt(1 + pow(cos_half_angle*z/beta_star, 2.0));

      //corrected for xing_angle dependence, 2015
      double sigma_yz = local_sigma_ystar*sqrt(1 + pow(cos_half_angle*z/beta_star, 2.0));

      double density_blue_z = f_z_profile_blue_->Eval(z*cos_half_angle-vel*t);
      double density_yell_z = f_z_profile_yell_->Eval(z*cos_half_angle+vel*t);

      if( ! new_fit_model_run ){
        if(ct == 40){ // set for t == 0
          gaus_z_blue->SetPoint( gaus_z_blue->GetN(), z, density_blue_z);
          gaus_z_yell->SetPoint( gaus_z_yell->GetN(), z, density_yell_z);
        }

        if(!new_model_run && !simple_gaus_model_run ) {
          if(ct == 35){
            new_model_z_blue[0]->SetPoint( new_model_z_blue[0]->GetN(), z, density_blue_z);
            new_model_z_yell[0]->SetPoint( new_model_z_yell[0]->GetN(), z, density_yell_z);
          }
          if(ct == 40){
            new_model_z_blue[1]->SetPoint( new_model_z_blue[1]->GetN(), z, density_blue_z);
            new_model_z_yell[1]->SetPoint( new_model_z_yell[1]->GetN(), z, density_yell_z);
          }
          if(ct == 45){
            new_model_z_blue[2]->SetPoint( new_model_z_blue[2]->GetN(), z, density_blue_z);
            new_model_z_yell[2]->SetPoint( new_model_z_yell[2]->GetN(), z, density_yell_z);
          }
        }
      }

      // Calcualte functional normalization
      double norm_x = sqrt(2.*PI)*sigma_xz;
      double norm_y = sqrt(2.*PI)*sigma_yz;
      double sigma_z1_blue = fabs(f_z_profile_blue_->GetParameter(2));
      double sigma_z2_blue = fabs(f_z_profile_blue_->GetParameter(5));
      double sigma_z3_blue = fabs(f_z_profile_blue_->GetParameter(8));
      double sigma_z1_yell = fabs(f_z_profile_yell_->GetParameter(2));
      double sigma_z2_yell = fabs(f_z_profile_yell_->GetParameter(5));
      double sigma_z3_yell = fabs(f_z_profile_yell_->GetParameter(8));

      double norm_z_blue = sqrt(2.*PI)*(sigma_z1_blue+sigma_z2_blue+sigma_z3_blue);
      double norm_z_yell = sqrt(2.*PI)*(sigma_z1_yell+sigma_z2_yell+sigma_z3_yell);

      luminosity_normalization_ = luminosity_constant_/(2.*norm_x*norm_y*norm_z_blue*norm_z_yell);
      
      for(int cx=0; cx<N_bin_x; cx++) {
        double x = x_position_[cx];
        double density_x1 = exp(-0.5*pow((x*cos_half_angle-beam_offset+half_angle*z)/sigma_xz, 2.0)); // only one bunch is offset
        double density_x2 = exp(-0.5*pow((x*cos_half_angle            -half_angle*z)/sigma_xz, 2.0));
        for(int cy=0; cy<N_bin_y; cy++) { 
          double y = y_position_[cy];
          double density_y1 = exp(-0.5*pow((y-yoff)/sigma_yz,2.0)); // offset bunch
          double density_y2 = density_y1; // fixed bunch

          double d_luminosity_ = (luminosity_normalization_
              *(density_y1 * density_x1 * density_blue_z) // bunch 1
              *(density_y2 * density_x2 * density_yell_z) // bunch 2
              ); // this is the summed value
          if(d_luminosity_ >= 0.0) {
            luminosity_tot_ += d_luminosity_; // this is literally just the value of the integral
            add_T += d_luminosity_; // this is the same thing as luminosity_tot_?
          } else {
            std::cout << "Luminosity integral generated negative value, there is an error in the code." << std::endl;
            break;
          }
          how_many_things++;
        }//end loop on y
      }//end loop on x
      gaussian_dist_[ct][cz] = luminosity_tot_; // storing prob for 2-D z-t grid summed over x,y
    }//end loop on z 
    z_norm_[ct] = add_T;
    sum_T += add_T;
    t_dist_[ct] = sum_T;//storing prob in t with z-prob summed over
  }//end loop on t
  // should match our actual luminosity when parameters are configured
  // correctly.
  std::cout << "Luminosity = " << luminosity_tot_ << std::endl;
  std::cout << "done accumulating Gaussian distbns." << std::endl;

  if( ! new_fit_model_run ) {
    gaus_compare = new TCanvas("gaus_compare","Comparing Real WCM Dist to Gaus Model",800,1200);
    gaus_compare->Divide(1,2);
    gaus_compare->cd(1);
    gaus_z_blue->SetLineWidth(2);
    gaus_z_yell->SetLineWidth(2);
    gaus_z_blue->SetLineColor(kBlue+2);
    gaus_z_yell->SetLineColor(kOrange+2);
    gaus_z_blue->Draw("AL");
    gaus_z_yell->Draw("L");
    gaus_compare->cd(2);
    z_bunch_profile_blue_->Draw("AL");
    z_bunch_profile_yell_->Draw("L");
    z_bunch_profile_blue_->GetXaxis()->SetRangeUser(-300.,300.);
    gaus_compare->Update();
    save_registry_.push_back(gaus_compare);
  }
  new_fit_model_run = true;
}

int HourglassSimulation::GenerateZVertexProfile() {
  //variables for random numbers
  srand(time(NULL));
  bool condz, condt, condp; // condition t, condition z, condition poisson, what are these, wp
  double rand_prob_z, rand_prob_t, temp_prob_poisson;
  double rand_prob_smear;

  //coll_no dependent arrays
  int N_coll = 0; // default initialization
  int temp_t_index[20];
  double coll_pos[20];
  double coll_time[20]; 
  double zpos = -999; // default initialization
  double smeared_zpos;
  double ztime;
  //====================running over a no. of bunch crossings==================================
  cross_count = 0;
  event_limit_count = 0;
  int nonzero_event_count = 0;//debugging
  while((cross_count < 20000000) && (event_limit_count < (unsigned int)count_norm)) { 
    //=================== counting crossing number required to reach event limit =============
    cross_count++; 
    //====================choosing no. of events for a crossing from Poisson distbn=======================
    temp_prob_poisson = fabs(static_cast<double>(rand())/RAND_MAX); // (0,1)
    condp = false;
    int l = 0;
    while((condp == false) && (l < (MAX_COLL - 1))) {
      if(temp_prob_poisson <= poisson_dist_[0]) {
        N_coll = 0;
        condp = true;
      } else if((temp_prob_poisson > poisson_dist_[l]) && (temp_prob_poisson <= poisson_dist_[l + 1])) { // this is why we set max_coll to 5 + 1
        N_coll = l + 1;
        condp = true;
      } else {
        condp = false;
      }
      l++;
    }
    if(N_coll == 0) continue;
    if(N_coll > 20) continue;
    nonzero_event_count += N_coll;   //debug
    //std::cout << "no. of collisions = " << N_coll << std::endl;//debug
    //==========intialisations of arrays depending on no. of events in a crossing============
    for(int h = 0; h < N_coll; h++) {
      coll_pos[h] = 0.0;
      coll_time[h] = 0.0;
      temp_t_index[h] = 0;    
    }
    //========================produce as many collisions at diff positions=================
    for(int coll_no = 0; coll_no < N_coll; coll_no++) {
      //=================CHOSSING TIME T FROM GAUSSIAN DISTRIBUTION=====================
      rand_prob_t = luminosity_tot_*fabs(static_cast<double>(rand())/RAND_MAX);
      condt = false;
      int find_t = 0;
      while((condt == false) && (find_t<(N_bin_t - 1))) {
        if(rand_prob_t <= t_dist_[0]) {
          coll_time[coll_no] = t_position_[0];
          temp_t_index[coll_no] = 0;
          condt = true;
        } else if((rand_prob_t > t_dist_[find_t]) && (rand_prob_t <= t_dist_[find_t + 1])) {
          coll_time[coll_no] = t_position_[find_t + 1];
          temp_t_index[coll_no] = find_t + 1;
          condt = true;
        } else {
          condt = false;
        }
        find_t ++;
      }
      ztime = coll_time[coll_no];

      //=================CHOSSING Z POSITION FROM GAUSSIAN DISTRIBUTION=====================
      for(int cz=0; cz<N_bin_z; cz++) {
        z_dist_[cz] = gaussian_dist_[temp_t_index[coll_no]][cz];
      }//selecting the z-column with t-row fixed

      rand_prob_z = gaussian_dist_[temp_t_index[coll_no]][0] + (z_norm_[temp_t_index[coll_no]])*fabs(static_cast<double>(rand())/RAND_MAX);
      condz = false;
      int find_z = 0;
      while((condz == false) && (find_z < (N_bin_z - 1))) {
        if(rand_prob_z <= z_dist_[0]) {
          coll_pos[coll_no] = z_position_[0];
          condz = true;
        } else if((rand_prob_z > z_dist_[find_z]) && (rand_prob_z <= z_dist_[find_z + 1])) {
          coll_pos[coll_no] = z_position_[find_z + 1];
          condz = true;
        } else {
          condz = false;
        }
        find_z ++;
      }
      zpos = coll_pos[coll_no];
      if(fabs(zpos) > 300.0) {
        std::cout << "generated an event out of position range. check code. z = " << zpos << std::endl;
      } 
      if(fabs(ztime) > 10) {
        std::cout << "generated an event out of time range. check code. t = " << ztime << std::endl;
      }
    }// end loop over no. of colls in each event

    //getting smearing corrected positions after 150 cm zdcwide online cut
    rand_prob_smear = fabs(static_cast<double>(rand())/RAND_MAX);
    smeared_zpos = SmearZVertex(rand_prob_smear, zpos);
    //std::cout << "rand_smear: " << rand_prob_smear << " zpos: " << zpos << " smeared: " << smeared_zpos << std::endl;

    if(smeared_zpos>-500.0) {//checking for junk value from smearing function
      event_limit_count++;//counting final recorded vertices
      zdc_zvertex_sim->Fill(smeared_zpos+z_vtx_off);
    } else {
      std::cout << "could not get smeared posn, check code." << std::endl;
    }
  }// end loop over no. of events
  return 0;
}

int HourglassSimulation::Init(
    const std::string& cfg_file, 
    const std::string& compare_file_name, 
    const std::string& z_profile_blue, 
    const std::string& z_profile_yell,
    const std::string& fit_file_name
) {
  InitFromConfig(cfg_file);
  InitSpacetime();
  InitProbabilityVariables();
  CreateCumulativePoissonDistribution();
  
  zdc_zvertex_sim = 
    new TH1F("zdc_zvertex_sim",
        "Z Vertex ZDC Profile Simulation;z vertex;counts",
        100,-300,300);  // this must match the HourglassData zvtx
  zdc_zvertex_sim->SetDirectory(0);
  // histogram binning and range if we are 
  // to compare.
  zdc_zvertex_sim->SetLineColor(kRed);
  zdc_zvertex_sim->SetLineWidth(2);
  zdc_zvertex_sim->Sumw2();
  save_registry_.push_back(zdc_zvertex_sim);;
  z_lookup_diff_blue = new TH1F("z_lookup_diff_blue",
      "Distribution of Differences between Lookup Value and Real Value in Z-Profile (Blue)",
      100,0,1.6);
  z_lookup_diff_blue->SetDirectory(0);
  z_lookup_diff_yell = new TH1F("z_lookup_diff_yell",
      "Distribution of Differences between Lookup Value and Real Value in Z-Profile (Yellow)",
      100,0,1.6);
  z_lookup_diff_yell->SetDirectory(0);
  save_registry_.push_back(z_lookup_diff_blue);
  save_registry_.push_back(z_lookup_diff_yell);

  TFile* data = new TFile(compare_file_name.c_str(), "READ");
  if(!data) {
    std::cout << compare_file_name << " could not be opened by " << this_name
      << " for comparison." << std::endl;
    return 1;
  }
  zdc_zvertex_dat = (TH1F*)data->Get(zdc_compare_histo_name_.c_str())->Clone("zdc_zvertex_data");
  zdc_zvertex_dat->SetDirectory(0);
  zdc_zvertex_dat->SetName("zdc_zvertex_dat");
  if( !zdc_zvertex_dat ) {
    std::cout << "opened " << compare_file_name << " but couldn't find " 
      << zdc_compare_histo_name_ << " to extract." << std::endl;
    return 1;
  }
  zdc_zvertex_dat->SetLineColor(kBlue);
  zdc_zvertex_dat->SetDirectory(0);
  save_registry_.push_back(zdc_zvertex_dat);
  zvertex_comparison_canvas = 
    new TCanvas("zvertex_comparison_canvas",
        "Data/Simulation Comparison",
        800,600);
  simulation_config_canvas = 
    new TCanvas("simulation_config_canvas",
        "Configuration For Simulation",
        800,600);
  config_and_vertex_compare =
    new TCanvas("config_and_vertex_compare",
        "Configuration Parameters + ZVertex",
        1600,800);
  config_and_vertex_compare->Divide(2,1);
  save_registry_.push_back(zvertex_comparison_canvas);
  save_registry_.push_back(simulation_config_canvas);
  save_registry_.push_back(config_and_vertex_compare);
  z_bunch_profile_blue_ = new TGraph();
  z_bunch_profile_yell_ = new TGraph();
  z_bunch_profile_blue_->SetName("blue_z_bunch_profile");
  z_bunch_profile_yell_->SetName("yell_z_bunch_profile");
  z_bunch_profile_blue_->SetLineWidth(2);
  z_bunch_profile_yell_->SetLineWidth(2);
  z_bunch_profile_blue_->SetLineColor(kBlue+2);
  z_bunch_profile_yell_->SetLineColor(kOrange+2);
  amaresh_z_blue = new TGraph();
  amaresh_z_blue->SetName("amaresh_z_blue");
  amaresh_z_blue->SetTitle("amaresh_z_blue");
  amaresh_z_yell = new TGraph();
  amaresh_z_yell->SetName("amaresh_z_yell");
  amaresh_z_yell->SetTitle("amaresh_z_yell");
  save_registry_.push_back(amaresh_z_blue);
  save_registry_.push_back(amaresh_z_yell);
  
  gaus_z_blue = new TGraph();
  gaus_z_blue->SetName("gaus_z_blue");
  gaus_z_blue->SetTitle("gaus_z_blue");
  gaus_z_yell = new TGraph();
  gaus_z_yell->SetName("gaus_z_yell");
  gaus_z_yell->SetTitle("gaus_z_yell");
  save_registry_.push_back(gaus_z_blue);
  save_registry_.push_back(gaus_z_yell);
  
  for(int i = 0; i < 3; i++) {
    std::stringstream name_b;
    name_b << "new_model_z_blue_" << i;
    std::stringstream name_y;
    name_y << "new_model_z_yell_" << i;
    new_model_z_blue[i] = new TGraph();
    new_model_z_blue[i] ->SetName (name_b.str().c_str());
    new_model_z_blue[i] ->SetTitle(name_b.str().c_str());
    new_model_z_blue[i] ->SetLineWidth(2.);
    new_model_z_blue[i] ->SetLineColor(kBlue);
    new_model_z_yell[i] = new TGraph();
    new_model_z_yell[i] ->SetName (name_y.str().c_str());
    new_model_z_yell[i] ->SetTitle(name_y.str().c_str());
    new_model_z_yell[i] ->SetLineWidth(2.);
    new_model_z_yell[i] ->SetLineColor(kOrange+2);
    save_registry_.push_back(new_model_z_blue[i]);
    save_registry_.push_back(new_model_z_yell[i]);
  }
  config_text = new TPaveText(0.1,.1,.9,.9,"br");
  save_registry_.push_back(config_text);

  std::cout << "Done creating plots." << std::endl;
  LoadZProfile(z_profile_blue, z_profile_yell,fit_file_name);
  NormalizeDensity();
  CreateCumulativePoissonDistribution();
  first_init = true;
  delete data;
  return 0;
}

int HourglassSimulation::ResetModel() {
  InitProbabilityVariables();
  CreateCumulativePoissonDistribution();
  zdc_zvertex_sim->Reset();
  zdc_zvertex_sim->SetDirectory(0);
  zdc_zvertex_dat->SetDirectory(0);
  zdc_zvertex_dat->SetName("zdc_zvertex_dat");
  zdc_zvertex_sim->SetName("zdc_zvertex_sim");
  zvertex_comparison_canvas->Clear("D");
  simulation_config_canvas->Clear("D");
  config_and_vertex_compare->Clear("D");
  config_text->Clear("D");
  return 0;
}

int HourglassSimulation::InitProbabilityVariables() {
  poisson_dist_.clear();
  gaussian_dist_.clear();
  poisson_dist_.resize(MAX_COLL,0.);
  gaussian_dist_.resize(N_bin_t);
  t_dist_.resize(N_bin_t,0.);
  z_dist_.resize(N_bin_z,0.);
  z_norm_.resize(N_bin_t,0.); // maybe a bug? Is N_bin_t the same as N_bin_z ? wp
  for(auto i = gaussian_dist_.begin(); i != gaussian_dist_.end(); ++i){
    (*i).resize(N_bin_z,0.);
  }
  luminosity_tot_ = 0.;

  std::cout << "Probability density values initialized!" << std::endl;
  return 0;
}

int HourglassSimulation::SaveFigures() {
  gErrorIgnoreLevel = kWarning;
  std::stringstream root_out_name;
  std::stringstream config_out_name;
  std::stringstream pdf_out_name;


  if(save_all) {
    root_out_name   << save_directory_ << "/" << save_file_stub_ <<  "_runs_" << how_many_runs << ".root";
    pdf_out_name    << save_directory_ << "/" << save_file_stub_ <<  "_runs_" << how_many_runs << ".pdf";
    config_out_name << save_directory_ << "/" << save_file_stub_ <<  "_runs_" << how_many_runs << "_final.conf";
  } else {
    root_out_name   << save_directory_ << "/" << save_file_stub_ << ".root";
    pdf_out_name    << save_directory_ << "/" << save_file_stub_ << ".pdf";
    config_out_name << save_directory_ << "/" << save_file_stub_ << "_final.conf";
  }

  std::string out_file_name = pdf_out_name.str() + "[";
  TCanvas* booklet = new TCanvas("booklet","BBC Efficiency Plots");
  TFile* root_out = new TFile(root_out_name.str().c_str(), "RECREATE");
  booklet->Print(out_file_name.c_str());
  std::string TH1F_type = "TH1F";
  std::string TCanvas_type = "TCanvas";
  std::cout << "Saving figures to: " << save_directory_ << "/" << save_file_stub_ << "*" << std::endl;
  for(auto plot_i = save_registry_.begin(); plot_i != save_registry_.end(); ++plot_i) {
    auto draw_obj = *plot_i;
    if(!draw_obj) continue;
    root_out->cd();
    draw_obj->Write();
    if(!draw_obj) continue;
    TCanvas* c = new TCanvas();
    if(TH1F_type.compare(draw_obj->ClassName()) == 0){
      c->cd();
      ((TH1F*)draw_obj)->DrawCopy();
      c->Print(pdf_out_name.str().c_str());
    } else if (TCanvas_type.compare(draw_obj->ClassName()) == 0) {
      ((TCanvas*)draw_obj)->Print(); 
    } else {
      c->cd();
      draw_obj->Draw();
      c->Print(pdf_out_name.str().c_str());
    }
    delete c;
  }
  out_file_name = pdf_out_name.str() + "]";
  booklet->Print(out_file_name.c_str());
  root_out->Write();
  root_out->Close();
  if(root_out) delete root_out;
  delete booklet;
  gErrorIgnoreLevel = kInfo;

  std::ofstream config_out(config_out_name.str().c_str());
  for(auto i = config_.par_.begin(); i != config_.par_.end(); ++i) {
    auto conf_name = i->first;
    auto conf_val  = i->second;
    config_out << conf_name << " " << conf_val << std::endl;
  }
  return 0;
}
