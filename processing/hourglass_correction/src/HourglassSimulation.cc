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

#include "TF1.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveText.h"
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
}

HourglassSimulation::~HourglassSimulation() {
  std::cout << "Destroying: " << this_name << std::endl;
  for(auto i = save_registry_.begin(); i != save_registry_.end(); ++i){
    if(*i) delete *i;
  }
}

int HourglassSimulation::InitSpacetime() {
  // Defining spatial corrdinates, to be used in arrays
  z_low = -300.0; // this must match the data histogram
  z_high = 300.0; // this must match the data histogram
  z_range = z_high - z_low;
  x_low = -0.3;
  x_high = 0.3;
  x_range = x_high - x_low;
  y_low = -0.3;
  y_high = 0.3;
  y_range = y_high - y_low;
  t_low =  -10.0e-8;
  t_high = 10.0e-8;
  t_range = t_high - t_low;
  vel = 3.0e10; 
  N_bin_t = 80;
  N_bin_z = 600;
  N_bin_x = 60;
  N_bin_y = 60;
  binsizeZ = z_range/N_bin_z;
  binsizeX = x_range/N_bin_x;
  binsizeY = y_range/N_bin_y;
  binsizeT = t_range/N_bin_t;
  return 0;
}

int HourglassSimulation::ShowConfig() {
  std::cout << "HOURGLASS SIMULATION CONFIGURATION" << std::endl;
  std::cout
    << std::setw(30) << "RUN_NUMBER" << std::setw(20) <<  run_number_ << std::endl
    << std::setw(30) << "ZDC_COUNTS" << std::setw(20) <<  count_norm << std::endl
    << std::setw(30) << "X_OFFSET" << std::setw(20) <<  xoff << std::endl
    << std::setw(30) << "Y_OFFSET" << std::setw(20) <<  yoff << std::endl
    << std::setw(30) << "HORIZONTAL_BEAM_WIDTH" << std::setw(20) << sigma_x << std::endl
    << std::setw(30) << "VERTICAL_BEAM_WIDTH" << std::setw(20)   << sigma_y << std::endl
    << std::setw(30) << "AVG_NUMBER_IONS_BLUE_BEAM" << std::setw(20) <<  N_blue << std::endl
    << std::setw(30) << "AVG_NUMBER_IONS_YELLOW_BEAM" << std::setw(20) <<  N_yell << std::endl
    << std::setw(30) << "BBC_ZDC_Z_VERTEX_OFFSET" << std::setw(20) <<  z_vtx_off << std::endl
    << std::setw(30) << "BETA_STAR" << std::setw(20) <<  beta_star << std::endl
    << std::setw(30) << "CROSSING_ANGLE_XZ" << std::setw(20) <<  angle << std::endl
    << std::setw(30) << "FILLED_BUNCHES" << std::setw(20) <<  n_bunch << std::endl
    << std::setw(30) << "BUNCH_CROSSING_FREQUENCY" << std::setw(20) <<  freq << std::endl
    << std::setw(30) << "Z_PROFILE_SCALE_VALUE" << std::setw(20) <<  scale << std::endl
    << std::setw(30) << "MULTIPLE_COLLISION_RATE" << std::setw(20) <<  rate << std::endl
    << std::setw(30) << "MAX_COLLISIONS" << std::setw(20) <<  max_coll << std::endl
    << std::setw(30) << "Z_BUNCH_WIDTH_LEFT_GAUSSIAN" << std::setw(20) <<  sigma_zl << std::endl
    << std::setw(30) << "Z_BUNCH_WIDTH_RIGHT_GAUSIAN" << std::setw(20) <<  sigma_zr << std::endl
    << std::setw(30) << "Z_BUNCH_WIDTH_CENTRAL_GAUSIAN" << std::setw(20) <<  sigma_zc << std::endl
    << std::setw(30) << "Z_BUNCH_WIDTH_LEFT_OFFSET" << std::setw(20) <<  mu_zl  << std::endl
    << std::setw(30) << "Z_BUNCH_WIDTH_RIGHT_OFFSET" << std::setw(20) <<  mu_zr << std::endl;
  return 0;
}

int HourglassSimulation::InitConfig() {
  run_number_ = config_.GetPar("RUN_NUMBER");
  count_norm  = std::stoi(config_.GetPar("ZDC_COUNTS"));
  xoff        = std::stod(config_.GetPar("X_OFFSET"));
  yoff        = std::stod(config_.GetPar("Y_OFFSET"));
  sigma_x     = std::stod(config_.GetPar("HORIZONTAL_BEAM_WIDTH"));
  sigma_y     = std::stod(config_.GetPar("VERTICAL_BEAM_WIDTH")  );
  sigma_xstar = sigma_x/pow(2,0.5); // applying constants in lumi calculation
  sigma_ystar = sigma_x/pow(2,0.5); // applying constants in lumi calculation
  N_blue      = std::stod(config_.GetPar("AVG_NUMBER_IONS_BLUE_BEAM"));
  N_yell      = std::stod(config_.GetPar("AVG_NUMBER_IONS_YELLOW_BEAM"));
  z_vtx_off   = std::stod(config_.GetPar("BBC_ZDC_Z_VERTEX_OFFSET"));
  beta_star   = std::stod(config_.GetPar("BETA_STAR" ));
  angle       = std::stod(config_.GetPar("CROSSING_ANGLE_XZ"));
  n_bunch     = std::stod(config_.GetPar("FILLED_BUNCHES"));
  freq        = std::stod(config_.GetPar("BUNCH_CROSSING_FREQUENCY"));
  scale       = std::stod(config_.GetPar("Z_PROFILE_SCALE_VALUE"));
  rate        = std::stod(config_.GetPar("MULTIPLE_COLLISION_RATE"));
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
  return 0;
}

int HourglassSimulation::InitFromConfig(const std::string& config_file) {
  HourglassConfiguration config;
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

int HourglassSimulation::Compare(
    const std::string& compare_file_name
) {
  TFile* data = new TFile(compare_file_name.c_str(), "READ");
  if(!data) {
    std::cout << compare_file_name << " could not be opened by " << this_name
        << " for comparison." << std::endl;
    return 1;
  }
  zdc_zvertex_dat = (TH1F*)data->Get(zdc_compare_histo_name_.c_str())->Clone("zdc_zvertex_data");
  if( !zdc_zvertex_dat ) {
    std::cout << "opened " << compare_file_name << " but couldn't find " 
        << zdc_compare_histo_name_ << " to extract." << std::endl;
    return 1;
  }
  zdc_zvertex_dat->SetLineColor(kBlue);
  zdc_zvertex_dat->SetDirectory(0);
  zdc_zvertex_dat_norm = (TH1F*)zdc_zvertex_dat->Clone("zdc_zvertex_data_norm");
  zdc_zvertex_dat_norm->SetDirectory(0);
  double data_scale = 1.0/zdc_zvertex_dat_norm->Integral();
  zdc_zvertex_dat_norm->Scale(data_scale);

  zdc_zvertex_sim_norm = (TH1F*)zdc_zvertex_sim->Clone("zdc_zvertex_sim_norm");
  zdc_zvertex_sim_norm->SetDirectory(0);
  double sim_scale = 1.0/zdc_zvertex_sim_norm->Integral();
  zdc_zvertex_sim_norm->Scale(sim_scale);

  std::cout << "Normalization of data: " << data_scale << ", and sim: " 
      << sim_scale << std::endl;
  std::cout << "Norm Difference (should be small or zero): " << data_scale - sim_scale << std::endl;

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
  norm_config_and_vertex_compare = 
      new TCanvas("norm_config_and_vertex_compare",
                  "Normalized Vertex Distribution for Simulation and Data",
                  1600,800);
  norm_config_and_vertex_compare->Divide(2,1);

  TPaveText* config_text = new TPaveText(0.1,.1,.9,.9,"br");
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
  zdc_zvertex_dat->Draw();
  zdc_zvertex_sim->Draw("same");

  norm_config_and_vertex_compare->cd(1);
  config_text->Draw();
  norm_config_and_vertex_compare->cd(2);
  zdc_zvertex_dat_norm->Draw();
  zdc_zvertex_sim_norm->Draw("same");

  simulation_config_canvas->cd();
  config_text->Draw();

  zvertex_comparison_canvas->cd();
  zdc_zvertex_dat->Draw();
  zdc_zvertex_sim->Draw("same");

  // Compare distributions

  // Get Chi2Test
  chi2_test = zdc_zvertex_dat->Chi2Test(zdc_zvertex_sim,"UW P");
  
  // Generate average residual
  double average_res = 0;
  for(int i = 0; i < zdc_zvertex_dat->GetNbinsX(); i++) {
    average_res +=  pow(zdc_zvertex_dat->GetBinContent(i) - zdc_zvertex_sim->GetBinContent(i),2.0);
  }
  squares_residual =  average_res/(zdc_zvertex_dat->GetNbinsX());

  // save canvases and plots
  save_registry_.push_back(zdc_zvertex_dat);
  save_registry_.push_back(zdc_zvertex_sim_norm);
  save_registry_.push_back(zdc_zvertex_dat_norm);
  save_registry_.push_back(zvertex_comparison_canvas);
  save_registry_.push_back(simulation_config_canvas);
  save_registry_.push_back(config_and_vertex_compare);
  save_registry_.push_back(norm_config_and_vertex_compare);
  delete data;
  return 0;
}
// Hardcode Configuration here, used for last-resort debugging
int HourglassSimulation::InitDefault() {
// DEFAULT Configuration - max x offset for 359711
  run_number_ = "359711";
  xoff = -0.1; 
  yoff = 0.0;  
  count_norm = 592;
  rate = 0.001; 
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
  angle = -0.08e-3; // wp was -0.08e-3
  return 0;
}

int HourglassSimulation::InitPlots() {
  zdc_zvertex_sim = 
      new TH1F("zdc_zvertex_sim",
               "Z Vertex ZDC Profile Simulation;z vertex;counts",
               100,z_low,z_high);  // this must match the HourglassData zvtx
                                   // histogram binning and range if we are 
                                   // to compare.
  zdc_zvertex_sim->SetLineColor(kRed);
  zdc_zvertex_sim->SetLineWidth(2);
  zdc_zvertex_sim->Sumw2();
  save_registry_.push_back(zdc_zvertex_sim);
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

int HourglassSimulation::Run() {
  // Final initialization
  InitSpacetime();
  InitPlots();
  ShowConfig();

  // store various timestamps through out the execuation of the code, for 
  // bottleneck tracking.
  std::map<std::string, long long unsigned int > time_tracker;
  time_tracker["start"] = GetTime().count();
  long long unsigned int how_many_things = 0; 

  double poisson_dist[MAX_COLL];
  double g;
  double norm;
  double ex_1l, ex_1c, ex_1r, ex_2l, ex_2c, ex_2r; // dividing luminosity integral into separable gaussian pieces
  double sigma_xz, sigma_yz;
  double sum_prob;
  double sum_T;
  double add_T;

  double gaussian_dist[N_bin_t][N_bin_z];//made it static to get over the size prob.
  /** these coordinates actually represent discrete space */
  double position[N_bin_z]; 
  double xposition[N_bin_x]; 
  double yposition[N_bin_y]; 
  double input_time[N_bin_t];
  /** end 4-space coordinates */

  //variables for random numbers
  srand(time(NULL));
  bool condz, condt, condp; // condition t, condition z, condition poisoon? wp
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

  double tot_prob = 0.0; 
  //========================== creating an array with cumulative Poisson Disribution =================================
  for(int no_count = 0; no_count < MAX_COLL; ++no_count) { 
    double p = ((exp(-rate))*(pow(rate, static_cast<double>(no_count))))/Factorial(no_count);
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
      poisson_dist[no_count] = tot_prob;
    }
    //std::cout << poisson_dist[no_count] << std::endl;//debug
  }
  std::cout << "done accumulating Poisson distbn." <<std::endl;

  time_tracker["phase_1"] = GetTime().count();
  std::cout << "phase_1" << std::endl;

  // CREATING DISCREET SPACIAL ARRAYS REPRESENTING SPACETIME CORRDINATES
  for(int t_ct = 0; t_ct < N_bin_t; t_ct++) { // dimension: time
    input_time[t_ct] = (t_low + static_cast<double>(t_ct)*binsizeT + binsizeT/2.0);
  }
  for(int z_ct = 0; z_ct < N_bin_z; z_ct++) { // dimension: z
    position[z_ct] = (z_low + static_cast<double>(z_ct)*binsizeZ + binsizeZ/2.0); 
  }
  for(int x_ct = 0; x_ct < N_bin_x; x_ct++) { // dimension: x
    xposition[x_ct] = (x_low + static_cast<double>(x_ct)*binsizeX + binsizeX/2.0); 
  }
  for(int y_ct = 0; y_ct < N_bin_y; y_ct++) { // dimension: y
    yposition[y_ct] = (y_low + static_cast<double>(y_ct)*binsizeY + binsizeY/2.0); 
  }

  time_tracker["phase_2"] = GetTime().count();
  std::cout << "phase_2" << std::endl;

  //====================================separate t and z dist========================
  sum_T = 0.0;
  sum_prob = 0.0;
  double t_dist[N_bin_t];
  double z_dist[N_bin_z];
  double z_norm[N_bin_t];
  for(int ct=0; ct<N_bin_t; ct++) {
    add_T = 0.0;
    z_norm[ct] = 0.0;
    for(int cz=0; cz<N_bin_z; cz++) {
      // Why is angle seemingly dependant on position?
      double half_angle = angle/2.;
      double cos_half_angle = cos(half_angle);
      sigma_xz = sigma_xstar*sqrt(1 + pow(cos_half_angle*position[cz]/beta_star, 2.0));//corrected for angle dependence, 2015
      sigma_yz = sigma_ystar*sqrt(1 + pow(cos_half_angle*position[cz]/beta_star, 2.0));//corrected for angle dependence, 2015
      norm = ((n_bunch*freq*N_blue*N_yell)/pow(2*PI, 1.5))/pow(sigma_xz*sigma_yz, 2.0);
      double spacetime_norm = norm*binsizeX*binsizeY*binsizeZ*vel*binsizeT;
      for(int cx=0; cx<N_bin_x; cx++) {
        double ex_1l_term1 = exp(-0.5*pow((xposition[cx]*cos_half_angle-xoff+half_angle*position[cz])/sigma_xz, 2.0));
        double ex_2l_term1 = exp(-0.5*pow((xposition[cx]*cos_half_angle-half_angle*position[cz])/sigma_xz, 2.0));
        for(int cy=0; cy<N_bin_y; cy++) {    
          // ALL COMPUTATION TIME SPENT HERE 
          // 
          // What have I checked to speed up?
          // 1. I replaced the computation with a root TF1 object. This actually
          //    slowed it down. Maybe worthwhile to bring this back if we are to
          //    implement gradient-descent regression?
          // 2. Are any expressions in here getting evaluated more often than
          //    they need to?  cos_half_angle is known before this loop, so we can
          //    evaluate it outside.  We have split the ex_1l equations up to get
          //    evaluated separately only when needed.  This will be a good test to
          //    see if the compiler was optimizating this part anyway.
          //
          //   This halfed the execution time.

          double ex_term2_common = exp(-0.5*pow(yposition[cy]/sigma_yz,2.0));
          
          // Density Function for First Bunch
          ex_1l = ex_1l_term1
              *ex_term2_common
              *1.522*exp(-0.5*(pow((position[cz]*cos_half_angle-sc_mu_zl-vel*input_time[ct])/sc_sigma_zl,2)))
              *(1.0/sc_sigma_zl);//corrected for rotaion in x-z plane, 2015
          ex_1c = ex_1l_term1
              *ex_term2_common
              *2.157*exp(-0.5*(pow((position[cz]*cos_half_angle-vel*input_time[ct])/sc_sigma_zc,2)))
              *(1.0/sc_sigma_zc);//corrected for rotaion in x-z plane, 2015
          ex_1r = ex_1l_term1
              *ex_term2_common
              *1.999*exp(-0.5*(pow((position[cz]*cos_half_angle-sc_mu_zr-vel*input_time[ct])/sc_sigma_zr,2)))
              *(1.0/sc_sigma_zr);
          
          // Density Function for Second Bunch
          ex_2l = ex_2l_term1
              *ex_term2_common
              *1.522*exp(-0.5*(pow((position[cz]*cos_half_angle+sc_mu_zl+vel*input_time[ct])/sc_sigma_zl,2)))
              *(1.0/sc_sigma_zl);
          ex_2c = ex_2l_term1
              *ex_term2_common
              *2.157*exp(-0.5*(pow((position[cz]*cos_half_angle+vel*input_time[ct])/sc_sigma_zc,2)))
              *(1.0/sc_sigma_zc);
          ex_2r = ex_2l_term1
              *ex_term2_common
              *1.999*exp(-0.5*(pow((position[cz]*cos_half_angle+sc_mu_zr+vel*input_time[ct])/sc_sigma_zr,2)))
              *(1.0/sc_sigma_zr);                       
          g = (spacetime_norm*(ex_1l+ex_1c+ex_1r)*(ex_2l+ex_2c+ex_2r)); // this is the integral value
          /* END EXPENSIVE PART */

          if(g >= 0.0) {
            sum_prob += g; // this is literally just the value of the integral
            add_T += g; // this is the same thing as sum_prob?
          } else {
            std::cout << "Gaussian prob. negative. Check code." <<  std::endl;
            return 0;
          }
          how_many_things++;
        }//end loop on y
      }//end loop on x
      gaussian_dist[ct][cz] = sum_prob; // storing prob for 2-D z-t grid summed over x,y
    }//end loop on z  
    z_norm[ct] = add_T;
    sum_T += add_T;
    t_dist[ct] = sum_T;//storing prob in t with z-prob summed over
  }//end loop on t

  time_tracker["phase_3"] = GetTime().count();
  std::cout << "phase_3" << std::endl;

  // should match our actual luminosity when parameters are configured
  // correctly.
  std::cout << "Luminosity = " << sum_prob << std::endl;
  std::cout << "done accumulating Gaussian distbns." << std::endl;

  //====================running over a no. of bunch crossings==================================
  unsigned long int cross_count = 0;
  unsigned long int event_limit_count = 0;
  int nonzero_event_count = 0;//debugging
  while((cross_count < 20000000) && (event_limit_count < (unsigned int)count_norm)) { 
    //=================== counting crossing number required to reach event limit =============
    cross_count++; 
    //====================choosing no. of events for a crossing from Poisson distbn=======================
    temp_prob_poisson = fabs(static_cast<double>(rand())/RAND_MAX);
    condp = false;
    int l = 0;
    while((condp == false) && (l < (MAX_COLL - 1))) {
      if(temp_prob_poisson <= poisson_dist[0]) {
        N_coll = 0;
        condp = true;
      } else if((temp_prob_poisson > poisson_dist[l]) && (temp_prob_poisson <= poisson_dist[l + 1])) { // this is why we set max_coll to 5 + 1
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
      rand_prob_t = sum_prob*fabs(static_cast<double>(rand())/RAND_MAX);
      condt = false;
      int find_t = 0;
      while((condt == false) && (find_t<(N_bin_t - 1))) {
        if(rand_prob_t <= t_dist[0]) {
          coll_time[coll_no] = input_time[0];
          temp_t_index[coll_no] = 0;
          condt = true;
        } else if((rand_prob_t > t_dist[find_t]) && (rand_prob_t <= t_dist[find_t + 1])) {
          coll_time[coll_no] = input_time[find_t + 1];
          temp_t_index[coll_no] = find_t + 1;
          condt = true;
        } else {
          condt = false;
        }
        find_t ++;
      }
      ztime = coll_time[coll_no];

      //=================CHOSSING Z POSITIONSFROM GAUSSIAN DISTRIBUTION=====================
      for(int cz=0; cz<N_bin_z; cz++) {
        z_dist[cz] = gaussian_dist[temp_t_index[coll_no]][cz];
      }//selecting the z-column with t-row fixed

      rand_prob_z = gaussian_dist[temp_t_index[coll_no]][0] + (z_norm[temp_t_index[coll_no]])*fabs(static_cast<double>(rand())/RAND_MAX);
      condz = false;
      int find_z = 0;
      while((condz == false) && (find_z < (N_bin_z - 1))) {
        if(rand_prob_z <= z_dist[0]) {
          coll_pos[coll_no] = position[0];
          condz = true;
        } else if((rand_prob_z > z_dist[find_z]) && (rand_prob_z <= z_dist[find_z + 1])) {
          coll_pos[coll_no] = position[find_z + 1];
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
  time_tracker["phase_4"] = GetTime().count();
  std::cout << "phase_4" << std::endl;
  std::cout << "job done" << std::endl;
  std::cout << "cross_count: " << cross_count << std::endl;
  std::cout << "event_limit_count: " << event_limit_count << std::endl;
  std::cout << "Finished profiling code." << std::endl;
  std::cout << "  raw time, phase_1: " << time_tracker["phase_1"] << std::endl;
  std::cout << "  raw time, phase_2: " << time_tracker["phase_2"] << std::endl;
  std::cout << "  raw time, phase_3: " << time_tracker["phase_3"] << std::endl;
  std::cout << "  raw time, phase_4: " << time_tracker["phase_4"] << std::endl;
  std::cout << "Milliseconds per phase: " << std::endl;
  std::cout << "  phase_1: " << time_tracker["phase_1"] - time_tracker["start"]   << " milliseconds " << "(" << (time_tracker["phase_1"]-time_tracker["start"]  )/60000.0 << " minutes)"<< std::endl;
  std::cout << "  phase_2: " << time_tracker["phase_2"] - time_tracker["phase_1"] << " milliseconds " << "(" << (time_tracker["phase_2"]-time_tracker["phase_1"])/60000.0 << " minutes)"<< std::endl;
  std::cout << "  phase_3: " << time_tracker["phase_3"] - time_tracker["phase_2"] << " milliseconds " << "(" << (time_tracker["phase_3"]-time_tracker["phase_2"])/60000.0 << " minutes)"<< std::endl;
  std::cout << "  phase_4: " << time_tracker["phase_4"] - time_tracker["phase_3"] << " milliseconds " << "(" << (time_tracker["phase_4"]-time_tracker["phase_3"])/60000.0 << " minutes)"<< std::endl;
  std::cout << "In phase_3, we did: " << how_many_things << " things." << std::endl;
  return 0;
}

int HourglassSimulation::OverrideSaveFile( const std::string& save_file) {
  override_save_file_name_ = save_file;
  override_save_file_ = true;
  return 0;
}

int HourglassSimulation::SaveFigures( const std::string& figure_output_dir = "./") {
  gErrorIgnoreLevel = kWarning;
  std::stringstream tfile_name;
  std::stringstream name;

  tfile_name << figure_output_dir << "/" << run_number_ << "_" << "h" << xoff << "_v" << yoff << "_beta" << beta_star << "_xing" << angle << "_HourglassSimulation.root"; 
  name       << figure_output_dir << "/" << run_number_ << "_" << "h" << xoff << "_v" << yoff << "_beta" << beta_star << "_xing" << angle << "_HourglassSimulation.pdf";
  
  if(override_save_file_) {
    tfile_name.str("");
    tfile_name << override_save_file_name_ << ".root";
    name.str("");
    name << override_save_file_name_ << ".pdf";
  }
  std::string out_file_name = name.str() + "[";
  TCanvas* booklet = new TCanvas("booklet","BBC Efficiency Plots");
  TFile* root_out = new TFile(tfile_name.str().c_str(), "RECREATE");
  booklet->Print(out_file_name.c_str());
  for(auto plot_i = save_registry_.begin(); plot_i != save_registry_.end(); ++plot_i) {
    auto draw_obj = *plot_i;
    if(!draw_obj) continue;
    std::cout << "Saving/Drawing: " << draw_obj->GetName() << std::endl;
    root_out->cd();
    draw_obj->Write();
    if(!draw_obj) continue;
    TCanvas* c = new TCanvas();
    c->cd();
    draw_obj->Draw();
    c->Print(name.str().c_str());
    delete c;
  }
  out_file_name = name.str() + "]";
  booklet->Print(out_file_name.c_str());
  root_out->Write();
  root_out->Close();
  if(root_out) delete root_out;
  delete booklet;
  gErrorIgnoreLevel = kInfo;
  return 0;
}
