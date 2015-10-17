#include "HourglassSimulation.h"

#include<iostream>
#include <cstdlib>
#include<fstream>
#include<cmath>
#include<ctime>
#include<math.h>
#include <chrono> 
#include <map>
#include <vector> 
#include <sstream>
#include <iostream>
#include <sstream>
#include <iostream>

#include "TF1.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"

#define PI 3.14159265

// Local function definition to avoid dealing with ROOT's immense stupidity
std::chrono::milliseconds GetTime() {
  std::chrono::milliseconds the_time;
  the_time = std::chrono::duration_cast < std::chrono::milliseconds >(std::chrono::system_clock::now().time_since_epoch());
  return the_time;
}

HourglassSimulation::HourglassSimulation() {
  std::stringstream ss;
  ss << "HourglassSimulation_0x" << std::hex << this;
  this_name = ss.str();
  std::cout << "Instantiating " << this_name << std::endl;
}

HourglassSimulation::~HourglassSimulation() {
  std::cout << "Destroying: " << this_name << std::endl;
  for(auto i = save_registry_.begin(); i != save_registry_.end(); ++i){
    if(*i) delete *i;
  }
}

int HourglassSimulation::Init() {
  /** scan configuration */
  run_number_ = "359711";
  xoff = 0.1; 
  yoff = 0.0;  
  count_norm = 592;
  rate = 0.001; 
  max_coll = 5 + 1;

  /** variables for z-t ~Gaussian distribution of intial beam profile */
  n_bunch = 109;
  freq = 78213.0;
  scale = 1.5; // What is this for..? wp
  Np = 120.0e9; // blue
  Nm = 88.0e9;  // yellow
  sigma_xstar = 0.0243/sqrt(2.0); // From run 12 scan (no weighting on beam width)
  sigma_ystar = 0.0237/sqrt(2.0); 
  sigma_zl = 35.15*scale; // need to use WCM profile directly.
  sigma_zc = 27.65*scale; 
  sigma_zr = 55.95*scale;  
  mu_zl = -70.2*scale ; 
  mu_zr =  56.7*scale; 
  z_vtx_off = 9.38; // wp was -35.25, maybe use difference at maxmial overlap for everything??
  beta_star = 85; // wp was 85
  angle = 8.e-5; // wp was -0.08e-3

  /** Defining spatial corrdinates, to be used in arrays */
  z_low = -300.0;
  z_high = 300.0;
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
  binsizeG = z_range/N_bin_z;
  binsizeX = x_range/N_bin_x;
  binsizeY = y_range/N_bin_y;
  binsizeT = t_range/N_bin_t;
  /** End spatial discreetization */
  
  zvertex = new TH1F("zvertex_simulation","Z Vertex ZDC Profile Simulation;z vertex;counts",50,z_low,z_high);
  save_registry_.push_back(zvertex);
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
  // store various timestamps through out the execuation of the code, for 
  // bottleneck tracking.
  std::map<std::string, long long unsigned int > time_tracker;
  time_tracker["start"] = GetTime().count();
  long long unsigned int how_many_things = 0; 

  std::stringstream zvtx_histo_data;
  zvtx_histo_data << "MakeHisto.txt";

  double poisson_dist[max_coll];
  std::cout << "multiple collisions rate = " << rate << std::endl;
  std::cout << "zdc accumulated counts   = " << count_norm << std::endl;
  std::cout << "horizontal beam offset   = " << xoff << std::endl;
  std::cout << "vertical beam offset     = " << yoff << std::endl;
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
  for(int no_count = 0; no_count < max_coll; ++no_count) { 
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
  //====================== creating an array with Gaussian Distribution in 'z' and 't' ===================================
  for(int t_ct = 0; t_ct < N_bin_t; t_ct++) {
    input_time[t_ct] = (t_low + static_cast<double>(t_ct)*binsizeT + binsizeT/2.0);
  }
  for(int z_ct = 0; z_ct < N_bin_z; z_ct++) {
    position[z_ct] = (z_low + static_cast<double>(z_ct)*binsizeG + binsizeG/2.0); 
    // check if I understand binning
    // std::cout << z_ct << ": " << z_low << " + " << z_ct << " * " << binsizeG << " + " << binsizeG/2.0 << " = " << position[z_ct] << std::endl;
  }
  for(int x_ct = 0; x_ct < N_bin_x; x_ct++) {
    xposition[x_ct] = (x_low + static_cast<double>(x_ct)*binsizeX + binsizeX/2.0); 
  }
  for(int y_ct = 0; y_ct < N_bin_y; y_ct++) {
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
      norm = ((n_bunch*freq*Np*Nm)/pow(2*PI, 1.5))/pow(sigma_xz*sigma_yz, 2.0);
      double spacetime_norm = norm*binsizeX*binsizeY*binsizeG*vel*binsizeT;
      for(int cx=0; cx<N_bin_x; cx++) {
        double ex_1l_term1 = exp(-0.5*pow((xposition[cx]*cos_half_angle-xoff+half_angle*position[cz])/sigma_xz, 2.0));
        double ex_2l_term1 = exp(-0.5*pow((xposition[cx]*cos_half_angle-half_angle*position[cz])/sigma_xz, 2.0));
        for(int cy=0; cy<N_bin_y; cy++) {    
          /* ALL COMPUTATION TIME SPENT HERE */
          // What have I checked to speed up?
          // 1. I replaced the computation with a root TF1 object. This actuall slowed it down.
          // 2. Are any expressions in here getting evaluated more often than they need to?
          //   cos_half_angle is known before this loop, so we can evaluate it outside.
          //   We have split the ex_1l equations up to get evaluated separately only when needed.
          //   THis will be a good test to see if the compiler was optimizating this part anyway.
          //
          //   This halfed the execution time.

          double ex_term2_common = exp(-0.5*pow(yposition[cy]/sigma_yz,2.0));
          ex_1l = ex_1l_term1
              *ex_term2_common
              *1.522*exp(-0.5*(pow((position[cz]*cos_half_angle-mu_zl-vel*input_time[ct])/sigma_zl,2)))
              *(1.0/sigma_zl);//corrected for rotaion in x-z plane, 2015
          ex_1c = ex_1l_term1
              *ex_term2_common
              *2.157*exp(-0.5*(pow((position[cz]*cos_half_angle-vel*input_time[ct])/sigma_zc,2)))
              *(1.0/sigma_zc);//corrected for rotaion in x-z plane, 2015
          ex_1r = ex_1l_term1
              *ex_term2_common
              *1.999*exp(-0.5*(pow((position[cz]*cos_half_angle-mu_zr-vel*input_time[ct])/sigma_zr,2)))
              *(1.0/sigma_zr);
          // Integrand second bunch
          ex_2l = ex_2l_term1
              *ex_term2_common
              *1.522*exp(-0.5*(pow((position[cz]*cos_half_angle+mu_zl+vel*input_time[ct])/sigma_zl,2)))
              *(1.0/sigma_zl);
          ex_2c = ex_2l_term1
              *ex_term2_common
              *2.157*exp(-0.5*(pow((position[cz]*cos_half_angle+vel*input_time[ct])/sigma_zc,2)))
              *(1.0/sigma_zc);
          ex_2r = ex_2l_term1
              *ex_term2_common
              *1.999*exp(-0.5*(pow((position[cz]*cos_half_angle+mu_zr+vel*input_time[ct])/sigma_zr,2)))
              *(1.0/sigma_zr);                       
          g = (spacetime_norm*(ex_1l+ex_1c+ex_1r)*(ex_2l+ex_2c+ex_2r));
          /* END EXPENSIVE PART */

          if(g >= 0.0) {
            sum_prob += g;
            add_T += g;
          } else {
            std::cout << "Gaussian prob. negative. Check code." <<  std::endl;
            return 0;
          }
          how_many_things++;
        }//end loop on y
      }//end loop on x
      gaussian_dist[ct][cz] = sum_prob;//storing prob for 2-D z-t grid summed over x,y
    }//end loop on z  
    z_norm[ct] = add_T;
    sum_T += add_T;
    t_dist[ct] = sum_T;//storing prob in t with z-prob summed over
  }//end loop on t
  time_tracker["phase_3"] = GetTime().count();
  std::cout << "phase_3" << std::endl;
  std::cout << "Luminosity = " << sum_prob << std::endl;//debug
  std::cout << "done accumulating Gaussian distbns." << std::endl;

  //====================running over a no. of bunch crossings==================================
  unsigned long int cross_count = 0;
  unsigned long int event_limit_count = 0;
  int nonzero_event_count = 0;//debugging
  std::ofstream raw_zvtx_out(zvtx_histo_data.str().c_str());
  while((cross_count < 20000000) && (event_limit_count < (unsigned int)count_norm)) { 
    //=================== counting crossing number required to reach event limit =============
    cross_count++; 
    //====================choosing no. of events for a crossing from Poisson distbn=======================
    temp_prob_poisson = fabs(static_cast<double>(rand())/RAND_MAX);
    condp = false;
    int l = 0;
    while((condp == false) && (l < (max_coll - 1))) {
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
      raw_zvtx_out << smeared_zpos+z_vtx_off << std::endl;
      zvertex->Fill(smeared_zpos+z_vtx_off);
    } else {
      std::cout << "could not get smeared posn, check code." << std::endl;
    }
  }// end loop over no. of events
  time_tracker["phase_4"] = GetTime().count();
  std::cout << "phase_4" << std::endl;
  std::cout << "job done" << std::endl;
  raw_zvtx_out.close();
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

int HourglassSimulation::SaveFigures( const std::string& figure_output_dir = "./") {
  std::stringstream tfile_name;
  std::stringstream name;

  tfile_name << figure_output_dir << "/" << run_number_ << "_" << "h" << xoff << "_v" << yoff << "_beta" << beta_star << "_xing" << angle << "_HourglassSimulation.root"; 
  name       << figure_output_dir << "/" << run_number_ << "_" << "h" << xoff << "_v" << yoff << "_beta" << beta_star << "_xing" << angle << "_HourglassSimulation.pdf";

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
  return 0;
}
