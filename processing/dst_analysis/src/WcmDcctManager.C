#include "WcmDcctManager.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TLatex.h"
#include "TError.h"
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <set>
const int WcmDcctManager::NUMBER_OF_BUNCHES = 120;

WcmDcctManager::WcmDcctManager() {
  std::stringstream ss;
  ss << "WcmDcctManager_0x" << std::hex << this;
  this_name_ = ss.str();
  run_number_ = "RUN_NOT_SET";
  wcm_b_fname_= "FILE_NOT_SET";
  wcm_y_fname_= "FILE_NOT_SET";
  dcct_fname_= "FILE_NOT_SET";
  TIME_OFFSET = 0;
  wcm_blue_max_ = 0.0;
  wcm_yell_max_ = 0.0;
  dcct_blue_max_ = 0.0;
  dcct_yell_max_ = 0.0;
  wcm_blue_min_ = 0.0;
  wcm_yell_min_ = 0.0;
  dcct_blue_min_ = 0.0;
  dcct_yell_min_ = 0.0;
  calib_blue_min_ = 0.0;
  calib_blue_max_ = 0.0;
  calib_yell_min_ = 0.0;
  calib_yell_max_ = 0.0;
}

int WcmDcctManager::SetTimeOffset(time_t offset) {
  TIME_OFFSET = offset; 
  std::cout << "Using new time offset of " << TIME_OFFSET << std::endl;
  return 0;
}

WcmDcctManager::~WcmDcctManager() {

}

int WcmDcctManager::Init( 
    const std::string& run_number_, 
    const std::string& wcm_b_fname_, 
    const std::string& wcm_y_fname_, 
    const std::string& dcct_fname_
    ) {
  this->run_number_  = run_number_;
  this->wcm_b_fname_ = wcm_b_fname_;
  this->wcm_y_fname_ = wcm_y_fname_;
  this->dcct_fname_  = dcct_fname_   ;
  LoadDcct();
  LoadWcmBlue();
  LoadWcmYellow();
  auto begin_map = blue_pop_.begin();
  auto end_map = blue_pop_.end();
  end_map--;

  // The WCM/DCCT data approximately overlaps the vernier scan.
  TIME_OFFSET = begin_map->first;
  BEGIN_SCAN = TIME_OFFSET;
  END_SCAN = end_map->first;
  std::cout << "Using time offset of " << TIME_OFFSET << std::endl;
  std::cout << "Scan started at " << BEGIN_SCAN << " offset: " << BEGIN_SCAN-TIME_OFFSET << std::endl;
  std::cout << "Scan ended at " << END_SCAN << " offset: " << END_SCAN-BEGIN_SCAN << std::endl;
  return 0;
}

int WcmDcctManager::GenerateBounds() {
  wcm_blue_max_ = 0.0;
  dcct_blue_max_ = 0.0;
  wcm_blue_min_ = 0.0;
  dcct_blue_min_ = 0.0;
  calib_blue_min_ = 0.0;
  calib_blue_max_ = 0.0;

  wcm_yell_max_ = 0.0;
  dcct_yell_max_ = 0.0;
  wcm_yell_min_ = 0.0;
  dcct_yell_min_ = 0.0;
  calib_yell_min_ = 0.0;
  calib_yell_max_ = 0.0;

  bool first_dcct = true;
  bool first_wcm = true;
  for(auto i = blue_pop_.begin(); i != blue_pop_.end(); ++i) {
    BeamPopulation p = i->second;
    if(first_dcct) {
      dcct_blue_max_ = p.dcct_;
      dcct_blue_min_ = p.dcct_;
      first_dcct = false;
    } else {
      if(p.dcct_ > dcct_blue_max_ ) dcct_blue_max_ = p.dcct_;
      if(p.dcct_ < dcct_blue_min_ ) dcct_blue_min_ = p.dcct_;
    }
    if(p.WcmAvailable()){
      std::vector<float> wcm_sorted = p.wcm_calib_;
      std::sort(wcm_sorted.begin(),wcm_sorted.end());
      wcm_blue_min_ = *std::lower_bound(wcm_sorted.begin(),wcm_sorted.end(),10.0); // first element NOT lower than 10.0 (i.e. lowest non-empty bunch)
      wcm_blue_max_ = *std::max_element(wcm_sorted.begin(),wcm_sorted.end());
      if(fabs(p.calib_ < 0.1)) std::cout << p.calib_ << std::endl;
      if(first_wcm) {
        calib_blue_min_ = p.calib_;
        calib_blue_max_ = p.calib_;
        first_wcm = false;
      } else {
        if(p.calib_ > calib_blue_max_ ) calib_blue_max_ = p.calib_;
        if(p.calib_ < calib_blue_min_ ) calib_blue_min_ = p.calib_;
      }
    }
  }
  first_dcct = true;
  first_wcm = true;
  for(auto i = yell_pop_.begin(); i != yell_pop_.end(); ++i) {
    BeamPopulation p = i->second;
    if(first_dcct) {
      dcct_yell_max_ = p.dcct_;
      dcct_yell_min_ = p.dcct_;
      first_dcct = false;
    } else {
      if(p.dcct_ > dcct_yell_max_ ) dcct_yell_max_ = p.dcct_;
      if(p.dcct_ < dcct_yell_min_ ) dcct_yell_min_ = p.dcct_;
    }
    if(p.WcmAvailable()){
      std::vector<float> wcm_sorted = p.wcm_calib_;
      std::sort(wcm_sorted.begin(),wcm_sorted.end());
      wcm_yell_min_ = *std::lower_bound(wcm_sorted.begin(),wcm_sorted.end(),10.0); // first element NOT lower than 10.0 (i.e. lowest non-empty bunch)
      wcm_yell_max_ = *std::max_element(wcm_sorted.begin(),wcm_sorted.end());
      if(first_wcm) {
        calib_yell_min_ = p.calib_;
        calib_yell_max_ = p.calib_;
        first_wcm = false;
      } else {
        if(p.calib_ > calib_yell_max_ ) calib_yell_max_ = p.calib_;
        if(p.calib_ < calib_yell_min_ ) calib_yell_min_ = p.calib_;
      }
    }
  }
  std::cout << "wcm_blue_max_  : " << wcm_blue_max_  << "e9 ions/bunch" << std::endl;
  std::cout << "wcm_blue_min_  : " << wcm_blue_min_  << "e9 ions/bunch" << std::endl;
  std::cout << "dcct_blue_max_ : " << dcct_blue_max_ << "e11 ions in beam" <<  std::endl;
  std::cout << "dcct_blue_min_ : " << dcct_blue_min_ << "e11 ions in beam" <<  std::endl;
  std::cout << "calib_blue_min_: " << calib_blue_min_ << " DCCT/WCM_TOT" << std::endl;
  std::cout << "calib_blue_max_: " << calib_blue_max_ << " DCCT/WCM_TOT" << std::endl;
  std::cout << "wcm_yell_max_  : " << wcm_yell_max_  << "e9 ions/bunch" << std::endl;
  std::cout << "wcm_yell_min_  : " << wcm_yell_min_  << "e9 ions/bunch" << std::endl;
  std::cout << "dcct_yell_max_ : " << dcct_yell_max_ << "e11 ions in beam" <<  std::endl;
  std::cout << "dcct_yell_min_ : " << dcct_yell_min_ << "e11 ions in beam" <<  std::endl;
  std::cout << "calib_yell_min_: " << calib_yell_min_ << " DCCT/WCM_TOT" << std::endl;
  std::cout << "calib_yell_max_: " << calib_yell_max_ << " DCCT/WCM_TOT" << std::endl;
  return 0;
}

int WcmDcctManager::InitPlots(const int beam_choice) {
  std::string beam_name = "";
  std::string obj_name = "";
  std::stringstream name;
  std::stringstream title;
  if(beam_choice != 0 && beam_choice != 1 ) {
    std::cerr << " Invalid choice for " << this_name_ << "::InitPlots. Choose 0 = blue beam, or 1 = yellow beam" << std::endl;
    return 1;
  }
  // Set up storage containers to point in the right place
  auto& wcm_dist         = (beam_choice == 0) ? wcm_dist_blue_ : wcm_dist_yell_;
  auto& wcm_dist_vs_time = (beam_choice == 0) ? wcm_dist_vs_time_blue_ : wcm_dist_vs_time_yell_;
  TGraph* calib_vs_time;
  TH1F*   calib;
  TGraph* dcct_vs_time;
  TGraph* wcm_tot_vs_time;
  float   wcm_min;
  float   wcm_max;
  float   dcct_min;
  float   dcct_max;
  float   calib_min;
  float   calib_max;
  // Setup Names, Domains
  switch(beam_choice) {
    case 0:
      beam_name = "Blue Beam";
      obj_name = "blue";
      wcm_min = wcm_blue_min_;
      wcm_max = wcm_blue_max_;
      dcct_min = dcct_blue_min_;
      dcct_max = dcct_blue_max_;
      calib_min = calib_blue_min_;
      calib_max = calib_blue_max_;
      break;
    case 1:
      beam_name = "Yellow Beam";
      obj_name = "yell";
      wcm_min = wcm_yell_min_;
      wcm_max = wcm_yell_max_;
      dcct_min = dcct_yell_min_;
      dcct_max = dcct_yell_max_;
      calib_min = calib_yell_min_;
      calib_max = calib_yell_max_;
      break;
    default:
      // Should never get here.
      std::cerr << " Invalid choice for " << this_name_ << "::InitPlots. Choose 0 = blue beam, or 1 = yellow beam" << std::endl;
      break;
  }
  for(int i = 0; i < 120; i++) { 
    name.str("");
    title.str("");
    name << "wcm_dist_" << obj_name << "_bunch" << i;
    title << "WCM distribution for Bunch:" << i << ";Bunch Population(#times 10^{9});counts";
    wcm_dist[i] = new TH1F(name.str().c_str(),title.str().c_str(),30,wcm_min,wcm_max);
    bunch_registry_.push_back(wcm_dist[i] );

    name.str("");
    title.str("");
    name << "wcm_dist_vs_time_" << obj_name << "_bunch" << i;
    title << "WCM Population, " << beam_name << " Bunch:" << i << " vs Time;time(s);population(#times 10^{9})";
    wcm_dist_vs_time[i] = new TGraph();
    wcm_dist_vs_time[i]->SetName(name.str().c_str());
    wcm_dist_vs_time[i]->SetTitle(title.str().c_str());
    bunch_registry_.push_back(wcm_dist_vs_time[i] );
  }

  // Initialize names, etc for everything 
  calib_vs_time = new TGraph();
  name.str("");
  title.str("");
  name << "calib_vs_time_" << obj_name;
  title << "Calibration vs Time;time(s);DCCT_SUM/WCM_SUM";
  calib_vs_time->SetName(name.str().c_str());
  calib_vs_time->SetTitle(title.str().c_str());
  save_registry_.push_back(calib_vs_time );

  name.str("");
  title.str("");
  name << "calib_" << obj_name;
  title << "Calibration Constant Distribution, " << beam_name << ";DCCT_SUM/WCM_TOT;Counts";
  calib = new TH1F(name.str().c_str(),title.str().c_str(),30,calib_min,calib_max);
  calib->SetName(name.str().c_str());
  calib->SetTitle(name.str().c_str());
  save_registry_.push_back(calib );

  name.str("");
  title.str("");
  name << "dcct_" << obj_name << "_vs_time";
  title << "DCCT Population vs Time, " << beam_name << ";time(s);Population(#times 10^{11})";
  dcct_vs_time = new TGraph();
  dcct_vs_time->SetName(name.str().c_str());
  dcct_vs_time->SetTitle(title.str().c_str());
  save_registry_.push_back(dcct_vs_time );

  name.str("");
  title.str("");
  name << "wcm_tot_vs_time_" << obj_name;
  title << "WCM Total vs Time, " << beam_name << ";time(s);Population(#times 10^{11})"; 
  wcm_tot_vs_time = new TGraph();
  wcm_tot_vs_time->SetName(name.str().c_str());
  wcm_tot_vs_time->SetTitle(title.str().c_str());
  save_registry_.push_back(wcm_tot_vs_time );

  // Now that we have actual memory allocated to pointers, we can point them to the
  // right place.
  switch(beam_choice) {
    case 0:
      calib_vs_time_blue_ = calib_vs_time;
      calib_blue_ = calib;
      dcct_vs_time_blue_ = dcct_vs_time;
      wcm_tot_vs_time_blue_ = wcm_tot_vs_time;
      break;
    case 1:
      calib_vs_time_yell_ = calib_vs_time;
      calib_yell_ = calib;
      dcct_vs_time_yell_ = dcct_vs_time;
      wcm_tot_vs_time_yell_ = wcm_tot_vs_time;
      break;
    default:
      std::cerr << " Invalid choice for " << this_name_ << "::InitPlots. Choose 0 = blue beam, or 1 = yellow beam" << std::endl;
      break;
  }

  // remove this when we finally use dcct_min/dcct_max for sometihng
  std::cout << "DCCT range: " << dcct_min << ", " << dcct_max << std::endl;
  return 0;
}

int WcmDcctManager::SetScanRange(time_t begin, time_t end) {
  BEGIN_SCAN = begin;
  END_SCAN = end;
  std::cout << "Updated: Using time offset of " << TIME_OFFSET << std::endl;
  std::cout << "Scan started at " << BEGIN_SCAN << " offset: " << BEGIN_SCAN-TIME_OFFSET << std::endl;
  std::cout << "Scan ended at " << END_SCAN << " offset: " << END_SCAN-BEGIN_SCAN << std::endl;
  return 0;
}

int WcmDcctManager::Run() {
  GenerateBounds();
  int blue_beam = 0;
  int yellow_beam = 1;
  InitPlots(blue_beam);
  InitPlots(yellow_beam);
  MakePlots(blue_beam);
  MakePlots(yellow_beam);
  FitRateLoss();
  CalculateLuminosityLoss();
  GetSummaryStatistics();
  return 0;
}

int WcmDcctManager::LoadDcct() {
  std::ifstream inFile(dcct_fname_.c_str());
  if(inFile) {
    std::string line = "";
    // Load only DCCT data, ignore the WCM summed data.
    while(getline(inFile, line)) {
      if(line[0] == '#') continue;
      time_t c_time;
      BeamPopulation blue;
      BeamPopulation yell;
      blue.Reset();
      yell.Reset();
      std::stringstream ss;
      ss.str(line);
      ss >> c_time >> blue.dcct_ >> yell.dcct_;
      blue_pop_[c_time] = blue;
      yell_pop_[c_time] = yell;
    }
  } else {
    std::cerr << "DCCT FILE NOT FOUND: " << dcct_fname_ << std::endl;
  }
  return 0; 
}

int WcmDcctManager::GetSummaryStatistics() {
  std::cout << "WCM/DCCT Summary: " << std::endl;
  float average_blue_bunch_pop = 0.;
  float average_yell_bunch_pop = 0.;
  float n_filled_blue = 0.;
  float n_filled_yell = 0.;
  std::set<int> empty_bunches;
  for(int i = 0; i < 120; i++) {
    float bunch_blue = wcm_dist_blue_[i]->GetMean();
    float bunch_yell = wcm_dist_yell_[i]->GetMean();
  // Diagnostic Output
  // std::cout << "Bunch_" << std::setw(3) << i 
  //         << " blue: " << std::setw(8) << bunch_blue 
  //         << " yell: " << std::setw(8) << bunch_yell 
  //         << " probable empty crossing (if zero): " << std::setw(8) << bunch_blue * bunch_yell << std::endl;
    if(fabs(bunch_blue) >  5.0) {
      average_blue_bunch_pop += bunch_blue;
      n_filled_blue++;
    } else {
      empty_bunches.insert(i);
    }
    if(fabs(bunch_yell) > 5.0) {
      average_yell_bunch_pop += bunch_yell;
      n_filled_yell++;
    } else {
      empty_bunches.insert(i);
    }
  }
  std::cout << "Blue or Yellow Empty Bunches cause effective empty bunch crossing for bunches: ";
  for(auto i = empty_bunches.begin(); i!=empty_bunches.end(); ++i) {
    std::cout << *i << ", ";
  }
  std::cout << std::endl;
  std::cout << "Average blue beam bunch population: " << average_blue_bunch_pop/n_filled_blue << "e9 ions" << std::endl;
  std::cout << "Filled Bunches Blue: " << n_filled_blue << std::endl;
  std::cout << "Average yell beam bunch population: " << average_yell_bunch_pop/n_filled_yell << "e9 ions" << std::endl;
  std::cout << "Filled Bunches Yell: " << n_filled_yell << std::endl;
  std::cout << "Effective filled bunches: " << 120 - empty_bunches.size() << std::endl;

  blue_beam_population_ = average_blue_bunch_pop/n_filled_blue;
  yellow_beam_population_ = average_yell_bunch_pop/n_filled_yell;
  filled_bunch_crossings_ = 120-empty_bunches.size();

  return 0;
}

int WcmDcctManager::FitRateLoss() {
  rate_loss_ = new TF1("rate_loss_","pol1(0)",(float)(BEGIN_SCAN-TIME_OFFSET),(float)(END_SCAN-TIME_OFFSET));
  save_registry_.push_back(rate_loss_ );
  wcm_product_vs_time_->Fit(rate_loss_);
  return 0;
}

int WcmDcctManager::MakePlots(const int beam_choice) {
  std::string beam_name = "";
  std::string obj_name = "";
  std::stringstream name;
  std::stringstream title;
  if(beam_choice != 0 && beam_choice != 1 ) {
    std::cerr << " Invalid choice for " << this_name_ << "::InitPlots. Choose 0 = blue beam, or 1 = yellow beam" << std::endl;
    return 1;
  }
  // Set up storage containers to point in the right place
  auto& wcm_dist         = (beam_choice == 0) ? wcm_dist_blue_ : wcm_dist_yell_;
  auto& wcm_dist_vs_time = (beam_choice == 0) ? wcm_dist_vs_time_blue_ : wcm_dist_vs_time_yell_;
  auto& beam_population  = (beam_choice == 0) ? blue_pop_ : yell_pop_;
  TGraph* calib_vs_time;
  TH1F*   calib;
  TGraph* dcct_vs_time;
  TGraph* wcm_tot_vs_time;
  float   wcm_min;
  float   wcm_max;
  float   dcct_min;
  float   dcct_max;
  switch(beam_choice) {
    case 0:
      calib_vs_time = calib_vs_time_blue_;
      calib = calib_blue_;
      dcct_vs_time = dcct_vs_time_blue_;
      wcm_tot_vs_time = wcm_tot_vs_time_blue_;
      beam_name = "Blue Beam";
      obj_name = "blue";
      wcm_min = wcm_blue_min_;
      wcm_max = wcm_blue_max_;
      dcct_min = dcct_blue_min_;
      dcct_max = dcct_blue_max_;
      break;
    case 1:
      calib_vs_time = calib_vs_time_yell_;
      calib = calib_yell_;
      dcct_vs_time = dcct_vs_time_yell_;
      wcm_tot_vs_time = wcm_tot_vs_time_yell_;
      beam_name = "Yellow Beam";
      obj_name = "yell";
      wcm_min = wcm_yell_min_;
      wcm_max = wcm_yell_max_;
      dcct_min = dcct_yell_min_;
      dcct_max = dcct_yell_max_;
      break;
    default:
      std::cerr << " Invalid choice for " << this_name_ << "::InitPlots. Choose 0 = blue beam, or 1 = yellow beam" << std::endl;
      break;
  }
  bool wcm_product_exists = false;
  if(!wcm_product_vs_time_){
    wcm_product_vs_time_ = new TGraph();
    wcm_product_vs_time_->SetName("wcm_product_vs_time");
    wcm_product_vs_time_->SetTitle("WCM_{blue total} #times WCM_{yellow total} vs time;time(s);Product");
    save_registry_.push_back(wcm_product_vs_time_);
  } else {
    wcm_product_exists = true;
  }
  // wp remove this when finally using these variables
  std::cout << "wcm_min  :" << wcm_min << std::endl;
  std::cout << "wcm_max  :" << wcm_max << std::endl;
  std::cout << "dcct_min :" << dcct_min << std::endl;
  std::cout << "dcct_max :" << dcct_max << std::endl;
  
  for(auto i = beam_population.begin(); i != beam_population.end(); ++i) {
    time_t t = i->first;
    BeamPopulation p = i->second;
    time_t t_rel = t - BEGIN_SCAN;
    dcct_vs_time->SetPoint(dcct_vs_time->GetN(), t_rel, p.dcct_);
    if(p.WcmAvailable()){ 
      wcm_tot_vs_time->SetPoint(wcm_tot_vs_time->GetN(), t_rel, p.wcm_tot_calib_);
      calib_vs_time->SetPoint(calib_vs_time->GetN(), t_rel, p.calib_);
      calib->Fill(p.calib_);
      if(wcm_product_exists) { // if we just created the object, this is true
        auto t_blue_lookup = blue_pop_.find(t);
        auto t_yell_lookup = yell_pop_.find(t);
        if( (t_blue_lookup != blue_pop_.end() ) && ( t_yell_lookup != yell_pop_.end() ) ){
          auto blue = t_blue_lookup->second;
          auto yell = t_yell_lookup->second;
          wcm_product_vs_time_->SetPoint(wcm_product_vs_time_->GetN(), t_rel, blue.wcm_tot_calib_ * yell.wcm_tot_calib_);
        }
      }
      for(int bunch_i = 0; bunch_i < 120; bunch_i++){
        wcm_dist[bunch_i]->Fill(p.wcm_calib_[bunch_i]);
        wcm_dist_vs_time[bunch_i]->SetPoint(wcm_dist_vs_time[bunch_i]->GetN(),t_rel,p.wcm_calib_[bunch_i]);
      }
    }
  }
  return 0;
}

int WcmDcctManager::LoadWcmBlue() {
  std::ifstream inFile(wcm_b_fname_.c_str());
  if(inFile) {
    std::string line = "";
    while(getline(inFile, line)) {
      if(line[0] == '#') continue;
      time_t c_time;
      std::vector<float> wcm;
      std::stringstream ss;
      ss.str(line);
      ss >> c_time;
      double population;
      while(ss >> population) {
        wcm.push_back(population);
      }
      if(wcm.size() != 120){ 
        std::cerr << "Could not parse bunch population from file!! Check Parsing for class defintion of " << this_name_ << std::endl;
        std::cerr << "WCM BLUE BEAM FILE PARSE ERROR: " << wcm_b_fname_ << std::endl;
        return 1;
      }
      auto search = blue_pop_.find(c_time);
      if( search != blue_pop_.end() ) {
        (*search).second.wcm_ = wcm;
        (*search).second.SumWcm();
        (*search).second.Calibrate();
      } else {
        std::cerr << "Did not find time index : " << c_time << " in " << wcm_b_fname_ << std::endl;
        std::cerr << "skipping, but if this happens too often, you may have mistmached dcct or wcm files." << std::endl;
      }
    }
  } else {
    std::cerr << "WCM BLUE BEAM FILE NOT FOUND: " << wcm_b_fname_ << std::endl;
  }
  return 0;
}

int WcmDcctManager::LoadWcmYellow() {
  std::ifstream inFile(wcm_y_fname_.c_str());
  if(inFile) {
    std::string line = "";
    while(getline(inFile, line)) {
      if(line[0] == '#') continue;
      time_t c_time;
      std::vector<float> wcm;
      std::stringstream ss;
      ss.str(line);
      ss >> c_time;
      double population;
      while(ss >> population) {
        wcm.push_back(population);
      }
      if(wcm.size() != 120){ 
        std::cerr << "Could not parse bunch population from file!! Check Parsing for class defintion of " << this_name_ << std::endl;
        std::cerr << "WCM YELLOW BEAM FILE PARSE ERROR: " << wcm_y_fname_ << std::endl;
        return 1;
      }
      auto search = yell_pop_.find(c_time);
      if( search != yell_pop_.end() ) {
        (*search).second.wcm_ = wcm;
        (*search).second.SumWcm();
        (*search).second.Calibrate();
      } else {
        std::cerr << "Did not find time index : " << c_time << " in " << wcm_y_fname_ << std::endl;
        std::cerr << "skipping, but if this happens too often, you may have mistmached dcct or wcm files." << std::endl;
      }
    }
  } else {
    std::cerr << "WCM YELLOW BEAM FILE NOT FOUND: " << wcm_y_fname_ << std::endl;
  }
  return 0;
}

void WcmDcctManager::PrintBlueDCCT() {
  for(auto i = blue_pop_.begin(); i != blue_pop_.end(); ++i) {
    auto t = i->first;
    auto p = i->second;
    std::cout << t << ": " << p.dcct_ << ", ";
  }
  std::cout << std::endl;
}

void WcmDcctManager::PrintYellowDCCT() {
  for(auto i = yell_pop_.begin(); i != yell_pop_.end(); ++i) {
    auto t = i->first;
    auto p = i->second;
    std::cout << t << ": " << p.dcct_ << ", ";
  }
  std::cout << std::endl;
}

void WcmDcctManager::PrintBlueWcmTotal() {
  for(auto i = blue_pop_.begin(); i != blue_pop_.end(); ++i) {
    auto t = i->first;
    auto p = i->second;
    if(p.wcm_tot_ > 0) std::cout << t << ": " << p.wcm_tot_ << ", ";
  }
  std::cout << std::endl;
}

void WcmDcctManager::PrintYellowWcmTotal() {
  for(auto i = yell_pop_.begin(); i != yell_pop_.end(); ++i) {
    auto t = i->first;
    auto p = i->second;
    if(p.wcm_tot_ > 0) std::cout << t << ": " << p.wcm_tot_ << ", ";
  }
  std::cout << std::endl;
}

bool WcmDcctManager::CheckIfTimeExists(time_t time) {
  auto dcct_blue_find = blue_pop_.find(time);
  if( dcct_blue_find == blue_pop_.end() ) return false;
  else return true;
}

int WcmDcctManager::PrintDataTimeIndex(time_t time) {
  std::cout << "Data at time " << time << std::endl; 
  auto find_blue = blue_pop_.find(time);
  if( find_blue == blue_pop_.end() ) {
    std::cerr << " Time index not found in DCCT blue. Check if it is in your files.\n"; 
    return 1;
  }
  auto find_yell = yell_pop_.find(time);
  if( find_yell == yell_pop_.end() ) {
    std::cerr << " Time index not found in DCCT yellow. Check if it is in your files.\n"; 
    return 1;
  }
  BeamPopulation p_blue = find_blue->second; 
  BeamPopulation p_yell = find_yell->second; 
  std::cout << "DCCT Blue    Population: " << p_blue.dcct_    << std::endl;
  std::cout << "DCCT Yellow  Population: " << p_yell.dcct_    << std::endl;
  std::cout << "WcmB Bunched Population: " << p_blue.wcm_tot_ << std::endl;
  std::cout << "WcmY Bunched Population: " << p_yell.wcm_tot_ << std::endl;

  int bunch_counter = 0;
  std::cout << "  Blue Beam: " << std::endl;
  for(auto i = p_blue.wcm_.begin() ; i != p_blue.wcm_.end(); ++i) {
    if(bunch_counter%12 == 0 ) std::cout << "    ";
    std::cout << *i << ", ";
    if((bunch_counter+1)%12 == 0) std::cout << std::endl;
    bunch_counter++;
  }

  bunch_counter = 0;
  std::cout << "  Yellow Beam: " << std::endl;
  for(auto i = p_yell.wcm_.begin(); i != p_yell.wcm_.end(); ++i) {
    if(bunch_counter%12 == 0 ) std::cout << "    ";
    std::cout << *i << ", ";
    if((bunch_counter+1)%12 == 0) std::cout << std::endl;
    bunch_counter++;
  }
  return 0;
}

int WcmDcctManager::CalculateLuminosityLoss() {
  if(!rate_loss_) return 1;
  luminosity_loss_ = ((rate_loss_->Eval(BEGIN_SCAN-TIME_OFFSET)) - rate_loss_->Eval(END_SCAN-TIME_OFFSET)) / rate_loss_->Eval(BEGIN_SCAN-TIME_OFFSET);
  return 0;
}

double WcmDcctManager::GetLuminosityLoss() {
  return luminosity_loss_;
}

int WcmDcctManager::SaveFigures(const std::string& figure_output_dir) {
  gErrorIgnoreLevel = kWarning;
  std::string tfile_name = figure_output_dir + "/" + run_number_ + "_WcmDcctPlots.root"; 
  std::string pdf_name   = figure_output_dir + "/" + run_number_ + "_WcmDcctPlots.pdf";
  std::string out_file_name = pdf_name + "[";
  TCanvas* booklet = new TCanvas("booklet","Scaler Rate Plots");
  TFile* root_out = new TFile(tfile_name.c_str(), "RECREATE");
  booklet->Print(out_file_name.c_str());
  for(auto plot_i = save_registry_.begin(); plot_i != save_registry_.end(); ++plot_i) {
    auto draw_obj = *plot_i;
    //std::string name = draw_obj->GetName();
    if(!draw_obj) continue;
    root_out->cd();
    draw_obj->Write();
    TCanvas* c = new TCanvas();
    c->cd();
    draw_obj->Draw();
    c->Print(pdf_name.c_str()); // ->GetName() returns const char*
    delete c;
  }
  TDirectory* dir_bunches = root_out->mkdir("bunches");
  for(auto i = bunch_registry_.begin(); i != bunch_registry_.end(); ++i) {
    auto draw_obj = *i;
    if(!draw_obj) continue;
    dir_bunches->cd();
    draw_obj->Write();
    TCanvas* c = new TCanvas();
    c->cd();
    draw_obj->Draw();
    c->Print(pdf_name.c_str());
    delete c;
  }
  out_file_name = pdf_name + "]";
  booklet->Print(out_file_name.c_str());
  root_out->Write();
  root_out->Close();
  delete booklet;
  gErrorIgnoreLevel = kInfo;
  return 0;
}

int WcmDcctManager::SaveBeamPopulations(const std::string& out_dir) {
  std::stringstream out_file_name;
  out_file_name << out_dir << "/" << run_number_ << "_WCMDCCT_BeamPopulation.txt";
  std::ofstream out_file(out_file_name.str().c_str());
  out_file << "AVG_NUMBER_IONS_BLUE_BEAM "   << blue_beam_population_   << "e9" << std::endl;
  out_file << "AVG_NUMBER_IONS_YELLOW_BEAM " << yellow_beam_population_ << "e9" << std::endl;
  out_file << "FILLED_BUNCHES " << filled_bunch_crossings_ << std::endl;
  return 0;
}
