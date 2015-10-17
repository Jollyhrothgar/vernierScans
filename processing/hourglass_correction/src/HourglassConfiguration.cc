#include "HourglassConfiguration.h"

#include <iostream>
#include <sstream>
#include <iostream>
#include <fstream>
#include <map>

HourglassConfiguration::HourglassConfiguration() {
  std::stringstream ss;
  ss << "HourglassConfiguration_0x" << std::hex << this;
  this_name = ss.str();
  GenerateEmptyConfigFile();
  std::cout << "Instantiating " << this_name << std::endl;
  GenerateEmptyConfigFile();
}

HourglassConfiguration::~HourglassConfiguration() {
  std::cout << "Destroying: " << this_name << std::endl;
}
// determine the minimum number of libraries needed to fill params
int HourglassConfiguration::GenerateEmptyConfigFile() {
  config_["RUN_NUMBER"]                     = ""; // BeamWidthTime
  config_["X_OFFSET"]                       = ""; // BeamWidthTime 
  config_["Y_OFFSET"]                       = ""; // BeamWidthTime  
  config_["ZDC_COUNTS"]                     = ""; // HourglassData
  config_["AVG_NUMBER_IONS_BLUE_BEAM"]      = ""; // WcmDcct
  config_["AVG_NUMBER_IONS_YELLOW_BEAM"]    = ""; // WcmDcct
  config_["HORIZONTAL_BEAM_WIDTH"]          = ""; // BeamWidthTime
  config_["VERTICAL_BEAM_WIDTH"]            = ""; // BeamWidthTime
  config_["BBC_ZDC_Z_VERTEX_OFFSET"]        = ""; // HourglassData
  config_["BETA_STAR" ]                     = ""; // GOAL: find this
  config_["CROSSING_ANGLE_XZ"]              = ""; // GOAL: find this 
  config_["FILLED_BUNCHES"]                 = ""; // Fixed
  config_["BUNCH_CROSSING_FREQUENCY"]       = ""; // Fixed
  config_["Z_PROFILE_SCALE_VALUE"]          = ""; // Fixed
  config_["MULTIPLE_COLLISION_RATE"]        = ""; // wp SetDefaultValues
  config_["MAX_COLLISIONS"]                 = ""; // wp SetDefaultValues
  config_["Z_BUNCH_WIDTH_LEFT_GAUSSIAN"]    = ""; // wp SetDefaultValues
  config_["Z_BUNCH_WIDTH_RIGHT_GAUSIAN"]    = ""; // wp SetDefaultValues
  config_["Z_BUNCH_WIDTH_CENTRAL_GAUSIAN"]  = ""; // wp SetDefaultValues
  config_["Z_BUNCH_WIDTH_LEFT_OFFSET"]      = ""; // wp SetDefaultValues
  config_["Z_BUNCH_WIDTH_RIGHT_OFFSET"]     = ""; // wp SetDefaultValues
  return 0;
}

int HourglassConfiguration::SetDefaultValues() {
  config_["NUMBER_OF_BUNCHES"]              = "109";
  config_["BUNCH_CROSSING_FREQUENCY"]       = "78213.";
  config_["Z_PROFILE_SCALE_VALUE"]          = "1.5";
  config_["MAX_COLLISIONS"]                 = "5.";
  config_["BETA_STAR"]                      = "85.";
  config_["CROSSING_ANGLE_XZ"]              = "8.0e-5";
  config_["MULTIPLE_COLLISION_RATE"]        = "0.001";
  config_["MAX_COLLISIONS"]                 = "5";
  config_["Z_BUNCH_WIDTH_LEFT_GAUSSIAN"]    = "36.15";
  config_["Z_BUNCH_WIDTH_CENTRAL_GAUSIAN"]  = "27.65";
  config_["Z_BUNCH_WIDTH_RIGHT_GAUSIAN"]    = "55.95";
  config_["Z_BUNCH_WIDTH_LEFT_OFFSET"]      = "-70.2";
  config_["Z_BUNCH_WIDTH_RIGHT_OFFSET"]     = "56.7";
  return 0;
}

int HourglassConfiguration::ModifyConfigParameter(const std::string& name, const std::string& val) {
  auto search = config_.find(name);
  if(search != config_.end()){
    config_[name] = val;
  } else {
    std::cout << "Unrecognized configuration parameter, remember parameters are case-senstive: " << name << std::endl;
  }
  return 0;
}

int HourglassConfiguration::SaveConfigFile(const std::string& out_dir ) {
  std::stringstream out_file_name;
  out_file_name << out_dir << "/" << GetConfigName();
  std::ofstream out_file(out_file_name.str().c_str());
  for(auto i = config_.begin(); i != config_.end(); ++i){
    out_file << i->first << " " << i->second << std::endl;
  }
  out_file.close();
  return 0;
}

int HourglassConfiguration::LoadConfigFile(const std::string& in_file_name) {
  std::ifstream in_file(in_file_name.c_str());
  if(in_file) {
    std::string key;
    std::string val;
    while(in_file >> key >> val) {
      ModifyConfigParameter(key,val);
    }
  } else {
    std::cout << "Could not load config file: " << in_file_name << std::endl;
  }
  return 0;
}

int HourglassConfiguration::ShowConfigFile() {
  std::cout << "CONFIGURATION FILE: " << GetConfigName() << std::endl;
  for(auto i = config_.begin(); i != config_.end(); ++i) {
    std::cout << i->first << " = " << i->second << std::endl;
  }
  return 0;
}

std::string HourglassConfiguration::GetConfigName() {
  std::stringstream name;
  name                 << config_["RUN_NUMBER"]     << "_" 
      << "hoff"        << config_["X_OFFSET"]       << "_"
      << "voff"        << config_["Y_OFFSET"]       << "_"
      << "betastar"    << config_["BETA_STAR"]      << "_"
      << "crossing" << config_["CROSSING_ANGLE_XZ"] << "_" 
      << "hourglass_sim.conf";
  return name.str();
}
