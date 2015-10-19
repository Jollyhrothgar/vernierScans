#include "HourglassConfiguration.h"

#include <iostream>
#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

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
  config_["RUN_NUMBER"]                     = ""; // FileManagement 
  config_["ZDC_COUNTS"]                     = ""; // HourglassData * (done)
  config_["X_OFFSET"]                       = ""; // BeamWidthTime * (done)
  config_["Y_OFFSET"]                       = ""; // BeamWidthTime * (done)
  config_["HORIZONTAL_BEAM_WIDTH"]          = ""; // BeamWidthTime (done)
  config_["VERTICAL_BEAM_WIDTH"]            = ""; // BeamWidthTime (done)
  config_["AVG_NUMBER_IONS_BLUE_BEAM"]      = ""; // WcmDcct (todo) 
  config_["AVG_NUMBER_IONS_YELLOW_BEAM"]    = ""; // WcmDcct (todo)
  config_["BBC_ZDC_Z_VERTEX_OFFSET"]        = ""; // HourglassData (done)
  config_["BETA_STAR" ]                     = ""; // GOAL: find this, SetDefaultValues
  config_["CROSSING_ANGLE_XZ"]              = ""; // GOAL: find this, SetDefaultValues
  config_["FILLED_BUNCHES"]                 = ""; // Fixed, SetDefaultValues
  config_["BUNCH_CROSSING_FREQUENCY"]       = ""; // Fixed, SetDefaultValues
  config_["Z_PROFILE_SCALE_VALUE"]          = ""; // Fixed, SetDefaultValues
  config_["MULTIPLE_COLLISION_RATE"]        = ""; // wp SetDefaultValues *
  config_["MAX_COLLISIONS"]                 = ""; // wp SetDefaultValues
  config_["Z_BUNCH_WIDTH_LEFT_GAUSSIAN"]    = ""; // wp SetDefaultValues
  config_["Z_BUNCH_WIDTH_RIGHT_GAUSIAN"]    = ""; // wp SetDefaultValues
  config_["Z_BUNCH_WIDTH_CENTRAL_GAUSIAN"]  = ""; // wp SetDefaultValues
  config_["Z_BUNCH_WIDTH_LEFT_OFFSET"]      = ""; // wp SetDefaultValues
  config_["Z_BUNCH_WIDTH_RIGHT_OFFSET"]     = ""; // wp SetDefaultValues
  // *these terms change based on the beam displaement
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
  if(ParameterExists(name)){
    config_[name] = val;
  } else {
    std::cout << "Unrecognized configuration parameter, remember parameters are case-senstive: " << name << std::endl;
  }
  return 0;
}

bool HourglassConfiguration::ParameterExists(const std::string& par) {
  auto search = config_.find(par);
  if(search != config_.end() ) return true;
  return false;
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
  name              << config_["RUN_NUMBER"]     << "_" 
      << "hoff"     << config_["X_OFFSET"]       << "_"
      << "voff"     << config_["Y_OFFSET"]       << "_"
      << "betastar" << config_["BETA_STAR"]      << "_"
      << "crossing" << config_["CROSSING_ANGLE_XZ"] << "_" 
      << "hourglass_sim.conf";
  return name.str();
}

int HourglassConfiguration::BatchCreateConfigFiles(
  const std::string& run_number_,
  const std::string& zdc_counts_per_step_file,
  const std::string& x_offsets_file,
  const std::string& y_offsets_file,
  const std::string& h_width_file,
  const std::string& v_width_file,
  const std::string& beam_population_file, 
  const std::string& sim_config_out_dir
  ) {
  SetDefaultValues();
  config_["RUN_NUMBER"] = run_number_;

  std::vector< std::map < std::string, std::string > > xoff;
  std::vector< std::map < std::string, std::string > > yoff;
  std::vector< std::map < std::string, std::string > > zdc_count;

  std::ifstream in_zdc(zdc_counts_per_step_file.c_str());
  std::ifstream in_xoff(x_offsets_file.c_str());
  std::ifstream in_yoff(y_offsets_file.c_str());
  std::ifstream in_hwidth(h_width_file.c_str());
  std::ifstream in_vwidth(v_width_file.c_str());
  std::ifstream in_beam(beam_population_file.c_str());

  // Check for bad file handles
  if(!in_zdc) {
    std::cerr << "couldn't open " << zdc_counts_per_step_file << std::endl;
    return 1;
  } 
  if(!in_xoff) {
    std::cerr << "couldn't open " << x_offsets_file << std::endl;
    return 1;
  } 
  if(!in_yoff) {
    std::cerr << "couldn't open " << y_offsets_file << std::endl;
    return 1;
  } 
  if(!in_hwidth) { 
    std::cerr << "couldn't open " << h_width_file << std::endl;
    return 1;
  } 
  if(!in_vwidth) {
    std::cerr << "couldn't open " << v_width_file << std::endl;
    return 1;
  } 
  if(!in_beam) {
    std::cerr << "couldn't open " << beam_population_file << std::endl;
    return 1;
  } 
  // Load ZDC steps
  std::string key;
  std::string val;
  while(in_zdc >> key >> val) {
    std::map<std::string,std::string> m;
    m[key] = val;
    zdc_count.push_back(m);
  }
  in_zdc.close();

  // Load xoff steps
  while(in_xoff >> key >> val) {
    std::map<std::string,std::string> m;
    m[key] = val;
    xoff.push_back(m);
  }
  in_xoff.close();

  // Load yoff steps
  while(in_yoff >> key >> val ) {
    std::map<std::string,std::string> m;
    m[key] = val;
    yoff.push_back(m);
  }
  in_yoff.close();
  unsigned const int number_of_steps = zdc_count.size();
  if( (xoff.size() != number_of_steps) && (yoff.size() != number_of_steps) ) {
    std::cerr << "You have inconsistant numbers of steps in your config parameter files:" << std::endl;
    std::cerr << x_offsets_file << " : " << xoff.size() << std::endl;
    std::cerr << y_offsets_file << " : " << yoff.size() << std::endl;
    std::cerr << zdc_counts_per_step_file << " : " << zdc_count.size() << std::endl;
    return 1;
  } 
  std::cout << "Generating config files for run: " << config_["RUN_NUMBER"];
  std::cout << "Sending config files to: " << sim_config_out_dir << std::endl;
  
  while(in_hwidth >> key >> val ) {
    ModifyConfigParameter(key,val);
  }
  in_hwidth.close();
  while(in_vwidth >> key >> val ) {
    ModifyConfigParameter(key,val);
  }
  in_vwidth.close();
  while(in_beam >> key >> val) {
    ModifyConfigParameter(key,val);
  }
  in_beam.close();

  for(unsigned int i = 0; i < number_of_steps; i++) {
    auto xoff_step = xoff[i]; // get the map
    auto yoff_step = yoff[i];
    auto zdc_step  = zdc_count[i]; 

    ModifyConfigParameter(xoff_step.begin()->first,xoff_step.begin()->second);
    ModifyConfigParameter(yoff_step.begin()->first,yoff_step.begin()->second);
    ModifyConfigParameter(zdc_step.begin()->first ,zdc_step.begin()->second );
    auto config_name = GetConfigName();
  }

  return 0;
}
