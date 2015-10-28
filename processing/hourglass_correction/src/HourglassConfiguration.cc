#include "HourglassConfiguration.h"

#include <iostream>
#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <iomanip>

HourglassConfiguration::HourglassConfiguration() {
  std::stringstream ss;
  ss << "HourglassConfiguration_" << std::hex << this;
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
  par_["RUN_NUMBER"]                     = ""; // FileManagement  
  par_["ZDC_COUNTS"]                     = ""; // HourglassData * (done)
  par_["X_OFFSET"]                       = ""; // BeamWidthTime * (done)
  par_["Y_OFFSET"]                       = ""; // BeamWidthTime * (done)
  par_["HORIZONTAL_BEAM_WIDTH"]          = ""; // BeamWidthTime (done)
  par_["VERTICAL_BEAM_WIDTH"]            = ""; // BeamWidthTime (done)
  par_["AVG_NUMBER_IONS_BLUE_BEAM"]      = ""; // WcmDcct (todo) 
  par_["AVG_NUMBER_IONS_YELLOW_BEAM"]    = ""; // WcmDcct (todo)
  par_["BBC_ZDC_Z_VERTEX_OFFSET"]        = ""; // HourglassData (done)
  par_["ZDC_VERTEX_DISTRIBUTION_NAME"]   = ""; // HourglassData (done)
  par_["BETA_STAR" ]                     = ""; // GOAL: find this, SetDefaultValues
  par_["CROSSING_ANGLE_XZ"]              = ""; // GOAL: find this, SetDefaultValues
  par_["FILLED_BUNCHES"]                 = ""; // Fixed, SetDefaultValues
  par_["BUNCH_CROSSING_FREQUENCY"]       = ""; // Fixed, SetDefaultValues
  par_["Z_PROFILE_SCALE_VALUE"]          = ""; // Fixed, SetDefaultValues
  par_["MULTIPLE_COLLISION_RATE"]        = ""; // wp SetDefaultValues *
  par_["MAX_COLLISIONS"]                 = ""; // wp SetDefaultValues
  par_["Z_BUNCH_WIDTH_LEFT_GAUSSIAN"]    = ""; // wp SetDefaultValues
  par_["Z_BUNCH_WIDTH_RIGHT_GAUSIAN"]    = ""; // wp SetDefaultValues
  par_["Z_BUNCH_WIDTH_CENTRAL_GAUSIAN"]  = ""; // wp SetDefaultValues
  par_["Z_BUNCH_WIDTH_LEFT_OFFSET"]      = ""; // wp SetDefaultValues
  par_["Z_BUNCH_WIDTH_RIGHT_OFFSET"]     = ""; // wp SetDefaultValues
  // *these terms change based on the beam displaement
  return 0;
}

int HourglassConfiguration::SetDefaultValues() {
  // Default values are set to the first scan step of Run 359711, a vernier scan
  // from Run 12.
  ModifyConfigParameter("RUN_NUMBER"                    , "359711");
  ModifyConfigParameter("ZDC_COUNTS"                    , "891");
  ModifyConfigParameter("X_OFFSET"                      , "-0.1");
  ModifyConfigParameter("Y_OFFSET"                      , "0.0");
  ModifyConfigParameter("HORIZONTAL_BEAM_WIDTH"         , "0.0245674");
  ModifyConfigParameter("VERTICAL_BEAM_WIDTH"           , "0.0238342");
  ModifyConfigParameter("AVG_NUMBER_IONS_BLUE_BEAM"     , "120.029e9");
  ModifyConfigParameter("AVG_NUMBER_IONS_YELLOW_BEAM"   , "88.167e9");
  ModifyConfigParameter("BBC_ZDC_Z_VERTEX_OFFSET"       , "-9.38");
  ModifyConfigParameter("ZDC_VERTEX_DISTRIBUTION_NAME"  , "zdc_zvtx_step_0");
  ModifyConfigParameter("BETA_STAR"                     , "85.");
  ModifyConfigParameter("CROSSING_ANGLE_XZ"             , "-0.08e-3");
  ModifyConfigParameter("FILLED_BUNCHES"                , "107");
  ModifyConfigParameter("BUNCH_CROSSING_FREQUENCY"      , "78213.");
  ModifyConfigParameter("Z_PROFILE_SCALE_VALUE"         , "1.5");
  ModifyConfigParameter("MULTIPLE_COLLISION_RATE"       , "0.001");
  ModifyConfigParameter("MAX_COLLISIONS"                , "5");
  ModifyConfigParameter("Z_BUNCH_WIDTH_LEFT_GAUSSIAN"   , "35.15"); 
  ModifyConfigParameter("Z_BUNCH_WIDTH_RIGHT_GAUSIAN"   , "27.65"); 
  ModifyConfigParameter("Z_BUNCH_WIDTH_CENTRAL_GAUSIAN" , "55.95"); 
  ModifyConfigParameter("Z_BUNCH_WIDTH_LEFT_OFFSET"     , "-70.2"); 
  ModifyConfigParameter("Z_BUNCH_WIDTH_RIGHT_OFFSET"    , "56.7"); 
  return 0;
}

std::string HourglassConfiguration::GetPar( const std::string&  par_name ) {
  if(ParameterExists(par_name)) {
    return par_[par_name];
  }
  return "-9999";
}

int HourglassConfiguration::ModifyConfigParameter(const std::string& name, const std::string& val) {
  if(ParameterExists(name)){
    par_[name] = val;
  } else {
    std::cout << "Unrecognized configuration parameter, remember parameters are case-senstive: " << name << std::endl;
  }
  return 0;
}

bool HourglassConfiguration::ParameterExists(const std::string& par) {
  auto search = par_.find(par);
  if(search != par_.end() ) return true;
  return false;
}

int HourglassConfiguration::SaveConfigFile(const std::string& file_name ) {
  std::ofstream out_file(file_name.c_str());
  for(auto i = par_.begin(); i != par_.end(); ++i){
    out_file << i->first << " " << i->second << std::endl;
  }
  out_file.close();
  return 0;
}

int HourglassConfiguration::LoadConfigFile(const std::string& in_file_name) {
  std::ifstream in_file(in_file_name.c_str());
  if(in_file) {
    std::string line;
    while(getline(in_file,line)) {
      //debug std::cout << line << std::endl;
      std::vector<std::string> tokens;
      std::string token;
      std::istringstream iss(line.c_str());
      while(getline(iss,token,' ')) { tokens.push_back(token); }
      if(tokens.size() != 2 ) {
        std::cout << "FILE " << in_file_name << " IS NOT FULLY INITIALIZED. CANNOT USE FOR SIMULATIONS!" << std::endl;
	return 1;
      } else {
        ModifyConfigParameter(tokens[0],tokens[1]);
      }
    }
  } else {
    std::cout << "Could not load config file: " << in_file_name << std::endl;
  }
  return 0;
}

int HourglassConfiguration::ShowConfigFile() {
  std::cout << "CONFIGURATION FILE: " << std::endl;
  for(auto i = par_.begin(); i != par_.end(); ++i) {
    std::cout << i->first << " = " << i->second << std::endl;
  }
  return 0;
}

// Default name for interactively generated config files.
std::string HourglassConfiguration::GetConfigName() {
  std::stringstream name;
  name              << par_["RUN_NUMBER"]     << "_" 
      << "hoff"     << par_["X_OFFSET"]       << "_"
      << "voff"     << par_["Y_OFFSET"]       << "_"
      << "sim.conf";
  return name.str();
}

int HourglassConfiguration::BatchCreateConfigFiles(
  const std::string& run_number_,
  const std::string& zdc_bbc_offset_file,
  const std::string& zdc_counts_per_step_file,
  const std::string& x_offsets_file,
  const std::string& y_offsets_file,
  const std::string& h_width_file,
  const std::string& v_width_file,
  const std::string& beam_population_file, 
  const std::string& zdc_vertex_histo_name_file,
  const std::string& sim_config_out_dir
  ) {
  SetDefaultValues();
  par_["RUN_NUMBER"] = run_number_;

  std::vector< std::map < std::string, std::string > > xoff;
  std::vector< std::map < std::string, std::string > > yoff;
  std::vector< std::map < std::string, std::string > > zdc_count;
  std::vector< std::map < std::string, std::string > > zdc_hist_name;

  std::ifstream in_zdc_off(zdc_bbc_offset_file.c_str());
  std::ifstream in_zdc(zdc_counts_per_step_file.c_str());
  std::ifstream in_xoff(x_offsets_file.c_str());
  std::ifstream in_yoff(y_offsets_file.c_str());
  std::ifstream in_hwidth(h_width_file.c_str());
  std::ifstream in_vwidth(v_width_file.c_str());
  std::ifstream in_beam(beam_population_file.c_str());
  std::ifstream in_zdc_histo_name(zdc_vertex_histo_name_file.c_str());

  // Check for bad file handles
  if(!in_zdc_off) {
    std::cerr << "couldn't open " << zdc_bbc_offset_file << std::endl;
    return 1;
  }
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
  if(!in_zdc_histo_name) {
    std::cerr << "couldd't open " << zdc_vertex_histo_name_file << std::endl;
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
 
  // Load BBC ZDC offset
  while(in_zdc_off >> key >> val) {
    ModifyConfigParameter(key,val);
  }
  in_zdc_off.close();

  // Load ZDC histogram names for comparison to simulation
  while(in_zdc_histo_name >> key >> val) {
    std::map< std::string, std::string > m;
    m[key] = val;
    zdc_hist_name.push_back(m);
  }

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
  
  // Load horizontal width
  while(in_hwidth >> key >> val ) {
    ModifyConfigParameter(key,val);
  }
  in_hwidth.close();

  // Load vertical width
  while(in_vwidth >> key >> val ) {
    ModifyConfigParameter(key,val);
  }
  in_vwidth.close();
  
  // load beam populations for blue and yellow
  while(in_beam >> key >> val) {
    ModifyConfigParameter(key,val);
  }
  in_beam.close();

  unsigned const int number_of_steps = zdc_count.size();
  if( (xoff.size() != number_of_steps) && (yoff.size() != number_of_steps) ) {
    std::cerr << "You have inconsistant numbers of steps in your config parameter files:" << std::endl;
    std::cerr << x_offsets_file << " : " << xoff.size() << std::endl;
    std::cerr << y_offsets_file << " : " << yoff.size() << std::endl;
    std::cerr << zdc_counts_per_step_file << " : " << zdc_count.size() << std::endl;
    return 1;
  } 
  std::cout << "Generating config files for run: " << par_["RUN_NUMBER"];
  std::cout << "Sending config files to: " << sim_config_out_dir << std::endl;
  for(unsigned int i = 0; i < number_of_steps; i++) {
    auto xoff_step = xoff[i]; // get the map
    auto yoff_step = yoff[i];
    auto zdc_step  = zdc_count[i]; 
    auto zdc_hist  = zdc_hist_name[i];

    ModifyConfigParameter(xoff_step.begin()->first , xoff_step.begin()->second);
    ModifyConfigParameter(yoff_step.begin()->first , yoff_step.begin()->second);
    ModifyConfigParameter(zdc_step .begin()->first , zdc_step .begin()->second );
    ModifyConfigParameter(zdc_hist .begin()->first , zdc_hist .begin()->second );

    std::stringstream file_name;
    file_name << sim_config_out_dir << "/" << run_number_ << "_step_" << std::setw(2) << std::setfill('0') << i << ".conf";
    SaveConfigFile(file_name.str());
  }
  return 0;
}
