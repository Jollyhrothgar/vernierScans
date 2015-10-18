#ifndef __CONFIGURATION__
#define __CONFIGTUATION__

#include<string>
#include <map>
class HourglassConfiguration{
 public:
  std::string this_name;
  HourglassConfiguration();
  ~HourglassConfiguration();
  int GenerateEmptyConfigFile();
  
  //  LoadConfigFile:
  // no checking is done to see if all config parameters are present, so
  // partial files may be loaded to enable mutliple sources to generate config
  // parameters.
  int LoadConfigFile(const std::string& in_file_name); 
  
  int ModifyConfigParameter(const std::string& param, const std::string& val);
  int SaveConfigFile(const std::string& out_dir);
  int ShowConfigFile();
  int SetDefaultValues();
  int BatchCreateConfigFile();
  std::string GetConfigName();

  // Setup Funcitons obtain config parameters from analysis environment
  int LoadWcmDcct(const std::string& wcm_dcct_file);
  int LoadBeamWidthParameters(const std::string& beam_width_parameters );
  int LoadHourGlassDataParameters(const std::string& hourglass_data_parameters );

  // returns true if par is found to be a parameter in config_. Designed usage
  // is for checking input when loading partial parameter lists from files.
  bool ParameterExists(const std::string& par);

 private:
  std::map<std::string,std::string> config_;
  /** config_: maps config parameter to config value. Conversion handled with
   * stringstreams. */
};

#endif
