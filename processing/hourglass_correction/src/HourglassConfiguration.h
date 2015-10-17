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
  int LoadConfigFile(const std::string& in_file_name); 
  /** LoadConfigFile:
   * no checking is done to see if all config parameters are present, so
   * partial files may be loaded to enable mutliple sources to generate config
   * parameters. */
  int ModifyConfigParameter(const std::string& param, const std::string& val);
  int SaveConfigFile(const std::string& out_dir);
  int ShowConfigFile();
  int SetDefaultValues();
  std::string GetConfigName();

  /** Setup Funcitons 
   * Obtain config parameters from analysis environment */
  int SetupWcmDcct(); /** wp, load a text file with the parameters */
  int SetupBeamWidthTime(); /** wp, load a text file with the parameters */
  int SetupBeamPositionSteps(); /** wp, load a text file with the parameters */
  int SetupHourglassData(); /** wp, load a text file with the parameters */ 

 private:
  std::map<std::string,std::string> config_;
  /** config_: maps config parameter to config value. Conversion handled with
   * stringstreams. */
};

#endif
