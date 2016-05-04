#include "../../../FileManagement.h"
#include <map>

int Run_HourglassSimulation(
  const std::string& config_file,
  const std::string& plot_output = plots_dir,
  const std::string& file_stub,
  const int simulation_mode = "",
  const int zprofile_mode = "",
  const std::string lib = "libVernierHourglass.so"
) {
  
  std::map<std::string,int> my_map;

  int run_index = 0;
  gSystem->Load(lib.c_str());
  HourglassSimulation sim;
  sim.SetSaveFileStub(file_stub);
  sim.Init(
      config_file,
      hourglass_data_file[run_index],
      z_profile_density_blue[run_index],
      z_profile_density_yellow[run_index],
      z_profile_fit_file[run_index]
      );
  sim.SetSaveDirectory(plot_output);
  sim.SetSaveAll();
  switch (simulation_mode) {
    case 0:
      sim.Run(zprofile_mode);
      sim.Compare();
      break;
    case 1:
      sim.RunRootFinder(zprofile_mode,hourglass_data_file[run_index]);
      break;
    default:
      std::cout << "no adequate running mode supplied." << std::endl;
      std::cout << " choose 1 for RootFinder mode" << std::endl;
      std::cout << " choose 0 for Normal Running mode " << std::endl;
  } 
  sim.SaveFigures();
  return 0;
}
