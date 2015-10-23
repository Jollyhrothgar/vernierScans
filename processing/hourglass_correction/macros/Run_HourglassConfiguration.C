#include "../../../FileManagement.h"

int Run_HourglassConfiguration(
) {
  gSystem->Load("libVernierHourglass.so");
  HourglassConfiguration config;

  int run_index = 0;
  // Create All Config Files
  config.BatchCreateConfigFiles(
    run_number[run_index],
    zdc_bbc_offset_sim_config[run_index],
    zdc_counts_per_step[run_index],
    x_offset_sim_config[run_index],
    y_offset_sim_config[run_index],
    h_width_sim_config[run_index],
    v_width_sim_config[run_index],
    beam_population_sim_config[run_index],
    zdc_zvtx_histo_name_sim_config[run_index],
    simulation_config_dir
    );
  // Shows the last internal configuration of HourglassConfiguration
  config.ShowConfigFile();

  // Running one config file -- example
  // con.ModifyConfigParameter("RUN_NUMBER","359711");
  // con.SetDefaultValues();
  // con.ShowConfigFile();
  // con.SaveConfigFile("/direct/phenix+spin2/beaumim/vernierScans/data/run_12/simulation_config");
  return 0;
}
