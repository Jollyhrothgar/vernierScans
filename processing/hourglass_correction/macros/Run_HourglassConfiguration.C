#include "../../../FileManagement.h"

int Run_HourglassConfiguration(
) {
  gSystem->Load("libVernierHourglass.so");
  HourglassConfiguration config;

  int run_index = 0;
  // MODE 1:  Create All Config Files
  // Batch Configure
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
 
  // wp need to test range based config generation
  // Shows the last internal configuration of HourglassConfiguration
  config.ShowConfigFile();

  // MODE 2:  Config Range
  // config.GenerateConfigFileRange("359711_step_00.conf","./","359711_step0_var_",1.0,10);

  // MODE 3: Single File
  // Running one config file -- example
  // con.ModifyConfigParameter("RUN_NUMBER","359711");
  // con.SetDefaultValues();
  // con.ShowConfigFile();
  // con.SaveConfigFile("/direct/phenix+spin2/beaumim/vernierScans/data/run_12/simulation_config");
  return 0;
}
