#include "../../../FileManagement.h"

void Run_BeamWidth( 
    int run_index = 0,
    const std::string& bbcgl1p_histo_name          = "BBCGL1P", // or BBC_GL1P_NN 
    const std::string& clockgl1p_histo_name_sum    = "ClockGL1P",
    const std::string& clockgl1p_histo_name_bunch  = "ClockGL1P",
    const std::string& bunchnumber                 = "", // or NN
    int livetimeRatio                              = 1.0,
    const std::string& vernierAnalysisLib          = "libVernierAnalysis.so"
    )
{
  gSystem->Load(vernierAnalysisLib.c_str());
  BeamWidth bw;
  bw.Init(
      run_number[run_index],
      reduced_dst_file[run_index],
      bbcgl1p_histo_name,
      clockgl1p_histo_name_sum,
      clockgl1p_histo_name_bunch,
      bbc_bin_step_boundaries[run_index],
      bunchnumber,
      livetimeRatio,
      bpm_file[run_index],
      relative_time_step_boundaries[run_index]
      );
  bw.Run();
}
