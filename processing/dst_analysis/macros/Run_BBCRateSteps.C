#include "../../../FileManagement.h"
void Run_BBCRateSteps(
    int run_index                                  = 0,
    const std::string& bunchnumber                 = "",
    const std::string& bbcgl1p_histo_name          = "BBCGL1P",
    const std::string& clockgl1p_histo_name_sum    = "ClockGL1P",
    const std::string& clockgl1p_histo_name_bunch  = "ClockGL1P",
    float livetime_ratio                           = 1.0, 
    const std::string& lib                         = "libVernierAnalysis.so"
    ) 
{
  gSystem->Load(lib.c_str());
  gSystem->Load("/direct/phenix+u/workarea/beaumim/install/lib/Vernier_Dict_rdict.pcm");
  BBCRateSteps bbc;
  bbc.Init(
      run_number[run_index],
      reduced_dst_file[run_index],
      bbcgl1p_histo_name,
      clockgl1p_histo_name_sum,
      clockgl1p_histo_name_bunch,
      bbc_bin_step_boundaries[run_index],
      bunchnumber,
      livetime_ratio
      );
  bbc.Run();
}
