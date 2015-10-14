#include "../../../FileManagement.h"
int Run_BBCEfficiency(
    int run_index   = 0, 
    std::string lib = "libVernierAnalysis.so"
    ) {
  
  gSystem->Load("/direct/phenix+u/workarea/beaumim/install/lib/libVernierAnalysis.so");
  BBCEfficiency bbc_eff;
  bbc_eff.Init(
      run_number[run_index],
      reduced_dst_file[run_index],
      "VernierTree"
      );
  bbc_eff.Run();
  bbc_eff.FitAcceptance(-42.,-28.,25.,44.);
  bbc_eff.GetCorrectedEfficiency();
  bbc_eff.MakeFigures("/direct/phenix+spin2/beaumim/vernierScans/plots");
  return 0;
}
