#include "../../../FileManagement.h"
#include <sstream>
int Run_BBCRate(
    int run_index = 0,
    int bunch_index = -1,
    const std::string& lib = "libVernierTimeAnalysis.so"
    ) {
  gSystem->Load(lib.c_str());

  std::stringstream run;
  std::stringstream file;
  if( bunch_index >= 0 && bunch_index < 120) {
    run.str("");
    file.str("");
    run << run_number[run_index] << "_bunch_" << bunch_index;
    file << scaler_file_bunch_stub[run_index] << bunch_index << ".txt";
  } else {
    run.str("");
    file.str("");
    run << run_number[run_index];
    file << scaler_file[run_index]; 
  }
  BBCRate bbc_rate;
  bbc_rate.Init(
      file.str(),
      run.str()
      );
  bbc_rate.Run();
  bbc_rate.MakeFigures("/direct/phenix+spin2/beaumim/vernierScans/plots");
  return 0;
}
