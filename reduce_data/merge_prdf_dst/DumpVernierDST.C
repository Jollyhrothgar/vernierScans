#include "../../FileManagement.h"
void DumpVernierDST(
  int run_index = 0,
  const std::string& output_directory = ".",
  const std::string& vernier_dst_lib = "libVernierAnalysis.so",
  const std::string& vernier_tree_name = "VernierTree",
		){
  gSystem->Load(vernier_dst_lib.c_str());
  VernierTreeVariables v;

  TFile* f = new TFile(reduced_dst_file[run_index].c_str(),"READ");
  TTree* t = (TTree*)f->Get(vernier_tree_name.c_str());
  std::string out_file_name = run_number[run_index] + "_DST_dump.txt";
  std::ofstream out_file(out_file_name.c_str());
  v.LinkTree(t,"READ");

  Long64_t entries = t->GetEntries();
  for(Long64_t i = 0; i < entries; i++){
    t->GetEntry(i);
    out_file << v.evtnumber
        << " " << v.trigraw
        << " " << v.triglive
        << " " << v.trigscaled
        << " " << v.bbc_z 
        << " " << v.zdcll1_z
        << std::endl;
    /*
    if(i%10000 == 0) {
      std::cout << "Progress: " << i << "/" << entries << "                \r" << std::flush;
    }
    */
  }
  std::cout << std::endl << "Done" << std::endl;
}
