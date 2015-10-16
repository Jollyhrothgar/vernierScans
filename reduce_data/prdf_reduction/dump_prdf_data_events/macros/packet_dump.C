void packet_dump( 
    const std::string& lib = "/phenix/u/beaumim/VernierAnalysis/prdf_analysis/prdf_tools/full_prdf_scalers/scaler_dump/librun_scalers.so",
    const std::string& PRDF_name = "/gpfs02/phenix/spin3/beaumim/Run12VernierPRDF/files/VSCANDATA_P00-0000364636-0005.PRDFF",
    const std::string& out_file_name = "testout.txt"
    ) {
  gSystem->Load(lib.c_str());
  set_output_file(out_file_name);
  pfileopen(PRDF_name.c_str());
  prun(0);
  finish();
  std::cout << "Job done!" << std::endl;
}
