void read_trees(){
  // all trees in Sanghwa's star scalers area are called 'tree'
  TChain* ch = new TChain("tree");

  ch->Add("./ping/*.root");

  std::cout << "Got " << ch->GetEntries() << " trees." << std::endl;

  // Scaler ID:
  // 0: BBCLL1_0
  // 1: BBCLL1_1
  // 2: ZDCwide
  // 3: ZDCnarrow
  // 4: BBC south > 0 PMT
  // 5: BBC north > 0 PMT
  // 6-7: ERT Stuff, but not timed in for Run12
  // 8: ZDCS
  // 9: ZDCN
  // 10-15: MPC stuff, but not timed in for Run12
  // 16: BBC30 (BBCLL1_1 or both BBCLL1_0 and BBCLL_1)
  // 17: BBC15 (BBCLL1_1 and BBCLL1_0)
  // 18: BBCnoVtx (BBCLL1_1 or BBCLL_0)
  // 19: BUNCHCROSS
  // 20: BUNCHCROSS without busy flag
  // 21: BUNCHCROSS with busy flag
  // 22: BBCNandS (both south, north PMT > 0)
  // 23: BBCS only > 0 PMT (exclusive)
  // 24: BBCN only > 0 PMT (exclusive)
  // 25: BBCS only > 0 PMT && no BBCnoVtx flag (exclusive for xcheck)
  // 26: BBCN only > 0 PMT && no BBCnoVtx flag (exclusive for xcheck)
  // 27: BBC0 flag
  
  const int bbc_n_incl = 5; // BBC North inclusive live counts
  const int bbc_s_incl = 4; // BBC South inclusive live counts
  const int bbc_novertex_live = 18; // BBC(>0 tubes) Novertex live trigger
  const int clock_live = 20; // Clock scaler
  const int bbc_s_excl = 23; // BBC South exclusive trigger
  const int bbc_n_excl = 25; // BBC North exclusive trigger
  
  Int_t run_number;
  Long64_t scaler[28][120];

  ch->SetBranchAddress("Run_Num",&run_number);
  ch->SetBranchAddress("scaler",&scaler);

  for(int i = 0; i < ch->GetEntries(); i++){
    ch->GetEntry(i);
    std::cout << std::endl << std::endl << "RUN: " << run_number << std::endl;
    for(int j = 0; j < 120; j++){
      std::cout << scaler[bbc_novertex_live][j] << " ";
      if(j%10 == 0 && j > 0) std::cout << std::endl;
    }
  }
}
