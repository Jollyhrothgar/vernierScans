void check_tree(){
  std::string root_file = "/direct/phenix+spin2/beaumim/vernierScans/data/run_12/dst_data/360879_reduced.root";
  TFile* f = new TFile(root_file.c_str(),"READ");

  std::cout << f->GetName() << std::endl;

  TTree* t = (TTree*)f->Get("VernierTree");
 
  t->Show(10);

  Int_t    cross_id     ;
  Int_t    evtnumber    ;
  Int_t    trigraw      ;
  Int_t    triglive     ;
  Int_t    trigscaled   ;
  Double_t bbc_z        ;
  Double_t zdc_z        ;
  Int_t    BbcNorthTubes;
  Int_t    BbcSouthTubes;

  int bbc_wide_trig_bit = 0x00000002;
  int zdc_wide_trig_bit = 0x00000004;

  t->SetBranchAddress("cross_id"     , &cross_id      );
  t->SetBranchAddress("evtnumber"    , &evtnumber     );
  t->SetBranchAddress("trigraw"      , &trigraw       );
  t->SetBranchAddress("triglive"     , &triglive      );
  t->SetBranchAddress("trigscaled"   , &trigscaled    );
  t->SetBranchAddress("bbc_z"        , &bbc_z         );
  t->SetBranchAddress("zdcll1_z"     , &zdc_z         );
  t->SetBranchAddress("BbcNorthTubes", &BbcNorthTubes );
  t->SetBranchAddress("BbcSouthTubes", &BbcSouthTubes );

  t->GetEntry(10);
  std::cout << cross_id      << std::endl;
  std::cout << evtnumber     << std::endl;
  std::cout << trigraw       << std::endl;
  std::cout << triglive      << std::endl;
  std::cout << trigscaled    << std::endl;
  std::cout << bbc_z         << std::endl;
  std::cout << zdc_z         << std::endl;
  std::cout << BbcNorthTubes << std::endl;
  std::cout << BbcSouthTubes << std::endl;

  TFile* out = new TFile("out.root","RECREATE");
  TH1F* bbc_zvtx[26];
  TH1F* zdc_zvtx[26];
  TCanvas* c[13];

  Int_t evt_range[26][2] = {
    {1468000,1506000}, 
    {1534000,1563000}, 
    {1610000,1694000}, 
    {1742000,1902000}, 
    {1955000,2156000}, 
    {2288000,2514000}, 
    {2627000,2882000}, 
    {3004000,3230000}, 
    {3353000,3551000}, 
    {3636000,3796000}, 
    {3862000,3947000}, 
    {3994000,4060000}, 
    {4088000,4116000}, 
    {4521000,4559000}, 
    {4597000,4634000}, 
    {4682000,4757000}, 
    {4805000,4945000}, 
    {5030000,5228000}, 
    {5346000,5567000}, 
    {5709000,5972000}, 
    {6080000,6305000}, 
    {6399000,6615000}, 
    {6686000,6831000}, 
    {6892000,6939000}, 
    {7028000,7061000}, 
    {7099000,7132000} 
  };
  
  char name[256];
  char title[256];
  for(int i = 0; i < 26; i++){
    if( i < 13 ){
      sprintf(name,"canvas_step_%02d_%02d",i,i+13);
      sprintf(title,"canvas_step_%02d_%02d",i,i+13);
      c[i] = new TCanvas(name,title,1200,800);
      
    }
    sprintf(name,"bbc_zvtx_%02d",i);
    sprintf(title,"BBC Z Vertex, Step %02d;z vtx;counts",i);
    bbc_zvtx[i] = new TH1F(name,title,100,-300,300);

    sprintf(name,"zdc_zvtx_%02d",i);
    sprintf(title,"ZDC Z Vertex, Step %02d;z vtx;counts",i);
    zdc_zvtx[i] = new TH1F(name,title,100,-300,300);

    bbc_zvtx[i]->Sumw2();
    zdc_zvtx[i]->Sumw2();
  }

  for(Long64_t entry = 0; entry < t->GetEntries(); entry++){
    t->GetEntry(entry);
    for(int step = 0; step < 26; step++){
      if((evtnumber > evt_range[step][0]) && (evtnumber < evt_range[step][1])) {
        if( ( trigscaled & bbc_wide_trig_bit ) > 0) {
          bbc_zvtx[step]->Fill(bbc_z); 
        } 
        if( ( trigscaled & zdc_wide_trig_bit ) > 0) {
          zdc_zvtx[step]->Fill(zdc_z); 
        }
      }
    }
  }
  out->cd();
  char save_file[256];
  for(int ci = 0; ci < 13; ci++){
    c[ci]->Divide(2,1);
    c[ci]->cd(1);
    zdc_zvtx[ci]->Draw();
    c[ci]->cd(2);
    zdc_zvtx[ci+13]->Draw();
    c[ci]->Draw();
    c[ci]->Write();
    bbc_zvtx[ci]->Write();
    zdc_zvtx[ci]->Write();
    bbc_zvtx[ci+13]->Write();
    zdc_zvtx[ci+13]->Write();
    sprintf(save_file,"step_compare_%02d_%02d.png",ci,ci+13);
    c[ci]->SaveAs(save_file);
  }
  out->Write();
  out->Close();
}
