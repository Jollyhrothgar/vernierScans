#include "VernierTreeVariables.h"

VernierTreeVariables::VernierTreeVariables() {
  Reset();
}

int VernierTreeVariables::Reset() {
  cross_id        = -999 ;
  evtnumber       = -999 ;
  trigraw         = -999 ;
  triglive        = -999 ;
  trigscaled      = -999 ;
  bbc_z           = -999.;
  zdcll1_z        = -999.;
  BbcNorthTubes   = -999 ;
  BbcSouthTubes   = -999 ;
  choice["READ"] = 0;
  choice["WRITE"] = 1;
  return 0;
}

int VernierTreeVariables::LinkTree(TTree* thisTree, const std::string& optStr) {

  if( choice.find(optStr) == choice.end() ) {
    std::cout << "VernierTreeVariables unrecognized option" << optStr << std::endl;
    std::cout << "Options are: " << std::endl; 
    for(std::map<std::string,int>::iterator it = choice.begin(); it != choice.end(); ++it) {
      std::cout << "  " << it->first << std::endl;
    }
  } else {
    int choice_int = choice[optStr];
    switch (choice_int) {
      case 0:
        thisTree->SetBranchAddress("cross_id"        , &(this->cross_id      ) ); 
        thisTree->SetBranchAddress("evtnumber"       , &(this->evtnumber     ) ); 
        thisTree->SetBranchAddress("trigraw"         , &(this->trigraw       ) ); 
        thisTree->SetBranchAddress("triglive"        , &(this->triglive      ) );
        thisTree->SetBranchAddress("trigscaled"      , &(this->trigscaled    ) );
        thisTree->SetBranchAddress("bbc_z"           , &(this->bbc_z         ) ); 
        thisTree->SetBranchAddress("zdcll1_z"        , &(this->zdcll1_z      ) ); 
        thisTree->SetBranchAddress("BbcNorthTubes"   , &(this->BbcNorthTubes ) ); 
        thisTree->SetBranchAddress("BbcSouthTubes"   , &(this->BbcSouthTubes ) ); 
        break;
      case 1:
        thisTree->Branch("cross_id"     , &(this->cross_id      ), "cross_id/I"     ); 
        thisTree->Branch("evtnumber"    , &(this->evtnumber     ), "evtnumber/I"    ); 
        thisTree->Branch("trigraw"      , &(this->trigraw       ), "trigraw/I"      ); 
        thisTree->Branch("triglive"     , &(this->triglive      ), "triglive/I"     );
        thisTree->Branch("trigscaled"   , &(this->trigscaled    ), "trigscaled/I"   );
        thisTree->Branch("bbc_z"        , &(this->bbc_z         ), "bbc_z/F"        ); 
        thisTree->Branch("zdcll1_z"     , &(this->zdcll1_z      ), "zdcll1_z/F"     ); 
        thisTree->Branch("BbcNorthTubes", &(this->BbcNorthTubes ), "BbcNorthTubes/I"); 
        thisTree->Branch("BbcSouthTubes", &(this->BbcSouthTubes ), "BbcSouthTubes/I"); 
        break;
    }
  }
  return 0;
}
