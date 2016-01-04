#include<iostream>
#include<fstream>
void MakeHistogram(const std::string& file) {
  std::ifstream in_file(file.c_str());
  TH1F* h = new TH1F("h","title",100,-300,300);
  float zvtx = 0;
  while(in_file >> zvtx){
    h->Fill(zvtx);
  }
  h->Draw();
}
