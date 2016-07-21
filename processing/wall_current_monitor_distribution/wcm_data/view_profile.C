void view_profile(const std::string& profile_1,const std::string& profile_2){
  TGraph* g1 = new TGraph(profile_1.c_str());
  g1->SetName("profile_1");
  g1->SetMarkerStyle(kFullCircle);
  g1->SetMarkerStyle(kRed);
  g1->Draw("AP");
  TGraph* g2 = new TGraph(profile_2.c_str());
  g2->SetName("profile_2");
  g1->SetMarkerStyle(kFullCircle);
  g1->SetMarkerStyle(kGreen);
  g2->Draw("P");
}
