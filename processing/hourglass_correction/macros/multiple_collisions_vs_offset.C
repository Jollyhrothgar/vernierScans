void multiple_collisions_vs_offset(){
  TGraph* g_mc_guess = new TGraph();
  g_mc_guess->SetName("g_mc_guess" );
  g_mc_guess->SetTitle("Multiple Collisions Guess;Beam Offset;Multiple Collision Rate Per Bunch Crossing");
  g_mc_guess->SetPoint(0, 0.0  , 0.435);
  g_mc_guess->SetPoint(1, 0.01 , 0.402);
  g_mc_guess->SetPoint(2, 0.025, 0.267);
  g_mc_guess->SetPoint(3, 0.040, 0.126);
  g_mc_guess->SetPoint(4, 0.060, 0.027);
  g_mc_guess->SetPoint(5, 0.090, 0.001);
  g_mc_guess->SetPoint(6, 0.1  , 0.001);
  std::cout << g_mc_guess->Eval(0.1);
  g_mc_guess->Draw("APL");
}
