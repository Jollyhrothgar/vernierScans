void MakeBpmDisplacements() {
  TFile* f = new TFile("./360879_BeamPositionPlots.root","READ");
  TH1F* x[26];
  TH1F* y[26];
  char name[256];
  char title[256];
  double x_m[26];
  double y_m[26];

  for(int i = 0; i < 26; i++){
    sprintf(name,"bpm_x_step_%d",i);
    x[i] = (TH1F*)f->Get(name);
    sprintf(name,"bpm_y_step_%d",i);
    y[i] = (TH1F*)f->Get(name);

    x_m[i]=(double)x[i]->GetMean();
    y_m[i]=(double)y[i]->GetMean();
    x_m[i]/=10000.;
    y_m[i]/=10000.;

    double disp = sqrt(x_m[i]*x_m[i]+y_m[i]*y_m[i]);
    double fac = 1.0;
    if((i < 13) && (x_m[i] < 0)){
      fac = -1.0;
    } else if ((i >= 13)&&(y_m[i] < 0)){
      fac = -1.0;
    }

    disp = disp*fac;

    printf("%f\n",disp);
  }
}
