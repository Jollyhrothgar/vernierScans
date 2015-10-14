// Standard Template Library
#include <vector>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <map>
#include <iomanip>

// ROOT
#include "TStyle.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TPaveText.h"

// CUSTOM
#include "BBCRateSteps.h"
#include "BeamPositionSteps.h"
#include "BeamWidth.h"
#include "StyleFunctions.h"

BeamWidth::BeamWidth () { 
  std::cout << "Beam width calculation is underway." << std::endl;
}
BeamWidth::~BeamWidth () { 
  std::cout << "Beam width calculation is finished" << std::endl;
}

int BeamWidth::Init(
    const std::string& runNumber       ,
    const std::string& RootFileName    ,
    const std::string& BBCGL1PName     ,
    const std::string& ClockSumGL1PName,
    const std::string& ClockGL1PName   ,
    const std::string& StepsFileName   ,
    // const std::string& RateLossFileName,
    const std::string& BunchNumber     ,
    float livetimeRatio                ,
    const std::string& bpmData         ,
    const std::string& bpmSteps         
) {

  bpm.Init(
      runNumber,
      bpmData  ,
      bpmSteps
      );
  bpm.Run();
  
  bbc.Init(
      runNumber       ,
      RootFileName    ,
      BBCGL1PName     ,
      ClockSumGL1PName,
      ClockGL1PName   ,
      StepsFileName   ,
      //RateLossFileName,
      BunchNumber     ,
      livetimeRatio   
      );
  bbc.Run();

  std::string name;
  std::string title;
  NumberOfSteps = bbc.GetNumberOfSteps();

  HScan1 = new TGraphErrors();
  name = "HScan1_"+runNumber;
  title = "Horizontal Sweep - First Half of Scan, Run: " + runNumber + ";Horizontal Separation (microns);BBC Rate(Hz)";
  HScan1->SetName(name.c_str());
  HScan1->SetTitle(title.c_str());
  HScan1->SetMarkerStyle(kFullCircle);
  HScan1->SetMarkerColor(kRed);
  
  HScan2 = new TGraphErrors();
  name = "HScan2_"+runNumber;
  title = "Horizontal Sweep - Second Half of Scan, Run: " + runNumber + ";Horizontal Separation (microns);BBC Rate(Hz)";
  HScan2->SetName(name.c_str());
  HScan2->SetTitle(title.c_str());
  HScan2->SetMarkerStyle(kFullCircle);
  HScan2->SetMarkerColor(kRed+2);

  VScan1 = new TGraphErrors();
  name = "VScan1_"+runNumber;
  title = "Vertical Sweep - First Half of Scan, Run: " + runNumber + ";Vertical Separation (microns);BBC Rate(Hz)";
  VScan1->SetName(name.c_str());
  VScan1->SetTitle(title.c_str());
  VScan1->SetMarkerStyle(kFullCircle);
  VScan1->SetMarkerColor(kBlue);
  
  VScan2 = new TGraphErrors();
  name = "VScan2_"+runNumber;
  title = "Vertical Sweep - Second Half of Scan, Run: " + runNumber + ";Vertical Separation (microns);BBC Rate(Hz)";
  VScan2->SetName(name.c_str());
  VScan2->SetTitle(title.c_str());
  VScan2->SetMarkerStyle(kFullCircle);
  VScan2->SetMarkerColor(kBlue+2);

  BeamWidth2D = new TGraph2DErrors();
  name = "BeamWidth2D_"+runNumber;
  title = " Horizontal and Vertical Sweeps - Vernier Scan Run " + runNumber +";Horizontal Separation (microns);Vertical Separation (microns);BBC Rate(Hz)";
  BeamWidth2D->SetName(name.c_str());
  BeamWidth2D->SetTitle(title.c_str());
  BeamWidth2D->SetMarkerStyle(kFullCircle);
  BeamWidth2D->SetMarkerSize(1);
  BeamWidth2D->SetMarkerColor(21);

  return 0;
}
int BeamWidth::Run() {
  MakeBeamWidth();
  DrawBeamWidth();
  return 0;
}
int BeamWidth::MakeBeamWidth() {
  for(int i = 0; i < NumberOfSteps; i++) {
    float bbc_rate          = bbc.GetRateStep(i);  
    float bbc_rate_err      = bbc.GetRateStepErr(i);
    float bbc_ievent_center = bbc.GetRateStepEventIndexCenter(i);
    double bbc_ievent_min, bbc_ievent_max; 
    bbc.GetRateStepEventIndexRange(i, bbc_ievent_min, bbc_ievent_max);
    double bbc_min_rate, bbc_max_rate;
    bbc.GetRateStepRateRange(i, bbc_min_rate, bbc_max_rate);
    float bbc_step_rms      = bbc.GetRateStepRMS(i);
    float bpm_hsep          = bpm.GetHStep(i);
    float bpm_vsep          = bpm.GetVStep(i);
    float bpm_hsep_err      = bpm.GetHStepErr(i);
    float bpm_vsep_err      = bpm.GetVStepErr(i);

    rate.push_back(bbc_rate);  
    rate_step_stat_err.push_back(bbc_rate_err);
    rate_ievent_center.push_back(bbc_ievent_center);
    rate_max.push_back(bbc_max_rate); 
    rate_min.push_back(bbc_min_rate);
    rate_ievent_center.push_back(bbc_ievent_center);
    rate_ievent_min_range.push_back(bbc_ievent_min);
    rate_ievent_max_range.push_back(bbc_ievent_max);
    rate_step_rms.push_back(bbc_step_rms);
    hsep.push_back(bpm_hsep);
    vsep.push_back(bpm_vsep);
    hsep_err.push_back(bpm_hsep_err);
    vsep_err.push_back(bpm_vsep_err);

    BeamWidth2D->SetPoint     (BeamWidth2D->GetN()  , bpm_hsep    , bpm_vsep    , bbc_rate    );
    BeamWidth2D->SetPointError(BeamWidth2D->GetN()-1, bpm_hsep_err, bpm_vsep_err, bbc_rate_err);

    std::cout << " Step #" << i << "(x:" << bpm_hsep << ", y:" << bpm_vsep << ", " << bbc_rate << ") " << std::endl;
    
    if(i <= (NumberOfSteps/2 - 1) ) { // First Half of Steps
      std::cout << "1 : " << i << std::endl;
      HScan1 -> SetPoint     (HScan1->GetN()  , bpm_hsep, bbc_rate);
      VScan1 -> SetPoint     (VScan1->GetN()  , bpm_vsep, bbc_rate);
      HScan1 -> SetPointError(HScan1->GetN()-1, bpm_hsep_err, bbc_rate_err);
      VScan1 -> SetPointError(VScan1->GetN()-1, bpm_vsep_err, bbc_rate_err);
    } else { // Second Half of Steps
      std::cout << "2 : " << i << std::endl;
      HScan2 -> SetPoint     (HScan2->GetN()  , bpm_hsep, bbc_rate);
      VScan2 -> SetPoint     (VScan2->GetN()  , bpm_vsep, bbc_rate);
      HScan2 -> SetPointError(HScan2->GetN()-1, bpm_hsep_err, bbc_rate_err);
      VScan2 -> SetPointError(VScan2->GetN()-1, bpm_vsep_err, bbc_rate_err);
    }
  }
  float HScan1Integral = HScan1->Integral();
  float HScan2Integral = HScan2->Integral();
  
  if(HScan1Integral > HScan2Integral) {
    std::cout << "Horizontal Scan First" << std::endl;
    hScan = HScan1;
    vScan = VScan2;
    HorizontalScanFirst = true;
  } else {
    std::cout << "Vertical Scan First" << std::endl;
    hScan = HScan2;
    vScan = VScan1;
    HorizontalScanFirst = false;
  }

  // Determine Maximum and Minimum Ranges for Fit
  std::vector<double> hScanValues;
  std::vector<double> vScanValues;
  for(int i = 0; i < hScan->GetN(); i++) {
    double h = 0;
    double v = 0;
    double bbc = 0;
    hScan->GetPoint(i,h,bbc);
    vScan->GetPoint(i,v,bbc);
    hScanValues.push_back(h);
    vScanValues.push_back(v);
  }
  double min_range_h_sep = *std::min_element(hScanValues.begin(),hScanValues.end());
  double max_range_h_sep = *std::max_element(hScanValues.begin(),hScanValues.end());
  double min_range_v_sep = *std::min_element(vScanValues.begin(),vScanValues.end());
  double max_range_v_sep = *std::max_element(vScanValues.begin(),vScanValues.end());

  std::cout << "Fit range horizonal [" << min_range_h_sep << ", " << max_range_h_sep << "]" << std::endl;
  std::cout << "Fit range vertical  [" << min_range_v_sep << ", " << max_range_v_sep << "]" << std::endl;
  
  hScanFit1D = new TF1("hScanFit1D","gaus",min_range_h_sep,max_range_h_sep);
  vScanFit1D = new TF1("hScanFit1D","gaus",min_range_v_sep,max_range_v_sep);
  ScanFit2D  = new TF2("ScanFit2D" ,"xygaus",min_range_h_sep, max_range_h_sep, min_range_v_sep, max_range_v_sep);

  hScan->Fit(hScanFit1D);
  vScan->Fit(vScanFit1D);
  
  hScanMean     = hScanFit1D->GetParameter("Mean"     );
  hScanConstant = hScanFit1D->GetParameter("Constant" );
  hScanSigma    = hScanFit1D->GetParameter("Sigma"    );

  vScanMean     = vScanFit1D->GetParameter("Mean"     );
  vScanConstant = vScanFit1D->GetParameter("Constant" );
  vScanSigma    = vScanFit1D->GetParameter("Sigma"    );
                                                // ROOT Default Parameter Numbers for gausxy 
  ScanFit2D->SetParameter( 0 ,hScanConstant );  // 0  Constant
  ScanFit2D->SetParameter( 1 ,hScanMean     );  // 1  MeanX   
  ScanFit2D->SetParameter( 2 ,vScanSigma    );  // 2  SigmaX  
  ScanFit2D->SetParameter( 3 ,vScanMean     );  // 3  MeanY   
  ScanFit2D->SetParameter( 4 ,vScanSigma    );  // 4  SigmaY  
                                                                     // ROOT Default Parameter Numbers for gausxy 
  ScanFit2D->SetParLimits( 0 ,hScanConstant*0.6, hScanConstant*1.6); // 0  Constant
  ScanFit2D->SetParLimits( 1 ,-hScanMean*1.5, hScanMean*1.5);        // 1  MeanX   
  ScanFit2D->SetParLimits( 2 ,-hScanSigma*0.4, hScanSigma*1.6);      // 2  SigmaX  
  ScanFit2D->SetParLimits( 3 ,-vScanMean*1.5, vScanMean*1.5);        // 3  MeanY   
  ScanFit2D->SetParLimits( 4 ,vScanSigma*0.4, vScanSigma*1.6);       // 4  SigmaY  
  BeamWidth2D->Fit(ScanFit2D);

  ScanConstant = ScanFit2D->GetParameter("Constant");
  ScanSigmaX   = ScanFit2D->GetParameter("SigmaX"  );
  ScanSigmaY   = ScanFit2D->GetParameter("SigmaY"  );

  std::cout << "Nominal Beam Center (x,y): (" << hScanMean << ", " << vScanMean << ")" << std::endl;
  std::cout << "Max Overlap Rate Guess: " << hScanConstant << " Hz" << std::endl;
  std::cout << "Horizontal Beam Width Guess: " << hScanSigma << " microns" << std::endl;
  std::cout << "Vertical Beam Width Guess: " << vScanSigma << " microns" << std::endl;
  std::cout << "  Simultaneous Fit XSigma: " << ScanSigmaX << " microns" << std::endl;
  std::cout << "  Simultaneous Fit YSigma: " << ScanSigmaY << " microns" << std::endl;
  std::cout << "  Simultaneous Fit Max BBC Rate: " << ScanConstant << " Hz" << std::endl; 
  MakeCorrectedBeamWidth();
  return 0;
}

int BeamWidth::MakeCorrectedBeamWidth(){
  // Construct range of steps for each scan.
  std::cout << 0 << ", " << NumberOfSteps/2-1 << ", " << NumberOfSteps - 1 << std::endl;
  int scan_1_start = 0;
  int scan_1_end   = NumberOfSteps/2; // needs <
  //int scan_2_start = NumberOfSteps/2;
  //int scan_2_end   = NumberOfSteps;
  int j = NumberOfSteps/2;
  std::cout << std::setw(10) << "Step" << std::setw(10) << "hSep" << std::setw(10) << "vSep" << std::setw(10) << "Rate" << std::endl;
  for(int i = scan_1_start; i < scan_1_end; i++ ) { 
    //std::cout << i << ", " << j << std::endl;
    if(HorizontalScanFirst) {
      std::cout << std::setw(10) << i << ", " << j  << std::setw(10) << bpm.GetHStep(i) << std::setw(10) << bpm.GetVStep(j) << std::setw(10) << bbc.GetRateStep(i) << std::endl;
    } else {
    
    }
    j++;
  }
  /*
  if(HorizontalScanFirst) {
   
  } else {

  }
  */
  return 0;
}

int BeamWidth::DrawBeamWidth() {
  BeamWidthCanvas1D = new TCanvas("BeamWidthCanvas1D","Beam Width - Vernier Scan",1200,800);
  BeamWidthCanvas2D = new TCanvas("BeamWidthCanvas2D","Beam Width - Vernier Scan - Simultaneous Fit",1200,800);
  BeamWidthAll = new TCanvas("BeamWidthAll", " Beam Width - Vernier Scan", 1200, 800);

  setStyle();
  gStyle->SetOptFit();
  BeamWidthCanvas1D->Divide(2,2);
  BeamWidthCanvas1D->cd(1);
  HScan1->Draw("AP");
  BeamWidthCanvas1D->cd(2);
  HScan2->Draw("AP");
  BeamWidthCanvas1D->cd(3);
  VScan1->Draw("AP");
  BeamWidthCanvas1D->cd(4);
  VScan2->Draw("AP");

  BeamWidthCanvas2D->cd();
  BeamWidthCanvas2D->Divide(2,1);
  gStyle->SetOptFit();
  BeamWidthCanvas2D->cd(1);
  BeamWidth2D->Draw("pcolz");
  ScanFit2D->DrawCopy("surf same");
  BeamWidthCanvas2D->cd(2);
  ScanFit2D->DrawCopy("colz");

  BeamWidthAll->Divide(3,2);
  BeamWidthAll->cd(1);
  BeamWidth2D->Draw("pcolz");
  ScanFit2D->DrawCopy("surf same");
  BeamWidthAll->cd(2);
  hScan->Draw("AP");
  BeamWidthAll->cd(3);
  vScan->Draw("AP");
  BeamWidthAll->cd(4);
  gStyle->SetOptFit();
  gStyle->SetOptFit(1111);
  ScanFit2D->DrawCopy("colz");

  return 0;
}
