#ifndef __STYLE_FUNCTIONS__
#define __STYLE_FUNCTIONS__
#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TLine.h"

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetTextFont(132);
  gStyle->SetNdivisions(505,"X");
  gStyle->SetNdivisions(505,"Y");
  gStyle->SetTitleFont(132,"X");
  gStyle->SetLabelFont(132,"X");
  gStyle->SetTitleFont(132,"Y");
  gStyle->SetLabelFont(132,"Y");
  gStyle->SetTitleFont(132,"Z");
  gStyle->SetLabelFont(132,"Z");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleOffset(1.6,"X");
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetTitleOffset(1.8,"Y");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetTitleSize(0.05,"Z");
  gStyle->SetTitleOffset(1.8,"Z");
  gStyle->SetLabelSize(0.05,"Z");
  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadBottomMargin(0.2);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

void setStyleWide()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetTextFont(132);
  gStyle->SetNdivisions(505,"X");
  gStyle->SetNdivisions(505,"Y");
  gStyle->SetTitleFont(132,"X");
  gStyle->SetLabelFont(132,"X");
  gStyle->SetTitleFont(132,"Y");
  gStyle->SetLabelFont(132,"Y");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleOffset(1.2,"X");
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetTitleOffset(0.9,"Y");
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);
  
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTickLength(0.015,"X");
  gStyle->SetTickLength(0.015,"Y");
}

TH1F *drawFrame(Float_t xlow, Float_t ylow, Float_t xhigh, Float_t yhigh,
               TString xTitle="", TString yTitle="")
{
  TH1F *hFrame;

  if(gPad->GetLogy()==0) {
    hFrame = gPad->DrawFrame(xlow,ylow,xhigh,yhigh);
    hFrame->GetXaxis()->SetTitle(xTitle);
    hFrame->GetYaxis()->SetTitle(yTitle); }
  else {
    hFrame = gPad->DrawFrame(xlow,pow(10,ylow),xhigh,pow(10,yhigh));
    hFrame->GetXaxis()->SetTitle(xTitle);
    hFrame->GetYaxis()->SetTitle(yTitle); }

  if(ylow*yhigh<0.) {
    TLine *l = new TLine(xlow,0,xhigh,0);
    l->SetLineStyle(2);
    l->Draw(); }
  
  return hFrame;
}

#endif
