#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLine.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TEllipse.h"
#include "TStyle.h"
#include "TList.h"
#include "TObject.h"
#include "TColor.h"
#include "TPaveStats.h"
#include "TKey.h"
#include "TROOT.h"
#include "TVirtualFFT.h"
#include <cmath>
#include <iostream>
#include <vector>
#include "boost/format.hpp"
#include <sstream>
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TText.h"
#include "TClass.h"
#include "TLegendEntry.h"
#include "TGaxis.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include <cstdlib>

using std::string;
using std::stringstream;
using std::vector;

void setStyle();

// taken from FFT of number of tracks vs time plot
double fCBO = 415; // kHz

int main(int argc, char* argv[]) {
  
  setStyle();

  int xmin =0.0;
  int xmax = 0.0;
  int nbins = 0.0;

  if (argc != 3) {
    xmin = 80;
    xmax = 400;
  }
  else {
    xmin = atoi(argv[1]);
    xmax = atoi(argv[2]);
    std::cout << "xmin: " << xmin << ", xmax: " << xmax << "\n";
  }

  stringstream ss;
  ss << xmin << "_" << xmax << ".eps";
  string app = ss.str();
  
  //TFile* f  = new TFile("/gm2/data/g2be/Production/Plots/7861/trackRecoPlots.root");
  //TFile* f  = new TFile("/gm2/data/g2be/Production/Plots/SummaryPlots/trackRecoPlots_LostMuonRuns.root");
  //TFile* f  = new TFile("/gm2/data/g2be/Production/Plots/SummaryPlots/trackRecoPlots_GoodRuns.root");
  TFile* f = new TFile("/gm2/data/g2be/RFQuadRuns/Plots/8858/trackRecoPlots.root");

  TDirectoryFile* extrap = (TDirectoryFile*)f->Get("Extrapolation");
  TDirectoryFile* vertex = (TDirectoryFile*)extrap->Get("vertices");
  TDirectoryFile* allStations = (TDirectoryFile*)vertex->Get("allStations");
  TDirectoryFile* allEvents = (TDirectoryFile*)allStations->Get("allEvents");
  TDirectoryFile* passEvents = (TDirectoryFile*)allStations->Get("pValue>0.000_and_noVolumesHit");

  TCanvas* c1 = new TCanvas("c1","c1",800,600);
  c1->SetLeftMargin(0.13);
  c1->SetBottomMargin(0.13);

  TH2F* radPosVsTime = (TH2F*)allEvents->Get("h_radialPos_vs_time");
  
  double cbo_period = 1.0 / (fCBO * 1000);
  cbo_period *= 1e6;
  double binw = radPosVsTime->GetXaxis()->GetBinWidth(1);
  std::cout << "cbo_period = " << cbo_period << ", binw = " << binw << ", cbo_period/binw = " << cbo_period / binw << "\n";

  double x1 = -44.9;
  double x2 = 44.9;
  TH1F* radialPosHistAll = (TH1F*)radPosVsTime->ProjectionY();
  
  radialPosHistAll->Draw();
  c1->SaveAs("proj.eps");

  // plot of the mean radial pos vs time
  int nbin_low  = radialPosHistAll->FindBin(x1);
  int nbin_high = radialPosHistAll->FindBin(x2);
  //TProfile* avgRadPosVsTime = (TProfile*)radPosVsTime->ProfileX("prof1",nbin_low,nbin_high);
  TProfile* avgRadPosVsTime = (TProfile*)radPosVsTime->ProfileX();//"prof1",nbin_low,nbin_high);
  avgRadPosVsTime->Draw();
  c1->SaveAs("mean.eps");

  // plot the rms vs time 
  nbins = avgRadPosVsTime->GetXaxis()->GetNbins();
  x1    = avgRadPosVsTime->GetXaxis()->GetXmin();
  x2    = avgRadPosVsTime->GetXaxis()->GetXmax();
  TH1F* xRmsVsTime = new TH1F("xRmsVsTime","",nbins,x1,x2);
  
  for (int ibin(1); ibin < nbins; ibin++) {
    double binError = avgRadPosVsTime->GetBinError(ibin);
    double binc     = avgRadPosVsTime->GetBinEntries(ibin);
    double rms      = binError * pow(binc,0.5);
    if (rms > 0) xRmsVsTime->SetBinContent(ibin,rms);
  }

  xRmsVsTime->GetXaxis()->SetRangeUser(xmin,xmax);
  xRmsVsTime->Draw();
  c1->SaveAs(("rms"+app).c_str());
  
  TCanvas* c2 = new TCanvas("c2","",800,600);
  c2->SetBottomMargin(0.13);
  //gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  Float_t small = 1e-5;

  c2->Divide(1,2,small,small);
  c2->cd(1);
  gPad->SetBottomMargin(small);
  //gPad->SetLeftMargin(1.3);

  xRmsVsTime->SetLineWidth(2);
  xRmsVsTime->SetLineColor(kRed+1);
  xRmsVsTime->GetXaxis()->SetRangeUser(xmin,xmax);
  xRmsVsTime->GetYaxis()->SetTitle("#sigma(x) [mm]");
  xRmsVsTime->GetYaxis()->CenterTitle();
  //xRmsVsTime->GetYaxis()->SetRangeUser(21.01,25);
  xRmsVsTime->Draw();
  
  TLegend* legend1 = new TLegend(0.15,0.8,0.4,0.85);
  TLegendEntry* ent1 = legend1->AddEntry(xRmsVsTime,"#sigma(x)","");
  legend1->SetTextColor(kRed+1);
  legend1->Draw("SAME");

  c2->cd(2);
  gPad->SetTopMargin(small);
  gPad->SetBottomMargin(0.15); // change this number
  //gPad->SetLeftMargin(1.3);
  gPad->SetTickx();

  avgRadPosVsTime->SetLineWidth(2);
  avgRadPosVsTime->SetLineColor(kBlack);
  avgRadPosVsTime->SetMarkerColor(kBlack);
  avgRadPosVsTime->SetTitle("");
  avgRadPosVsTime->GetYaxis()->SetTitle("#bar{x} [mm]");
  avgRadPosVsTime->GetXaxis()->SetTitle("Time [#mus]");
  avgRadPosVsTime->GetXaxis()->CenterTitle();
  avgRadPosVsTime->GetXaxis()->SetRangeUser(xmin,xmax);
  //avgRadPosVsTime->GetYaxis()->SetRangeUser(-5,9.9);
  avgRadPosVsTime->GetYaxis()->CenterTitle();
  avgRadPosVsTime->Draw();
  TLegend* legend2 = new TLegend(0.15,0.85,0.4,0.9);
  TLegendEntry* ent2 = legend2->AddEntry(avgRadPosVsTime,"#bar{x}","");
  legend2->SetTextColor(kBlack);
  legend2->Draw("SAME");
  c2->SaveAs(("meanAndRms"+app).c_str());

}

void setStyle() {
  
  TStyle* myStyle = new TStyle("Plain", "myStyle");
  myStyle->cd();

  myStyle->SetCanvasColor(0);

  // stats box options
  //myStyle->SetOptStat("rem");
  myStyle->SetOptStat(0);
  myStyle->SetStatBorderSize(0);
  myStyle->SetStatY(0.85);
  myStyle->SetStatX(0.85);

  // legend options
  myStyle->SetLegendBorderSize(0);
  myStyle->SetLegendFont(62);
  myStyle->SetLegendTextSize(0.08);

  // title options
  myStyle->SetTitleTextColor(1); 
  //myStyle->SetTitleFont(62);
  //myStyle->SetTitleFontSize(0.05);
  myStyle->SetTitleBorderSize(0);
  myStyle->SetTitleX(0.5);
  myStyle->SetTitleAlign(23);

  // label options

  myStyle->SetTitleOffset(0.8,"x");
  myStyle->SetTitleOffset(0.8,"y");

  myStyle->SetTitleSize(0.06,"x");
  myStyle->SetTitleSize(0.06,"y");

  myStyle->SetTitleFont(62,"x");
  myStyle->SetTitleFont(62,"y");
  
  myStyle->SetLabelFont(62,"x");
  myStyle->SetLabelFont(62,"y");
  
  //myStyle->SetLabelSize(0.04,"x");
  //myStyle->SetLabelSize(0.04,"y");
  myStyle->SetLabelSize(0.05,"x");
  myStyle->SetLabelSize(0.05,"y");

  myStyle->SetTickLength(0.04,"XY");  

  // plotting options
  //myStyle->SetHistLineWidth(1.5);
  //myStyle->SetHistLineColor(kBlue+1);

  gROOT->ForceStyle();
  
}
