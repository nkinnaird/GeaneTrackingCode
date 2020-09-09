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

// taken from FT of number of tracks vs time plot
double fCBO = 415; // kHz

TH1F* shiftHisto(TH1F* hist, double xmin, double xmax);
TH1F* cutHisto(TH1F* hist, double xmin, double xmax);
TH1* FFTHisto(TProfile* hist);
TH1* FFTHisto(TH1F* hist);

int main(int argc, char* argv[]) {
  
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(62);
  gStyle->SetLegendTextSize(0.03);
  gStyle->SetPalette(55);
  //gStyle->SetOptStat("rem");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1112);
  gStyle->SetTitleOffset(1.2,"x");
  gStyle->SetTitleOffset(1.4,"y");
  
  int xmin =0.0;
  int xmax = 0.0;
  int nbins = 0.0;
  double binw = 0.0;
  int binshift = 0;
  double range = 0.0;
  int bin_range = 0;

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

  TH2F* xVsTime = (TH2F*)allEvents->Get("h_radialPos_vs_time");
  
  double cbo_period = 1.0 / (fCBO * 1000);
  cbo_period *= 1e6;
  binw = xVsTime->GetXaxis()->GetBinWidth(1);
  std::cout << "cbo_period = " << cbo_period << ", binw = " << binw << ", cbo_period/binw = " << cbo_period / binw << "\n";

  double x1 = -44.9;
  double x2 = 44.9;

  TH1F* radialPosHistAll = (TH1F*)xVsTime->ProjectionY();
  radialPosHistAll->Draw();
  c1->SaveAs("tmp1.eps");

  // plot of the mean radial pos vs time
  int nbin_low  = radialPosHistAll->FindBin(x1);
  int nbin_high = radialPosHistAll->FindBin(x2);

  int nbin_low2  = radialPosHistAll->FindBin(xmin);
  int nbin_high2 = radialPosHistAll->FindBin(xmax);

  TH1F* radialPosHistAll2 = (TH1F*)xVsTime->ProjectionY("proj1",nbin_low2,nbin_high2);
  radialPosHistAll2->Draw();
  c1->SaveAs("tmp2.eps");

  TProfile* xMeanVsTime = (TProfile*)xVsTime->ProfileX("prof1",nbin_low,nbin_high);

  // Project the TProfile onto x-axis in order to get a TH1 so we can shift to start at t=0
  TH1F* xMeanVsTimeHist = (TH1F*)xMeanVsTime->ProjectionX("px1","");
  xMeanVsTimeHist->GetXaxis()->SetRangeUser(50,150);

  // shift the histo to start at t = 0 and remove bins outside time range
  TH1F* xMeanVsTimeHistShift = shiftHisto(xMeanVsTimeHist,xmin,xmax);
  
  // plot the rms vs time 
  nbins = xMeanVsTime->GetXaxis()->GetNbins();
  x1    = xMeanVsTime->GetXaxis()->GetXmin();
  x2    = xMeanVsTime->GetXaxis()->GetXmax();
  TH1F* xRmsVsTime = new TH1F("xRmsVsTime","",nbins,x1,x2);
  
  for (int ibin(1); ibin < nbins; ibin++) {
    double binError = xMeanVsTime->GetBinError(ibin);
    double binc     = xMeanVsTime->GetBinEntries(ibin);
    double rms      = binError * pow(binc,0.5);
    if (rms > 0) xRmsVsTime->SetBinContent(ibin,rms);
  }

  // shift to start at t = 0 and remove bins outside time range
  TH1F* xRmsVsTimeNew = shiftHisto(xRmsVsTime,xmin,xmax);
  
  // perform FFT
  TH1* xRmsVsTime_FFT_scale = FFTHisto(xRmsVsTimeNew);
  xRmsVsTime_FFT_scale->GetXaxis()->SetTitle("Frequency [MHz]");
  xRmsVsTime_FFT_scale->GetXaxis()->SetRangeUser(xRmsVsTime_FFT_scale->GetBinLowEdge(2),1.5);
  xRmsVsTime_FFT_scale->GetXaxis()->CenterTitle();
  xRmsVsTime_FFT_scale->GetYaxis()->CenterTitle();
  xRmsVsTime_FFT_scale->SetTitle("RMS of radial position FFT");
  xRmsVsTime_FFT_scale->Draw("HIST");
  c1->SaveAs(("xRmsVsTime_FFT_scale"+app).c_str());
  
  xRmsVsTime_FFT_scale->GetXaxis()->SetRangeUser(0,3.5);
  xRmsVsTime_FFT_scale->GetYaxis()->SetRangeUser(0,30);
  xRmsVsTime_FFT_scale->GetXaxis()->SetTitle("Frequency [MHz]");
  xRmsVsTime_FFT_scale->GetXaxis()->CenterTitle();
  xRmsVsTime_FFT_scale->GetYaxis()->CenterTitle();
  xRmsVsTime_FFT_scale->GetYaxis()->SetTitle("");
  xRmsVsTime_FFT_scale->Draw("HIST");
  c1->SaveAs(("xRmsVsTime_FFT_scale_zoom"+app).c_str());
  xRmsVsTime_FFT_scale->Delete();
  
  // FFT of the rms plot to get frequency for modulo plot
  c1->Clear();

  TH1F* hist_new = (TH1F*)xMeanVsTimeHistShift->Clone();
  double total = 0.0;
  for (int ibin(0); ibin <= hist_new->GetNbinsX(); ibin++) {
    total+=xMeanVsTimeHistShift->GetBinContent(ibin);
  }
  
  double mean = 0.0;//total / hist_new->GetNbinsX();
  for (int ibin(0); ibin <= hist_new->GetNbinsX(); ibin++) {
    hist_new->SetBinContent(ibin,hist_new->GetBinContent(ibin)-mean);
  }

  //hist_new->Draw();
  //c1->SaveAs("hist_new.eps");
    
  //TH1* xMeanVsTime_FFT_scale = FFTHisto(xMeanVsTimeHistShift);
  TH1* xMeanVsTime_FFT_scale = FFTHisto(hist_new);
  
  xMeanVsTime_FFT_scale->GetXaxis()->SetTitle("Frequency [MHz]");
  xMeanVsTime_FFT_scale->GetXaxis()->CenterTitle();
  xMeanVsTime_FFT_scale->GetXaxis()->SetRangeUser(xMeanVsTime_FFT_scale->GetBinLowEdge(2),1.5);
  xMeanVsTime_FFT_scale->GetYaxis()->CenterTitle();
  xMeanVsTime_FFT_scale->SetTitle("Mean of radial position FFT");
  xMeanVsTime_FFT_scale->SetMinimum(0);
  xMeanVsTime_FFT_scale->Draw("HIST");
  c1->SaveAs(("xMeanVsTime_FFT_scale"+app).c_str());
  
  xMeanVsTime_FFT_scale->GetXaxis()->SetRangeUser(0,3.5);
  xMeanVsTime_FFT_scale->GetYaxis()->SetRangeUser(0,30);
  xMeanVsTime_FFT_scale->GetXaxis()->SetTitle("Frequency [MHz]");
  xMeanVsTime_FFT_scale->GetXaxis()->CenterTitle();
  xMeanVsTime_FFT_scale->GetYaxis()->CenterTitle();
  xMeanVsTime_FFT_scale->Draw();
  c1->SaveAs(("xMeanVsTime_FFT_scale_zoom"+app).c_str());
  xMeanVsTime_FFT_scale->Delete();

}

TH1F* cutHisto(TH1F* hist, double xmin, double xmax) {
  
  TH1F* hist2 = (TH1F*)hist->Clone();

  int bin1 = hist2->FindBin(xmin);
  std::cout << "bin1 = " << bin1 << "\n";

  int bin2 = hist2->FindBin(xmax);
  std::cout << "bin2 = " << bin2 << "\n";
  
  for (int ibin(1); ibin <= hist2->GetNbinsX(); ibin++) {
    if (ibin < bin1 || ibin > bin2) hist2->SetBinContent(ibin,0);
  }

  return hist2;
}

TH1F* shiftHisto(TH1F* hist, double xmin, double xmax) {

  // remove empty bins not in range
  int binshift = 0;
  for (int ibin(0); ibin < hist->GetXaxis()->GetNbins(); ibin++) {
    double binc = hist->GetBinContent(ibin);
    if (binc > 0) {
      binshift = ibin;
      break;
    }
  }

  double binw = hist->GetXaxis()->GetBinWidth(1);
  double range = xmax - xmin;
  int bin_range = int(range / binw);
  std::cout << bin_range;
  TH1F* histNew = new TH1F("hnew","",bin_range,0,range);
  for (int ibin(binshift); ibin < (binshift + bin_range); ibin++) {
    double binc = hist->GetBinContent(ibin);
    histNew->SetBinContent(ibin-binshift+1,binc);    
  }
  return histNew;
}

TH1* FFTHisto(TProfile* hist) {

  // FFT of the rms plot to get frequency for modulo plot
  TH1* hist_FFT = 0;
  TVirtualFFT::SetTransform(0);
  hist_FFT = (TH1*)hist->FFT(hist_FFT,"MAG");

  int nbins = hist_FFT->GetNbinsX();
  double binw  = hist->GetBinWidth(1);
  double freq = 1/binw;
  
  double xfft1 = hist_FFT->GetXaxis()->GetXmin();
  double xfft2 = hist_FFT->GetXaxis()->GetXmax();
  
  TH1F* hist_FFT_scale = new TH1F("hfft","",nbins/2 + 1,xfft1,xfft2);
  for (int i(1); i<=( nbins/2 + 1); i++) {
    double y0 = (hist_FFT->GetBinContent(i) - hist_FFT->GetBinContent(nbins+1 - i));
    double y = sqrt(y0*y0);
    double ynew = y/sqrt(nbins);
    double x = hist_FFT->GetXaxis()->GetBinCenter(i);
    hist_FFT_scale->Fill(x,ynew);
  }
  
  hist_FFT->Delete();
  hist_FFT_scale->GetXaxis()->SetLimits(0, freq);
  return hist_FFT_scale;

}

TH1* FFTHisto(TH1F* hist) {

  // FFT of the rms plot to get frequency for modulo plot
  TH1* hist_FFT = 0;
  TVirtualFFT::SetTransform(0);
  hist_FFT = (TH1*)hist->FFT(hist_FFT,"MAG");

  int nbins = hist_FFT->GetNbinsX();
  double binw  = hist->GetBinWidth(1);
  double freq = 1/binw;
  
  double xfft1 = hist_FFT->GetXaxis()->GetXmin();
  double xfft2 = hist_FFT->GetXaxis()->GetXmax();
  
  TH1F* hist_FFT_scale = new TH1F("hfft","",nbins/2 + 1,xfft1,xfft2);
  for (int i(1); i<=( nbins/2 + 1); i++) {
    double y0 = (hist_FFT->GetBinContent(i) - hist_FFT->GetBinContent(nbins+1 - i));
    double y = sqrt(y0*y0);
    double ynew = y/sqrt(nbins);
    double x = hist_FFT->GetXaxis()->GetBinCenter(i);
    hist_FFT_scale->Fill(x,ynew);
  }

  hist_FFT->Delete();
  hist_FFT_scale->GetXaxis()->SetLimits(0, freq);

  return hist_FFT_scale;

}
