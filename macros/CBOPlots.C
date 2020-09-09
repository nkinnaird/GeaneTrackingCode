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

TH2F* moduloPlot(TH2F* hist, double period, double xmin, double xmax) {

  int nbinsx = hist->GetXaxis()->GetNbins();
  int nbinsy = hist->GetYaxis()->GetNbins();
  double maxy = hist->GetYaxis()->GetBinUpEdge(nbinsy);
  double miny = hist->GetYaxis()->GetBinLowEdge(1);
  
  double binw = hist->GetXaxis()->GetBinWidth(1);
  double t_tot = nbinsx * binw;
  double frac = t_tot / period;
  int n_osc_tot = int(frac / binw);
  int nbins_mod = int(period / binw);
  
  int n_osc(1);
  int newBin = 1;
  
  TH2F* h2mod = new TH2F("h2mod","",nbins_mod,0,nbins_mod*binw,nbinsy,miny,maxy);

  for (int ibinx(1); ibinx < nbinsx; ibinx++) {
    for (int ibiny(1); ibiny < nbinsy; ibiny++) {
      int binc = hist->GetBinContent(ibinx,ibiny);
      double x = hist->GetXaxis()->GetBinCenter(ibinx);
      double y = hist->GetYaxis()->GetBinCenter(ibiny);

      if (x < xmin || x > xmax) continue;

      int current_binc = h2mod->GetBinContent(newBin,ibiny);
      h2mod->SetBinContent(newBin,ibiny,current_binc+binc);
    }
    if (ibinx % (nbins_mod * n_osc) == 0) {
      newBin = 0;
      n_osc++;
    }
    newBin++;
  }

  //  TProfile* h = (TProfile*)h2mod->ProfileX();
  return h2mod;

}

TH1F* moduloPlot(TH1F* hist, double period, double xmin, double xmax) {
  
  int nbins = hist->GetXaxis()->GetNbins();
  
  double binw = hist->GetXaxis()->GetBinWidth(1);
  double t_tot = nbins * binw;
  double frac = t_tot / period;
  int n_osc_tot = int(frac / binw);
  int nbins_mod = int(period / binw);
  
  TH1F* hmod = new TH1F("hmod","",nbins_mod,0,period);

  int n_osc(1);
  int newBin = 1;
  
  for (int ibin(1); ibin < nbins; ibin++) {
    int binc = hist->GetBinContent(ibin);
    double x = hist->GetXaxis()->GetBinCenter(ibin);

    if (x < xmin || x > xmax) continue;
    
    int current_binc = hmod->GetBinContent(newBin);
    hmod->SetBinContent(newBin,current_binc+binc);

    if (ibin % (nbins_mod * n_osc) == 0) {
      newBin = 0;
      n_osc++;
    }
    newBin++;    
  }
  return hmod;
}

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
  TFile* f  = new TFile("/gm2/data/g2be/Production/Plots/SummaryPlots/trackRecoPlots_GoodRuns.root");

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

  // project the radial pos vs time histo onto the position axis in the storage region only
  double x1 = -44.9;
  double x2 = 44.9;
  TH1F* radialPosHistAll = (TH1F*)radPosVsTime->ProjectionY();
  int nbin_low  = radialPosHistAll->FindBin(x1);
  int nbin_high = radialPosHistAll->FindBin(x2);
  std::cout << "nbin_low = " << nbin_low << ", nbin_high = " << nbin_high << "\n";
  std::cout << "radialPosHistAll->GetXaxis()->GetBinCenter(nbin_low)  = " << radialPosHistAll->GetXaxis()->GetBinCenter(nbin_low)  << "\n";
  std::cout << "radialPosHistAll->GetXaxis()->GetBinCenter(nbin_high) = " << radialPosHistAll->GetXaxis()->GetBinCenter(nbin_high) << "\n";
  TH1F* radialPosHist = (TH1F*)radPosVsTime->ProjectionY("proj1",nbin_low,nbin_high);
  
  radialPosHist->Draw();
  c1->SaveAs("tmp2.eps");

  radialPosHistAll->Draw();
  c1->SaveAs("tmp0.eps");
  
  TH1F* radialPosHist_rebin = (TH1F*)radialPosHist->Clone();
  radialPosHist_rebin->Rebin(5);
  int nbin_low_rebin  = radialPosHist_rebin->FindBin(x1);
  int nbin_high_rebin = radialPosHist_rebin->FindBin(x2);
  double y1 = radialPosHist_rebin->GetBinContent(nbin_low_rebin);
  double y2 = radialPosHist_rebin->GetBinContent(nbin_high_rebin);
  std::cout << "y1 = " << y1 << ", y2 = " << y2 << "\n";
  TLine* l1 = new TLine(x1,0,x1,y1);
  TLine* l2 = new TLine(x2,0,x2,y2);
  l1->SetLineColor(kRed+1);
  l1->SetLineStyle(2);
  l1->SetLineWidth(2);
  l2->SetLineColor(kRed+1);
  l2->SetLineStyle(2);
  l2->SetLineWidth(2);
  TH1F* radialPos_lostMuons = (TH1F*)radialPosHist_rebin->Clone();
  for (int ibin(0); ibin < radialPos_lostMuons->GetXaxis()->GetNbins(); ibin++) {
    if (ibin > nbin_low_rebin) radialPos_lostMuons->SetBinContent(ibin,0);
  }
  radialPos_lostMuons->SetFillStyle(3344);
  radialPos_lostMuons->SetFillColor(kBlue+1);
  
  TLegend* leg = new TLegend(0.15,0.6,0.4,0.8);
  TLegendEntry* e1 = leg->AddEntry(l1,"Storage region boundary","l");
  TLegendEntry* e2 = leg->AddEntry(radialPos_lostMuons,"Low radius","lf");
  e1->SetTextColor(kRed+1);  
  e2->SetTextColor(kBlue+1);  
  radialPosHist_rebin->GetXaxis()->SetTitle("x [mm]");
  radialPosHist_rebin->GetYaxis()->SetTitle("Entries");
  radialPosHist_rebin->SetTitle("Radial position");
  radialPosHist_rebin->Draw();
  leg->Draw("SAME");
  radialPos_lostMuons->GetXaxis()->SetTitle("x [mm]");
  radialPos_lostMuons->GetYaxis()->SetTitle("Entries");
  radialPos_lostMuons->SetTitle("Radial position");
  radialPos_lostMuons->Draw("SAME");
  l1->Draw("SAME");
  l2->Draw("SAME");
  c1->SaveAs(("radialPos_" + app).c_str());

  //TH2F* tmp1 = (TH2F*)radPosVsTime->Clone();
  //TH1F* tmp2 = (TH1F*)tmp1->ProjectionY("proj1",nbin_low,nbin_high);
  //tmp2->Rebin(10);
  //tmp2->Draw();
  //c1->SaveAs(("radialPos_"+app).c_str());

  TProfile* avgRadPosVsTime = (TProfile*)radPosVsTime->ProfileX("prof1",nbin_low,nbin_high);
  TProfile* avgRadPosVsTime_lostMuons = (TProfile*)radPosVsTime->ProfileX("prof2",0,nbin_low);

  TH2F* radPosVsTimeMod = moduloPlot(radPosVsTime,cbo_period,xmin,xmax);
  TProfile* avgRadPosVsTimeMod = (TProfile*)radPosVsTimeMod->ProfileX();
  
  // rebin in CBO period
  int rebin_factor = int(cbo_period/binw);
  TProfile* avgRadPosVsTime_rebin = (TProfile*)avgRadPosVsTime->Clone();

  avgRadPosVsTime_rebin->Rebin(rebin_factor);
  avgRadPosVsTime_rebin->GetXaxis()->SetRangeUser(xmin,xmax);
  avgRadPosVsTime->GetXaxis()->SetRangeUser(xmin,xmax);
  avgRadPosVsTime->GetYaxis()->SetRangeUser(-200,50);

  avgRadPosVsTime->SetLineColor(kRed+1);
  avgRadPosVsTime->SetMarkerColor(kRed+1);
  avgRadPosVsTime_lostMuons->GetXaxis()->SetRangeUser(xmin,xmax);
  avgRadPosVsTime_lostMuons->SetLineColor(kBlue+1);
  avgRadPosVsTime_lostMuons->SetMarkerColor(kBlue+1);

  TLegend* leg2 = new TLegend(0.15,0.5,0.4,0.6);
  TLegendEntry* e21 = leg2->AddEntry(avgRadPosVsTime,"Stored e+","");
  TLegendEntry* e22 = leg2->AddEntry(avgRadPosVsTime_lostMuons,"Low radius","");
  e21->SetTextColor(kRed+1);  
  e22->SetTextColor(kBlue+1);  
  avgRadPosVsTime->GetXaxis()->SetTitle("Time [#mus]");
  avgRadPosVsTime->GetYaxis()->SetTitle("");
  avgRadPosVsTime_lostMuons->GetXaxis()->SetTitle("Time [#mus]");
  avgRadPosVsTime_lostMuons->GetYaxis()->SetTitle("");
  avgRadPosVsTime->SetTitle("Average radial pos vs time, `good' runs");
  avgRadPosVsTime_lostMuons->SetTitle("Average radial pos vs time, 'good' runs");
  avgRadPosVsTime->Draw();
  leg2->Draw("SAME");
  avgRadPosVsTime_lostMuons->Draw("SAME");
  c1->SaveAs(("avgRadPosVsTime_posiAndMuon"+app).c_str());

  avgRadPosVsTime->SetLineColor(kBlue+1);
  avgRadPosVsTime->SetMarkerColor(kBlue+1);
  avgRadPosVsTime->GetXaxis()->SetRangeUser(xmin,xmax);
  avgRadPosVsTime->Draw();
  c1->SaveAs(("avgRadPosVsTime"+app).c_str());

  avgRadPosVsTime_lostMuons->GetXaxis()->SetRangeUser(xmin,xmax);
  avgRadPosVsTime_lostMuons->Draw();
  c1->SaveAs(("avgRadPosVsTime_lostMuons"+app).c_str());

  radPosVsTime->GetXaxis()->SetRangeUser(xmin,xmax);
  radPosVsTime->Draw();
  c1->SaveAs(("radPosVsTime"+app).c_str());

  //radPosVsTime_rebin->GetXaxis()->SetRangeUser(100,300);
  avgRadPosVsTime_rebin->Draw();
  c1->SaveAs(("avgRadPosVsTime_rebin"+app).c_str());

  radPosVsTimeMod->Draw();
  c1->SaveAs(("radPosVsTimeMod"+app).c_str());

  avgRadPosVsTimeMod->Draw();
  c1->SaveAs(("avgRadPosVsTimeMod"+app).c_str());
  
  avgRadPosVsTime_rebin->SetLineColor(kRed+1);
  avgRadPosVsTime->Draw();
  avgRadPosVsTime_rebin->Draw("SAME");
  c1->SaveAs(("both"+app).c_str());
  
  // plot diff as fn of time

  nbins = avgRadPosVsTime->GetXaxis()->GetNbins();
  x1 = 0;
  x2 = avgRadPosVsTime->GetXaxis()->GetBinUpEdge(nbins);
  TH1F* pos_diff = new TH1F("pos_diff","pos_diff",nbins,x1,x2);
  for (int ibin(0); ibin < nbins; ibin++) {
    double time = avgRadPosVsTime->GetBinCenter(ibin);
    if (time < cbo_period) continue;
    double pos1 = avgRadPosVsTime->GetBinContent(avgRadPosVsTime->FindBin(time));
    double pos2 = avgRadPosVsTime_rebin->GetBinContent(avgRadPosVsTime_rebin->FindBin(time));
    double diff = pos1 - pos2;
    //    std::cout << "pos1 = " << pos1 << ", pos2 = " << pos2 << ", diff = " << diff << ", time = " << time << "\n";
    pos_diff->SetBinContent(ibin,diff);
  }

  pos_diff->GetXaxis()->SetRangeUser(xmin,xmax);
  pos_diff->Draw();
  c1->SaveAs(("pos_diff"+app).c_str());
  
  TH1F* pos_diff_mod = moduloPlot(pos_diff,cbo_period,xmin,xmax);
  pos_diff_mod->Draw();
  c1->SaveAs(("pos_diff_mod"+app).c_str());

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

  myStyle->SetTitleOffset(1.1,"x");
  myStyle->SetTitleOffset(1.1,"y");

  myStyle->SetTitleSize(0.04,"x");
  myStyle->SetTitleSize(0.04,"y");

  myStyle->SetTitleFont(62,"x");
  myStyle->SetTitleFont(62,"y");
  
  myStyle->SetLabelFont(62,"x");
  myStyle->SetLabelFont(62,"y");
  
  //myStyle->SetLabelSize(0.04,"x");
  //myStyle->SetLabelSize(0.04,"y");
  myStyle->SetLabelSize(0.04,"x");
  myStyle->SetLabelSize(0.04,"y");

  myStyle->SetTickLength(0.02,"XY");  

  // plotting options
  //myStyle->SetHistLineWidth(1.5);
  //myStyle->SetHistLineColor(kBlue+1);

  gROOT->ForceStyle();
  
}
