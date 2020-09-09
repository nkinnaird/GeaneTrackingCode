#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <TFile.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TImage.h>
#include <TGraphErrors.h>
#include <TKey.h>
#include <TDirectory.h>
#include <sstream>


using namespace std;

int maxNumPlanes = 33;
std::vector<string> planeFolderNames_;

// To use input to root macro in command line do: root -l geanePlotsMacro.C+\(\"filepath.txt\"\) - the slashes indicate to bash which characters are used as given
int geanePlotsMacro(std::string filePath)
{


// open the file
   TFile *f = TFile::Open(filePath.c_str());

   if (f == 0) {
      // if we cannot open the file, print an error message and return immediatly
      printf("Error: cannot open file\n");
      return 0;
   }

  auto topDirNameChar = f->GetListOfKeys()->Last()->GetName(); 
  // this gets the top directory, first probably only works if there is one directory in the root file, last can work depending on how many plot folders there are and what you want to grab
  string topDirName = string(topDirNameChar);
  std::cout << "Top directory name: " << topDirName << std::endl;

  gStyle->SetOptStat(111111);
  gStyle->SetOptTitle(1);



  TH1F* tempChiSquaredHist = (TH1F*) f->Get("TrackSummary/FitResults/Chi2");
  std::cout << "Num entries: " << tempChiSquaredHist->GetEntries() << std::endl; 

  double UmpRMSarray[maxNumPlanes];
  double VmpRMSarray[maxNumPlanes];

  double UmtRMSarray[maxNumPlanes];
  double VmtRMSarray[maxNumPlanes];

  double UptRMSarray[maxNumPlanes];
  double VptRMSarray[maxNumPlanes];


  double truthPull1oPDistMean[maxNumPlanes];
  double truthPull1oPDistMeanErrors[maxNumPlanes];
  double truthPull1oPDistRMS[maxNumPlanes];
  double truthPull1oPDistRMSErrors[maxNumPlanes];

  double truthPullPuoPxDistMean[maxNumPlanes];
  double truthPullPuoPxDistMeanErrors[maxNumPlanes];
  double truthPullPuoPxDistRMS[maxNumPlanes];
  double truthPullPuoPxDistRMSErrors[maxNumPlanes];

  double truthPullPvoPxDistMean[maxNumPlanes];
  double truthPullPvoPxDistMeanErrors[maxNumPlanes];
  double truthPullPvoPxDistRMS[maxNumPlanes];
  double truthPullPvoPxDistRMSErrors[maxNumPlanes];

  double truthPullUDistMean[maxNumPlanes];
  double truthPullUDistMeanErrors[maxNumPlanes];
  double truthPullUDistRMS[maxNumPlanes];
  double truthPullUDistRMSErrors[maxNumPlanes];

  double truthPullVDistMean[maxNumPlanes];
  double truthPullVDistMeanErrors[maxNumPlanes];
  double truthPullVDistRMS[maxNumPlanes];
  double truthPullVDistRMSErrors[maxNumPlanes];
/////////////////////////////////////////////////////////////////////////////////////
  double measPullDistEntries[maxNumPlanes];
  double measPullDistMean[maxNumPlanes];
  double measPullDistMeanErrors[maxNumPlanes];
  double measPullDistRMS[maxNumPlanes];
  double measPullDistRMSErrors[maxNumPlanes];

  double measPullDistNumEntries[maxNumPlanes];
  double measPullDistNumMean[maxNumPlanes];
  double measPullDistNumMeanErrors[maxNumPlanes];
  double measPullDistNumRMS[maxNumPlanes];
  double measPullDistNumRMSErrors[maxNumPlanes];

  double measPullDistDenomEntries[maxNumPlanes];
  double measPullDistDenomMean[maxNumPlanes];
  double measPullDistDenomMeanErrors[maxNumPlanes];
  double measPullDistDenomRMS[maxNumPlanes];
  double measPullDistDenomRMSErrors[maxNumPlanes];

  double measPullDistX[33] = {0, 10, 15.195999999998, 30.196, 35.391999999998, 144.36310579611, 149.55910579611, 164.55910579611, 169.75510579611, 
    278.68219761525,  283.91345618201,  298.87819761525,  304.10945618201,  413.0710783996,  418.2670783996,  433.2670783996,  438.4630783996,  547.4207561905,  
    552.6167561905,  567.6167561905,  572.7537194832,  681.76061095307,  686.95661095307,  701.95661095307,  707.15261095307,  816.10295380575,  821.29895380575,  
    836.29895380575,  841.49495380575,  950.43279788281,  955.6287978828,  970.6287978828,  975.8247978828}; 
    // hardcoded for ease of use - if straw tracker geometry/placement changes then this will be off, but not by much - just for plotting so not really important

  double pointNo[maxNumPlanes];
  double zeros[maxNumPlanes];

  for (int i = 0; i < maxNumPlanes; ++i)
  {
    pointNo[i] = i;
    zeros[i] = 0;
  }


  for (int planeNum = 1; planeNum < maxNumPlanes; ++planeNum)
  {
    TH1F* tempMeasPull = (TH1F*) f->Get(("TrackSummary/TrackLength/Dist" + std::to_string(planeNum) + "/Measure Pulls/UV Measure Pull Dist " + std::to_string(planeNum)).c_str());
    TH1F* tempMeasPullNum = (TH1F*) f->Get(("TrackSummary/TrackLength/Dist" + std::to_string(planeNum) + "/Measure Pulls/Num Denom/Meas Pull Numerator Dist " + std::to_string(planeNum)).c_str());
    TH1F* tempMeasPullDenom = (TH1F*) f->Get(("TrackSummary/TrackLength/Dist" + std::to_string(planeNum) + "/Measure Pulls/Num Denom/Meas Pull Denominator Dist " + std::to_string(planeNum)).c_str());

/////////////////////////////////////////////////////////////////////////////////////    
    measPullDistEntries[planeNum] = tempMeasPull->GetEntries();
    measPullDistMean[planeNum] = tempMeasPull->GetMean();
    measPullDistMeanErrors[planeNum] = tempMeasPull->GetMeanError();
    measPullDistRMS[planeNum] = tempMeasPull->GetRMS();
    measPullDistRMSErrors[planeNum] = tempMeasPull->GetRMSError();
/////////////////////////////////////////////////////////////////////////////////////

    measPullDistNumEntries[planeNum] = tempMeasPullNum->GetEntries();

    measPullDistNumMean[planeNum] = tempMeasPullNum->GetMean();
    measPullDistNumMeanErrors[planeNum] = tempMeasPullNum->GetMeanError();

    measPullDistNumRMS[planeNum] = tempMeasPullNum->GetRMS();
    measPullDistNumRMSErrors[planeNum] = tempMeasPullNum->GetRMSError();

/////////////////////////////////////////////////////////////////////////////////////

    measPullDistDenomEntries[planeNum] = tempMeasPullDenom->GetEntries();

    measPullDistDenomMean[planeNum] = tempMeasPullDenom->GetMean();
    measPullDistDenomMeanErrors[planeNum] = tempMeasPullDenom->GetMeanError();

    measPullDistDenomRMS[planeNum] = tempMeasPullDenom->GetRMS();
    measPullDistDenomRMSErrors[planeNum] = tempMeasPullDenom->GetRMSError();

/////////////////////////////////////////////////////////////////////////////////////
  //code with slashes in names for now - hack to get around it, but then this seg faults on the second pass..

  // code for when running over plots without slashes in names
    TH1F* tempTruthPull1oP = (TH1F*) f->Get(("TrackSummary/TrackLength/Dist" + std::to_string(planeNum) + "/Truth Pulls/1P Truth Pull Dist " + std::to_string(planeNum)).c_str());
    TH1F* tempTruthPullPuoPz = (TH1F*) f->Get(("TrackSummary/TrackLength/Dist" + std::to_string(planeNum) + "/Truth Pulls/PuPz Truth Pull Dist " + std::to_string(planeNum)).c_str());
    TH1F* tempTruthPullPvoPz = (TH1F*) f->Get(("TrackSummary/TrackLength/Dist" + std::to_string(planeNum) + "/Truth Pulls/PvPz Truth Pull Dist " + std::to_string(planeNum)).c_str());

    TH1F* tempTruthPullU = (TH1F*) f->Get(("TrackSummary/TrackLength/Dist" + std::to_string(planeNum) + "/Truth Pulls/U Truth Pull Dist " + std::to_string(planeNum)).c_str());
    TH1F* tempTruthPullV = (TH1F*) f->Get(("TrackSummary/TrackLength/Dist" + std::to_string(planeNum) + "/Truth Pulls/V Truth Pull Dist " + std::to_string(planeNum)).c_str());

/////////////////////////////////////////////////////////////////////////////////////
    truthPull1oPDistMean[planeNum] = tempTruthPull1oP->GetMean();
    truthPull1oPDistMeanErrors[planeNum] = tempTruthPull1oP->GetMeanError();
    truthPull1oPDistRMS[planeNum] = tempTruthPull1oP->GetRMS();
    truthPull1oPDistRMSErrors[planeNum] = tempTruthPull1oP->GetRMSError();

    truthPullPuoPxDistMean[planeNum] = tempTruthPullPuoPz->GetMean();
    truthPullPuoPxDistMeanErrors[planeNum] = tempTruthPullPuoPz->GetMeanError();
    truthPullPuoPxDistRMS[planeNum] = tempTruthPullPuoPz->GetRMS();
    truthPullPuoPxDistRMSErrors[planeNum] = tempTruthPullPuoPz->GetRMSError();

    truthPullPvoPxDistMean[planeNum] = tempTruthPullPvoPz->GetMean();
    truthPullPvoPxDistMeanErrors[planeNum] = tempTruthPullPvoPz->GetMeanError();
    truthPullPvoPxDistRMS[planeNum] = tempTruthPullPvoPz->GetRMS();
    truthPullPvoPxDistRMSErrors[planeNum] = tempTruthPullPvoPz->GetRMSError();

    truthPullUDistMean[planeNum] = tempTruthPullU->GetMean();
    truthPullUDistMeanErrors[planeNum] = tempTruthPullU->GetMeanError();
    truthPullUDistRMS[planeNum] = tempTruthPullU->GetRMS();
    truthPullUDistRMSErrors[planeNum] = tempTruthPullU->GetRMSError();

    truthPullVDistMean[planeNum] = tempTruthPullV->GetMean();
    truthPullVDistMeanErrors[planeNum] = tempTruthPullV->GetMeanError();
    truthPullVDistRMS[planeNum] = tempTruthPullV->GetRMS();
    truthPullVDistRMSErrors[planeNum] = tempTruthPullV->GetRMSError();            

  // TH1F* tempChiSquaredHist = (TH1F*) f->Get(("iteration1/chiSquaredsPerPlane/ChiSquareds Planes Hit "+std::to_string(planeNum)).c_str()); 


  // tempChiSquaredHist->Draw();
  // TF1* chi2pdf = new TF1("chi2pdf","[2]*ROOT::Math::chisquared_pdf(x,[0],[1])",0,60);
  // chi2pdf->SetParameters(planeNum-5, 0., tempChiSquaredHist->Integral("WIDTH")); 

  // chi2pdf->Draw("SAME");
  // myCanvas->SaveAs(("Plots/chiSquaredsPerPlane/ChiSquaredsPlanesHit"+std::to_string(planeNum)+".png").c_str());



  } // planenum loop


/////////////////////////////////////////////////////////////////////////////////////

  TCanvas* myCanvas = new TCanvas("myCanvas","Canvas",200,10,1200,800);

  myCanvas->Divide(2,2);
  myCanvas->cd(1);


  TGraph* measEntriesPerDist = new TGraph(maxNumPlanes, measPullDistX, measPullDistEntries);
  measEntriesPerDist->SetTitle("meas pull Entries; Dist mm; Meas Pull Entries ");
  measEntriesPerDist->GetYaxis()->SetRangeUser(0,200000);
  measEntriesPerDist->SetMarkerStyle(20);
  measEntriesPerDist->SetMarkerColor(1);
  measEntriesPerDist->Draw("AP"); 

  myCanvas->cd(2);

  TGraphErrors* measMeanPerDist = new TGraphErrors(maxNumPlanes, measPullDistX, measPullDistMean,zeros,measPullDistMeanErrors);
  measMeanPerDist->SetTitle("meas pull Mean; Dist mm; Meas Pull Mean ");
  measMeanPerDist->GetYaxis()->SetRangeUser(-.05,.05);
  measMeanPerDist->SetMarkerStyle(20);
  measMeanPerDist->SetMarkerColor(1);
  measMeanPerDist->Draw("AP"); 

  TLine* zLine = new TLine(0,0,1100,0);
  zLine->SetLineStyle(2);
  zLine->SetLineWidth(1);
  zLine->Draw();

  myCanvas->cd(3);


  TGraphErrors* measRMSPerDist = new TGraphErrors(maxNumPlanes, measPullDistX, measPullDistRMS,zeros,measPullDistRMSErrors);
  measRMSPerDist->SetTitle("meas pull RMS; Dist mm; Meas Pull RMS ");
  measRMSPerDist->GetYaxis()->SetRangeUser(0.9,1.1);
  measRMSPerDist->SetMarkerStyle(20);
  measRMSPerDist->SetMarkerColor(1);
  measRMSPerDist->Draw("AP"); 

  TLine* oneLine = new TLine(0,1,1100,1);
  oneLine->SetLineStyle(2);
  oneLine->SetLineWidth(1);
  oneLine->Draw();


/////////////////////////////////////////////////////////////////////////////////////

  TCanvas* myCanvas2 = new TCanvas("myCanvas2","Canvas2",200,10,1200,800);


  myCanvas2->Divide(2,2);
  myCanvas2->cd(1);


  TGraph* measEntriesPerDistNum = new TGraph(maxNumPlanes, measPullDistX, measPullDistNumEntries);
  measEntriesPerDistNum->SetTitle("meas pull Num Entries; Dist mm; Meas Pull Num Entries ");
  measEntriesPerDistNum->GetYaxis()->SetRangeUser(0,200000);
  measEntriesPerDistNum->SetMarkerStyle(20);
  measEntriesPerDistNum->SetMarkerColor(1);
  measEntriesPerDistNum->Draw("AP"); 

  myCanvas2->cd(2);

  TGraphErrors* measMeanPerDistNum = new TGraphErrors(maxNumPlanes, measPullDistX, measPullDistNumMean,zeros,measPullDistNumMeanErrors);
  measMeanPerDistNum->SetTitle("meas pull Num Mean; Dist mm; Meas Pull Num Mean ");
  measMeanPerDistNum->GetYaxis()->SetRangeUser(-.01,.01);
  measMeanPerDistNum->SetMarkerStyle(20);
  measMeanPerDistNum->SetMarkerColor(1);
  measMeanPerDistNum->Draw("AP"); 

  zLine->Draw();

  myCanvas2->cd(3);


  TGraphErrors* measRMSPerDistNum = new TGraphErrors(maxNumPlanes, measPullDistX, measPullDistNumRMS,zeros,measPullDistNumRMSErrors);
  measRMSPerDistNum->SetTitle("meas pull Num RMS; Dist mm; Meas Pull Num RMS ");
  measRMSPerDistNum->GetYaxis()->SetRangeUser(.07,.17);
  measRMSPerDistNum->SetMarkerStyle(20);
  measRMSPerDistNum->SetMarkerColor(1);
  measRMSPerDistNum->Draw("AP"); 


/////////////////////////////////////////////////////////////////////////////////////

  TCanvas* myCanvas3 = new TCanvas("myCanvas3","Canvas3",200,10,1200,800);


  myCanvas3->Divide(2,2);
  myCanvas3->cd(1);


  TGraph* measEntriesPerDistDenom = new TGraph(maxNumPlanes, measPullDistX, measPullDistDenomEntries);
  measEntriesPerDistDenom->SetTitle("meas pull Denom Entries; Dist mm; Meas Pull Denom Entries ");
  measEntriesPerDistDenom->GetYaxis()->SetRangeUser(0,200000);
  measEntriesPerDistDenom->SetMarkerStyle(20);
  measEntriesPerDistDenom->SetMarkerColor(1);
  measEntriesPerDistDenom->Draw("AP"); 

  myCanvas3->cd(2);

  TGraphErrors* measMeanPerDistDenom = new TGraphErrors(maxNumPlanes, measPullDistX, measPullDistDenomMean,zeros,measPullDistDenomMeanErrors);
  measMeanPerDistDenom->SetTitle("meas pull Denom Mean; Dist mm; Meas Pull Denom Mean ");
  measMeanPerDistDenom->GetYaxis()->SetRangeUser(0,.2);
  measMeanPerDistDenom->SetMarkerStyle(20);
  measMeanPerDistDenom->SetMarkerColor(1);
  measMeanPerDistDenom->Draw("AP"); 

  myCanvas3->cd(3);


  TGraphErrors* measRMSPerDistDenom = new TGraphErrors(maxNumPlanes, measPullDistX, measPullDistDenomRMS,zeros,measPullDistDenomRMSErrors);
  measRMSPerDistDenom->SetTitle("meas pull Denom RMS; Dist mm; Meas Pull Denom RMS ");
  measRMSPerDistDenom->GetYaxis()->SetRangeUser(0,.04);
  measRMSPerDistDenom->SetMarkerStyle(20);
  measRMSPerDistDenom->SetMarkerColor(1);
  measRMSPerDistDenom->Draw("AP"); 

/////////////////////////////////////////////////////////////////////////////////////
  myCanvas3->cd(4);

  TGraphErrors* measRMSPerDistNumcopy = new TGraphErrors(maxNumPlanes, measPullDistX, measPullDistNumRMS,zeros,measPullDistNumRMSErrors);
  measRMSPerDistNumcopy->SetTitle("Num RMS and Denom Mean; Dist mm; Num RMS and Denom Mean");
  measRMSPerDistNumcopy->GetYaxis()->SetRangeUser(.07,.17);
  measRMSPerDistNumcopy->SetMarkerStyle(20);
  measRMSPerDistNumcopy->SetMarkerColor(2);
  measRMSPerDistNumcopy->SetLineColor(2);
  measRMSPerDistNumcopy->Draw("AP"); 

  measMeanPerDistDenom->Draw("PSAME"); 

/////////////////////////////////////////////////////////////////////////////////////

  TCanvas* myCanvas4 = new TCanvas("myCanvas4","Canvas4",200,10,1200,800);

  myCanvas4->Divide(2,2);
  myCanvas4->cd(1);

  TGraphErrors* truth1oPMeanPerDist = new TGraphErrors(maxNumPlanes, measPullDistX, truthPull1oPDistMean,zeros,truthPull1oPDistMeanErrors);
  truth1oPMeanPerDist->SetTitle("truth pull 1/P Mean; Dist mm; Truth Pull 1/P Mean ");
  truth1oPMeanPerDist->GetYaxis()->SetRangeUser(-.05,.05);
  truth1oPMeanPerDist->SetMarkerStyle(20);
  truth1oPMeanPerDist->SetMarkerColor(1);
  truth1oPMeanPerDist->Draw("AP"); 

  zLine->Draw();

  myCanvas4->cd(2);

  TGraphErrors* truth1oPRMSPerDist = new TGraphErrors(maxNumPlanes, measPullDistX, truthPull1oPDistRMS,zeros,truthPull1oPDistRMSErrors);
  truth1oPRMSPerDist->SetTitle("truth pull 1/P RMS; Dist mm; Truth Pull 1/P Pull RMS ");
  truth1oPRMSPerDist->GetYaxis()->SetRangeUser(0.9,1.1);
  truth1oPRMSPerDist->SetMarkerStyle(20);
  truth1oPRMSPerDist->SetMarkerColor(1);
  truth1oPRMSPerDist->Draw("AP"); 

  oneLine->Draw();

/////////////////////////////////////////////////////////////////////////////////////

  TCanvas* myCanvas5 = new TCanvas("myCanvas5","Canvas5",200,10,1200,800);

  myCanvas5->Divide(2,2);
  myCanvas5->cd(1);

  TGraphErrors* truthPuPxMeanPerDist = new TGraphErrors(maxNumPlanes, measPullDistX, truthPullPuoPxDistMean,zeros,truthPullPuoPxDistMeanErrors);
  truthPuPxMeanPerDist->SetTitle("truth pull Pu/Px Mean; Dist mm; Truth Pull Pu/Px Mean ");
  truthPuPxMeanPerDist->GetYaxis()->SetRangeUser(-.05,.05);
  truthPuPxMeanPerDist->SetMarkerStyle(20);
  truthPuPxMeanPerDist->SetMarkerColor(1);
  truthPuPxMeanPerDist->Draw("AP"); 

  zLine->Draw();

  myCanvas5->cd(2);

  TGraphErrors* truthPuPxRMSPerDist = new TGraphErrors(maxNumPlanes, measPullDistX, truthPullPuoPxDistRMS,zeros,truthPullPuoPxDistRMSErrors);
  truthPuPxRMSPerDist->SetTitle("truth pull Pu/Px RMS; Dist mm; Truth Pull Pu/Px Pull RMS ");
  truthPuPxRMSPerDist->GetYaxis()->SetRangeUser(0.9,1.1);
  truthPuPxRMSPerDist->SetMarkerStyle(20);
  truthPuPxRMSPerDist->SetMarkerColor(1);
  truthPuPxRMSPerDist->Draw("AP"); 

  oneLine->Draw();

  myCanvas5->cd(3);

  TGraphErrors* truthPvPxMeanPerDist = new TGraphErrors(maxNumPlanes, measPullDistX, truthPullPvoPxDistMean,zeros,truthPullPvoPxDistMeanErrors);
  truthPvPxMeanPerDist->SetTitle("truth pull Pv/Px Mean; Dist mm; Truth Pull Pv/Px Mean ");
  truthPvPxMeanPerDist->GetYaxis()->SetRangeUser(-.05,.05);
  truthPvPxMeanPerDist->SetMarkerStyle(20);
  truthPvPxMeanPerDist->SetMarkerColor(1);
  truthPvPxMeanPerDist->Draw("AP"); 

  zLine->Draw();

  myCanvas5->cd(4);

  TGraphErrors* truthPvPxRMSPerDist = new TGraphErrors(maxNumPlanes, measPullDistX, truthPullPvoPxDistRMS,zeros,truthPullPvoPxDistRMSErrors);
  truthPvPxRMSPerDist->SetTitle("truth pull Pv/Px RMS; Dist mm; Truth Pull Pv/Px Pull RMS ");
  truthPvPxRMSPerDist->GetYaxis()->SetRangeUser(0.9,1.1);
  truthPvPxRMSPerDist->SetMarkerStyle(20);
  truthPvPxRMSPerDist->SetMarkerColor(1);
  truthPvPxRMSPerDist->Draw("AP"); 

  oneLine->Draw();

///////////////////////////////////////////////////////////////////////////////////// 

  TCanvas* myCanvas6 = new TCanvas("myCanvas6","Canvas6",200,10,1200,800);

  myCanvas6->Divide(2,2);
  myCanvas6->cd(1);

  TGraphErrors* truthUMeanPerDist = new TGraphErrors(maxNumPlanes, measPullDistX, truthPullUDistMean,zeros,truthPullUDistMeanErrors);
  truthUMeanPerDist->SetTitle("truth pull U Mean; Dist mm; Truth Pull U Mean ");
  truthUMeanPerDist->GetYaxis()->SetRangeUser(-.05,.05);
  truthUMeanPerDist->SetMarkerStyle(20);
  truthUMeanPerDist->SetMarkerColor(1);
  truthUMeanPerDist->Draw("AP"); 

  zLine->Draw();

  myCanvas6->cd(2);

  TGraphErrors* truthURMSPerDist = new TGraphErrors(maxNumPlanes, measPullDistX, truthPullUDistRMS,zeros,truthPullUDistRMSErrors);
  truthURMSPerDist->SetTitle("truth pull U RMS; Dist mm; Truth Pull U Pull RMS ");
  truthURMSPerDist->GetYaxis()->SetRangeUser(0.9,1.1);
  truthURMSPerDist->SetMarkerStyle(20);
  truthURMSPerDist->SetMarkerColor(1);
  truthURMSPerDist->Draw("AP"); 

  oneLine->Draw();

  myCanvas6->cd(3);

  TGraphErrors* truthVMeanPerDist = new TGraphErrors(maxNumPlanes, measPullDistX, truthPullVDistMean,zeros,truthPullVDistMeanErrors);
  truthVMeanPerDist->SetTitle("truth pull V Mean; Dist mm; Truth Pull V Mean ");
  truthVMeanPerDist->GetYaxis()->SetRangeUser(-.05,.05);
  truthVMeanPerDist->SetMarkerStyle(20);
  truthVMeanPerDist->SetMarkerColor(1);
  truthVMeanPerDist->Draw("AP"); 

  zLine->Draw();

  myCanvas6->cd(4);

  TGraphErrors* truthVRMSPerDist = new TGraphErrors(maxNumPlanes, measPullDistX, truthPullVDistRMS,zeros,truthPullVDistRMSErrors);
  truthVRMSPerDist->SetTitle("truth pull V RMS; Dist mm; Truth Pull V Pull RMS ");
  truthVRMSPerDist->GetYaxis()->SetRangeUser(0.9,1.1);
  truthVRMSPerDist->SetMarkerStyle(20);
  truthVRMSPerDist->SetMarkerColor(1);
  truthVRMSPerDist->Draw("AP"); 

  oneLine->Draw();

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

/*
  // separate loop starting at 0 just copied from other temp macro
  for (int planeNum = 0; planeNum < maxNumPlanes; ++planeNum)
  {


    stringstream ss;
    ss << "Plane" << planeNum;
    planeFolderNames_.push_back( ss.str() );


    TH1F* tempHistUmp = (TH1F*) f->Get(("0nCut/PlanePlots/" + planeFolderNames_.at(planeNum)+ "/Measure Residuals/UresidualsMeasPred Plane " + std::to_string(planeNum)).c_str()); // Things can seg fault if names are not correct - be careful.
    TH1F* tempHistVmp = (TH1F*) f->Get(("0nCut/PlanePlots/" + planeFolderNames_.at(planeNum)+ "/Measure Residuals/VresidualsMeasPred Plane " + std::to_string(planeNum)).c_str()); 

    UmpRMSarray[planeNum] = tempHistUmp->GetRMS();
    VmpRMSarray[planeNum] = tempHistVmp->GetRMS();

    TH1F* tempHistUmt = (TH1F*) f->Get(("0nCut/PlanePlots/" + planeFolderNames_.at(planeNum)+ "/Measure Residuals/UresidualsMeasTruth Plane " + std::to_string(planeNum)).c_str()); // Things can seg fault if names are not correct - be careful.
    TH1F* tempHistVmt = (TH1F*) f->Get(("0nCut/PlanePlots/" + planeFolderNames_.at(planeNum)+ "/Measure Residuals/VresidualsMeasTruth Plane " + std::to_string(planeNum)).c_str()); 

    UmtRMSarray[planeNum] = tempHistUmt->GetRMS();
    VmtRMSarray[planeNum] = tempHistVmt->GetRMS();


    TH1F* tempHistUpt = (TH1F*) f->Get(("0nCut/PlanePlots/" + planeFolderNames_.at(planeNum)+ "/Truth Residuals/Absolute/U Truth Residual Abs Plane " + std::to_string(planeNum)).c_str()); // Things can seg fault if names are not correct - be careful.
    TH1F* tempHistVpt = (TH1F*) f->Get(("0nCut/PlanePlots/" + planeFolderNames_.at(planeNum)+ "/Truth Residuals/Absolute/V Truth Residual Abs Plane " + std::to_string(planeNum)).c_str()); 

    UptRMSarray[planeNum] = tempHistUpt->GetRMS();
    VptRMSarray[planeNum] = tempHistVpt->GetRMS();

  } // planenum loop

  

  TCanvas* myCanvasRMS = new TCanvas("myCanvasRMS","Canvas",200,10,1200,800);

  TGraph* RMSPerPlaneUmp = new TGraph(maxNumPlanes, pointNo, UmpRMSarray);
  RMSPerPlaneUmp->SetTitle("Meas-Pred RMS; Plane Number; RMS (mm)");
  // RMSPerPlaneUmp->GetXaxis()->SetRangeUser(0,8);
  RMSPerPlaneUmp->GetYaxis()->SetRangeUser(0.,.3);
  RMSPerPlaneUmp->SetMarkerStyle(20);
  RMSPerPlaneUmp->SetMarkerColor(1);
  RMSPerPlaneUmp->Draw("AP"); 

  TGraph* RMSPerPlaneVmp = new TGraph(maxNumPlanes, pointNo, VmpRMSarray);
  RMSPerPlaneVmp->SetMarkerStyle(20);
  RMSPerPlaneVmp->SetMarkerColor(2);
  RMSPerPlaneVmp->Draw("PSAME");

  TLegend* legUVmp = new TLegend(0.663,0.70,0.8,0.8);
  legUVmp->AddEntry(RMSPerPlaneUmp,"U plane","p");
  legUVmp->AddEntry(RMSPerPlaneVmp,"V plane","p");
  legUVmp->SetFillStyle(0);
  legUVmp->SetBorderSize(0);
  legUVmp->Draw();

  // myCanvasRMS->SaveAs("Plots/UVPlaneRMSmp.png");

  TCanvas* myCanvasRMS2 = new TCanvas("myCanvasRMS2","Canvas",200,10,1200,800);

  TGraph* RMSPerPlaneUmt = new TGraph(maxNumPlanes, pointNo, UmtRMSarray);
  RMSPerPlaneUmt->SetTitle("Meas-Truth RMS; Plane Number; RMS (mm)");
  // RMSPerPlaneUmt->GetXaxis()->SetRangeUser(0,8);
  RMSPerPlaneUmt->GetYaxis()->SetRangeUser(0.,.3);
  RMSPerPlaneUmt->SetMarkerStyle(20);
  RMSPerPlaneUmt->SetMarkerColor(1);
  RMSPerPlaneUmt->Draw("AP"); 

  TGraph* RMSPerPlaneVmt = new TGraph(maxNumPlanes, pointNo, VmtRMSarray);
  RMSPerPlaneVmt->SetMarkerStyle(20);
  RMSPerPlaneVmt->SetMarkerColor(2);
  RMSPerPlaneVmt->Draw("PSAME");

  TLegend* legUVmt = new TLegend(0.663,0.70,0.8,0.8);
  legUVmt->AddEntry(RMSPerPlaneUmt,"U plane","p");
  legUVmt->AddEntry(RMSPerPlaneVmt,"V plane","p");
  legUVmt->SetFillStyle(0);
  legUVmt->SetBorderSize(0);
  legUVmt->Draw();


  TCanvas* myCanvasRMS3 = new TCanvas("myCanvasRMS3","Canvas",200,10,1200,800);

  TGraph* RMSPerPlaneUpt = new TGraph(maxNumPlanes, pointNo, UptRMSarray);
  RMSPerPlaneUpt->SetTitle("Pred-Truth RMS; Plane Number; RMS (mm)");
  // RMSPerPlaneUpt->GetXaxis()->SetRangeUser(0,8);
  RMSPerPlaneUpt->GetYaxis()->SetRangeUser(0.,.3);
  RMSPerPlaneUpt->SetMarkerStyle(20);
  RMSPerPlaneUpt->SetMarkerColor(1);
  RMSPerPlaneUpt->Draw("AP"); 

  TGraph* RMSPerPlaneVpt = new TGraph(maxNumPlanes, pointNo, VptRMSarray);
  RMSPerPlaneVpt->SetMarkerStyle(20);
  RMSPerPlaneVpt->SetMarkerColor(2);
  RMSPerPlaneVpt->Draw("PSAME");

  TLegend* legUVpt = new TLegend(0.663,0.70,0.8,0.8);
  legUVpt->AddEntry(RMSPerPlaneUpt,"U plane","p");
  legUVpt->AddEntry(RMSPerPlaneVpt,"V plane","p");
  legUVpt->SetFillStyle(0);
  legUVpt->SetBorderSize(0);
  legUVpt->Draw();
*/


  // delete myCanvas;


  return 1;
}
