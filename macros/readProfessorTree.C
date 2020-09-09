// An example macro showing how to read data created by the ProfessorTreeDumper_module
// James Mott, Aug 29th 2016

#include<iostream>
#include<TFile.h>
#include<TTree.h>
#include<TH1D.h>
#include<TCanvas.h>

//Struct for tree defintion
struct StrawData {
  unsigned int run;
  unsigned int event;
  double triggerTime;
  int layer;
  int wire;
  double hitTime;
};

void readProfessorTree(){

  // Get tree from file (need to change these hard-coded paths if upstream module is changed)
  TFile* inFile = new TFile("Lab3TreeDump.root");
  TTree* strawTree = (TTree*)inFile->Get("professorTreeDumper/strawHits");
  TTree* scintTree = (TTree*)inFile->Get("professorTreeDumper/scintHits");

  // Create holders for data
  StrawData strawHit, scintHit;
  strawTree->SetBranchAddress("strawHits",&strawHit.run);
  scintTree->SetBranchAddress("scintHits",&scintHit.run);

  // Example histograms that we might fill
  TH1D* hStrawWire = new TH1D("hStrawWire","Wire of Straw Hit",32,-0.5,31.5);
  TH1D* hScintHitTime = new TH1D("hScintHitTime","Scintillator Hit Times", 120, 0, 12);

  // Loop over straw tree and fill histogram
  for(int iEntry = 0; iEntry < strawTree->GetEntries(); iEntry++){
    strawTree->GetEntry(iEntry);
    hStrawWire->Fill(strawHit.wire);
  }

  // Loop over scint tree and fill histogram
  for(int iEntry = 0; iEntry < scintTree->GetEntries(); iEntry++){
    scintTree->GetEntry(iEntry);
    hScintHitTime->Fill(scintHit.hitTime/1e6); //ns to ms
  }

  // Draw histograms
  TCanvas* cStrawWire = new TCanvas("cStrawWire","Straw Wire");
  hStrawWire->SetTitle(";Hit Wire;No. Entries");
  hStrawWire->Draw();

  TCanvas* cScintHitTime = new TCanvas("cScintHitTime","Scint Hit Time");
  hScintHitTime->SetTitle(";Hit Time [ms];No. Entries / 0.1 ms");
  hScintHitTime->Draw();

  
  return;

}
