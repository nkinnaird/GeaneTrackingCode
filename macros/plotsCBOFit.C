TH1D* RescaleAxis(TH1* input, Double_t Scale) {
  int bins = input->GetNbinsX();
  TAxis* xaxis = input->GetXaxis();
  double* ba = new double[bins+1];
  xaxis->GetLowEdge(ba);
  ba[bins] = ba[bins-1] + xaxis->GetBinWidth(bins);
  for (int i = 0; i < bins+1; i++) {
    ba[i] *= Scale;
  }
  TH1D* out = new TH1D(input->GetName(), input->GetTitle(), bins, ba);
  for (int i = 0; i <= bins; i++) {
    out->SetBinContent(i, input->GetBinContent(i));
    out->SetBinError(i, input->GetBinError(i));
  }
  return out;
}

double CBOFinalFunc(double* x, double* p){
  double t = x[0];
  //the offset, allowing it to vary with time
  double offset = p[0];
  double offsetGradient = p[1];
  //the exponetial part
  double ampExp = p[2];
  double offsetExp = p[3];
  double tcbo = p[4];
  //the sinusoidal wave
  double amp = p[5];
  double phase = p[6];
  double omega = p[7];
  double freqPower = p[8];

  double linear = offset + offsetGradient * t;
  double exponential = ampExp * exp( (offsetExp - t) / tcbo);
  //double sinusoid = amp * sin(TMath::TwoPi() * (t-phase) / omega); 
  double sinusoid = amp * sin(TMath::TwoPi() * (t-phase) / pow(omega,freqPower)); 
  
  //return linear;
  return linear + (exponential * sinusoid);
  
}


double CBOAmpFunc(double* x, double* p){

  double t = x[0];
  double tMin = p[0];
  double tMax = p[1];
  double CBOPeriod = p[2];
  double phase = p[3];
  double offset = p[4];

  if(t < tMin or t > tMax){
    return 0;
  }

  int segment = (t-tMin)/CBOPeriod;
  double A = p[2*segment+5]*sin(TMath::TwoPi()*(t-phase)/(p[2*segment+6]*CBOPeriod));
  return A + offset;

}

void plotsCBOFit() {

  //string folderName = "11169-77/";
  //string sample = "superbowl"; 

  //string folderName = "11589-92/";
  //string sample = "65kV";

  //for just single runs
  //string folderName = "11169-77/";
  //string sample = "11170";

  //string folderName = "11589-92/";
  //string sample = "11590"; 

  string folderName = "/gm2/app/users/jprice/MC1/offline/gm2Dev_v8_04_00/srcs/gm2tracker/fcl/Collab/";
  string sample = "stPatrick"; 

  string outFolder = "CBOFits";
  string ext = ".eps";
  
  TCanvas *c1 = new TCanvas("c1","c1",600,500);
  c1->SetRightMargin(0.115);

  gStyle->SetOptStat(000000);
  
  TFile *tf = new TFile((folderName + "trackRecoPlots_"+sample+".root").c_str());
  vector<string> TrackSummaryFolders = {"TrackSummary", "TrackSumary12", "TrackSumary18"};

  TLegend* stationLegend = new TLegend(0.7, 0.75, 0.85, 0.89);
  stationLegend->SetBorderSize(0);
  TLine* s12 = new TLine();
  TLine* s18 = new TLine();
  s12->SetLineColor(4);
  s18->SetLineColor(2);
  stationLegend->AddEntry(s12, "station 12", "l");
  stationLegend->AddEntry(s18, "station 18", "l");

  string TrackFolder = TrackSummaryFolders.at(0);
  TDirectory* dir = gFile->GetDirectory(Form("%s/FitResults", TrackFolder.c_str()));

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Mean & RMS 
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  for (auto stat : {"12", "18"}){
    string station = stat;

    // Position vs time
    dir = gFile->GetDirectory(("Extrapolation/vertices/station"+station+"/pValue>0.005_and_noVolumesHit").c_str());
    TH2F* h11;
    dir->GetObject("h_radialPos_vs_time",h11);
    h11->SetTitleOffset(1.5,"Y");
    h11->SetTitleOffset(1.2,"X");
    h11->SetTitle("");
    h11->GetXaxis()->SetTitle("Track Time [#mus]");
    h11->GetYaxis()->SetTitle("Radial Position [mm]");
    h11->GetXaxis()->SetRangeUser(30,70);
    h11->GetYaxis()->SetRangeUser(-70,70);
    //  h11->Draw("COLZ");
    //  c1->Update();
    //  c1->Print((outFolder+"/"+sample+"/RadialPosVsTime_s18"+ext).c_str());
    
    // Get mean radial pos vs time
    h11->GetXaxis()->UnZoom();
    
    //shift based o
    //rootManager_->Add(subDir, new TH2F("h_radialPos_vs_time"  ,"Radial pos vs mean time;Mean island time [us];Radial pos [mm]",6000,0,6000*0.148936,500,-250,250));
    //TH2F* h11_shifted = new TH2F("shifted", ";time after injection [#mus]; Radial Position [mm]", h11->GetXaxis()->GetNbins(), 0., 200., h11->GetYaxis()->GetNbins(), -70, 70);
    TH2F* h11_shifted =  (TH2F*)h11->Clone();

    //set every bin to 0 first
    for (int ix(1); ix < h11_shifted->GetXaxis()->GetNbins(); ix++){
      for (int iy(1); iy < h11_shifted->GetYaxis()->GetNbins(); iy++){
	h11_shifted->SetBinContent(ix,iy,0.0);
      }
    }

    double tCaloTrackerOffset = -1.040;
    double shiftTime = ( (24930.0 ) * 1.25 * 0.001) + tCaloTrackerOffset;
    
    for (int ix(1); ix < h11->GetXaxis()->GetNbins(); ix++){
      double t1 = h11->GetXaxis()->GetBinCenter(ix);
      double t2 = t1 - shiftTime;
      if (t2 < 0) continue;
      int bx = h11_shifted->GetXaxis()->FindBin(t2);
      //cout << bx << " new time for shifted: " << t2 << " original time: " << t1 << " from bin " << ix << "\n";
      for (int iy(1); iy < h11->GetYaxis()->GetNbins(); iy++){
	double y = h11->GetYaxis()->GetBinCenter(iy);
	int by = h11_shifted->GetYaxis()->FindBin(y);
	h11_shifted->SetBinContent(bx,by,h11->GetBinContent(ix,iy));
      }
    }

    //to not shift in time uncomment these
    //TProfile* avgRadPosVsTime = (TProfile*)h11->ProfileX();
    //shiftTime = 0.0;

    TProfile* avgRadPosVsTime = (TProfile*)h11_shifted->ProfileX();
    TProfile* avgRadPosVsTimeClone = (TProfile*)avgRadPosVsTime->Clone();

    //make long term CBO plot
    c1->Clear();
    c1->cd();
    avgRadPosVsTime->SetTitle(";Track Time [us];Radial Position Mean [mm]");
    avgRadPosVsTime->GetYaxis()->SetRangeUser(-7.99,32.99);
    avgRadPosVsTime->GetXaxis()->SetRangeUser(40 - shiftTime, 250-shiftTime);
    avgRadPosVsTime->SetMarkerStyle(7);
    avgRadPosVsTime->Draw("Hist P");
    //c1->SaveAs((outFolder+"/"+sample+"/RadialPosvsTime_Mean_s18"+ext).c_str());
    
    //Do FFT
    TH1 *hmCBO =0;
    TVirtualFFT::SetTransform(0);
    TH1F* xMeanVsTimeHist = (TH1F*)avgRadPosVsTime->ProjectionX("px1","");
    hmCBO = xMeanVsTimeHist->FFT(hmCBO,"MAG");
    TH1D* fftMagCBO = RescaleAxis(hmCBO,1./(xMeanVsTimeHist->GetXaxis()->GetXmax()*1e-3));
    fftMagCBO->SetTitle(";Frequency (kHz);Magnitude [Arb Units]");
    fftMagCBO->GetYaxis()->SetTitleOffset(1.2);
    fftMagCBO->GetXaxis()->SetTitleFont(42);
    fftMagCBO->GetYaxis()->SetTitleFont(42);
    fftMagCBO->GetXaxis()->SetLabelFont(42);
    fftMagCBO->GetYaxis()->SetLabelFont(42);
    fftMagCBO->Scale(0.001);
    fftMagCBO->GetXaxis()->SetRangeUser(0,1500);
    fftMagCBO->SetStats(0);
    fftMagCBO->Draw("HIST");
    c1->Update();
    c1->Print((outFolder+"/"+sample+"/FFT_CBOS"+station+ext).c_str());
    
    // Find maximum after 50 kHz and use this to set range
    fftMagCBO->GetXaxis()->SetRangeUser(50,1500);
    fftMagCBO->SetMaximum(fftMagCBO->GetMaximum()*1.1);
    fftMagCBO->GetXaxis()->SetRangeUser(0,1500);
    fftMagCBO->Draw("HIST");
    c1->Update();
    c1->Print((outFolder+"/"+sample+"/FFTZoomed_CBOS"+station+ext).c_str());
    
    // Find maximum after 50 kHz and use this to set range
    fftMagCBO->GetXaxis()->SetRangeUser(395,430);
    fftMagCBO->Draw("HIST");
    c1->Update();
    c1->Print((outFolder+"/"+sample+"/FFTVeryZoomed_CBOS"+station+ext).c_str());

    //do FFT of early and late times separately
    TH1 *hmCBO_early =0;
    TH1 *hmCBO_late  =0;
    TVirtualFFT::SetTransform(0);
    double width = xMeanVsTimeHist->GetBinWidth(2);
    //low: 40.0745 high: 933.542 width: 0.148936
    TH1F* xMeanVsTimeHist_early = new TH1F("early", "", 1000 ,40 - shiftTime,   (40 - shiftTime)+ (1000*width));
    TH1F* xMeanVsTimeHist_late  = new TH1F("late",  "", 1000 ,200 - shiftTime, (200 - shiftTime)+(1000*width));
    int nBins = xMeanVsTimeHist->GetXaxis()->GetNbins();
    int iEarly = 1;
    int iLate = 1;

    //loop over bins
    for (int i(1); i < nBins; i++){
      double cont = xMeanVsTimeHist->GetBinContent(i);
      double time = xMeanVsTimeHist->GetBinCenter(i);
      if (time >= xMeanVsTimeHist_early->GetXaxis()->GetXmin() && time <= xMeanVsTimeHist_early->GetXaxis()->GetXmax()) { 
	xMeanVsTimeHist_early->SetBinContent(iEarly, cont);
	iEarly++;
      }
      if (time >= xMeanVsTimeHist_late->GetXaxis()->GetXmin() && time <= xMeanVsTimeHist_late->GetXaxis()->GetXmax() ) { 
	xMeanVsTimeHist_late->SetBinContent(iLate, cont);
	iLate++;
      }
    }

    hmCBO_early = xMeanVsTimeHist_early->FFT(hmCBO_early,"MAG");
    //TH1D* fftMagCBO_early = RescaleAxis(hmCBO_early, 1./(xMeanVsTimeHist_early->GetXaxis()->GetXmax()*1e-3));
    TH1D* fftMagCBO_early = RescaleAxis(hmCBO_early, 1./( (xMeanVsTimeHist_early->GetXaxis()->GetXmax() - xMeanVsTimeHist_early->GetXaxis()->GetXmin() )*1e-3));
    fftMagCBO_early->SetTitle(";Frequency (kHz);Magnitude [Arb Units]");
    fftMagCBO_early->GetYaxis()->SetTitleOffset(1.2);
    fftMagCBO_early->GetXaxis()->SetTitleFont(42);
    fftMagCBO_early->GetYaxis()->SetTitleFont(42);
    fftMagCBO_early->GetXaxis()->SetLabelFont(42);
    fftMagCBO_early->GetYaxis()->SetLabelFont(42);
    fftMagCBO_early->Scale(0.001);
    fftMagCBO_early->GetXaxis()->SetRangeUser(100,800);
    fftMagCBO_early->SetStats(0);
    fftMagCBO_early->Draw("HIST");
    //xMeanVsTimeHist_early->Draw();
    c1->Update();
    c1->Print((outFolder+"/"+sample+"/FFT_CBOS"+station+"_early"+ext).c_str());
    
    fftMagCBO_early->GetXaxis()->SetRangeUser(395,430);
    fftMagCBO_early->Draw("HIST");
    c1->Update();
    c1->Print((outFolder+"/"+sample+"/FFT_CBOS"+station+"_early_zoomed"+ext).c_str());

    xMeanVsTimeHist_early->Draw();
    c1->Update();
    c1->Print((outFolder+"/"+sample+"/FFTInput_CBOS"+station+"_early"+ext).c_str());
    
    hmCBO_late = xMeanVsTimeHist_late->FFT(hmCBO_late,"MAG");
    //TH1D* fftMagCBO_late = RescaleAxis(hmCBO_late, 1./(xMeanVsTimeHist_late->GetXaxis()->GetXmax()*1e-3));
    TH1D* fftMagCBO_late = RescaleAxis(hmCBO_late, 1./( (xMeanVsTimeHist_late->GetXaxis()->GetXmax() - xMeanVsTimeHist_late->GetXaxis()->GetXmin() ) *1e-3));
    fftMagCBO_late->SetTitle(";Frequency (kHz);Magnitude [Arb Units]");
    fftMagCBO_late->GetYaxis()->SetTitleOffset(1.2);
    fftMagCBO_late->GetXaxis()->SetTitleFont(42);
    fftMagCBO_late->GetYaxis()->SetTitleFont(42);
    fftMagCBO_late->GetXaxis()->SetLabelFont(42);
    fftMagCBO_late->GetYaxis()->SetLabelFont(42);
    fftMagCBO_late->Scale(0.001);
    fftMagCBO_late->GetXaxis()->SetRangeUser(100,800);
    fftMagCBO_late->SetStats(0);
    fftMagCBO_late->Draw("HIST");
    //xMeanVsTimeHist_late->Draw();
    c1->Update();
    c1->Print((outFolder+"/"+sample+"/FFT_CBOS"+station+"_late"+ext).c_str());

    fftMagCBO_late->GetXaxis()->SetRangeUser(395,430);
    fftMagCBO_late->Draw("HIST");
    c1->Update();
    c1->Print((outFolder+"/"+sample+"/FFT_CBOS"+station+"_late_zoomed"+ext).c_str());

    xMeanVsTimeHist_late->Draw();
    c1->Update();
    c1->Print((outFolder+"/"+sample+"/FFTInput_CBOS"+station+"_late"+ext).c_str());
    
    ////////////////////////////////
    // Amplitude of CBO - intial fit
    ////////////////////////////////

    double CBOPeriod = 4 * 2.43168; //us - the factor of 4 determines how many oscillation to treat as constant amplitude for the initial fit
    //double CBOPeriod = 2.418379; //us
    //double CBOPeriod = 2.42; //us
    //double tMin(39.6), tMax(250);
    double tMin = 65.0 - shiftTime;
    double tMax = 250 - shiftTime;
    int nSegments = (tMax-tMin)/CBOPeriod;
    tMax = nSegments*CBOPeriod + tMin;
    TF1* CBOAmp = new TF1("CBOAmp", CBOAmpFunc, tMin, tMax, 2*nSegments + 5); // period, phase, start, stop, offset
    
    CBOAmp->FixParameter(0,tMin);
    CBOAmp->FixParameter(1,tMax);
    CBOAmp->FixParameter(2,CBOPeriod);
    
    // these phases were calculated for the unshifted hists, just intial parameters that get fit eventually
    double p = (station == "18")? 1.427 : 0.7618;
    double a = (station == "18")? 9.1 : 9.2;
    
    //cout << "PHASE shift due to time shift: " << fmod(shiftTime, (CBOPeriod / 4.0)) << "\n";
    p += (CBOPeriod/4.0) - fmod(shiftTime, (CBOPeriod / 4.0)); 

    // Station 18 phase and offset parameter
    CBOAmp->SetParameter(3,p);
    CBOAmp->SetParameter(4,a);
    
    for(int par = 5; par < CBOAmp->GetNpar(); par+=2){
      CBOAmp->SetParameter(par,8);
      CBOAmp->SetParameter(par+1,1.0 * 0.25);
      CBOAmp->SetParLimits(par+1,.9 * 0.25, 1.1 * 0.25);
    }
    CBOAmp->SetNpx(2000);
    avgRadPosVsTime->Draw("Hist P");
    CBOAmp->Draw("SAME");
    avgRadPosVsTime->Fit(CBOAmp,"R");
    
    cout << sample << " station " << stat << " " << CBOAmp->GetChisquare() << " ndf: " << CBOAmp->GetNDF() << " giving: " 
	 << CBOAmp->GetChisquare() / double(CBOAmp->GetNDF()) << "\n";
    
    //CBOAmpS1->Draw("SAME");;
    c1->SaveAs((outFolder+"/"+sample+"/CBOAmpS"+station+ext).c_str());
    avgRadPosVsTime->Draw("Hist E");
    CBOAmp->Draw("SAME");
    avgRadPosVsTime->GetXaxis()->SetRangeUser(70 - shiftTime, 80 - shiftTime);
    c1->SaveAs((outFolder+"/"+sample+"/CBOAmpS"+station+"_early"+ext).c_str());
    avgRadPosVsTime->GetXaxis()->SetRangeUser(200 - shiftTime, 210 - shiftTime);
    c1->SaveAs((outFolder+"/"+sample+"/CBOAmpS"+station+"_late"+ext).c_str());
    
    TGraphErrors* tgCBOAmp = new TGraphErrors();
    TGraphErrors* tgCBOPhase = new TGraphErrors();
    for(int i = 0; i < 2*nSegments; i+=2){
      int nPts = tgCBOAmp->GetN();
      tgCBOAmp->SetPoint(nPts, tMin+CBOPeriod*nPts, CBOAmp->GetParameter(5+i));
      tgCBOAmp->SetPointError(nPts, 0, CBOAmp->GetParError(5+i));
      tgCBOPhase->SetPoint(nPts, tMin+CBOPeriod*nPts, CBOAmp->GetParameter(6+i));
      tgCBOPhase->SetPointError(nPts, 0, CBOAmp->GetParError(6+i));
    }
    
    TCanvas* c2 = new TCanvas();
    tgCBOAmp->SetMarkerColor(1);
    tgCBOAmp->SetLineColor(1);
    tgCBOAmp->SetMarkerStyle(2);
    tgCBOAmp->GetYaxis()->SetRangeUser(0,20);
    tgCBOAmp->SetTitle(";Time [us];Radial CBO Amplitude [mm]");
    tgCBOAmp->Draw("AP");
    
    TF1* exp1 = new TF1("exp1","[0] * exp((x - [1])/(-1.0 *[2]))", 70.0 - shiftTime, 250.0 - shiftTime);
    exp1->SetParameter(0, 12.5);
    exp1->FixParameter(1, 70.0 - shiftTime);
    exp1->SetParameter(2, 250.0);
    exp1->SetParLimits(0, 5.0, 30.0);
    exp1->SetParLimits(2, 100,350.0);
    //exp1->FixParameter(0, 14.5);
    //exp1->FixParameter(2, 250.0);
    exp1->SetParName(0, "initial Amplitude");
    exp1->SetParName(1, "time offset");
    exp1->SetParName(2, "#tau_{CBO}");
    
    //tgCBOAmp->GetXaxis()->SetRangeUser(tMin - 10.0, 200.0);
    tgCBOAmp->Fit(exp1,"R");
    gStyle->SetOptFit(1111);
    //gStyle->SetOptStat(1);
    //tgCBOAmp->SetOptStat();
    c2->SaveAs((outFolder+"/"+sample+"/CBOAmpVsTimeS"+station+ext).c_str()); 
    
    TCanvas* c3 = new TCanvas();
    tgCBOPhase->SetMarkerColor(2);
    tgCBOPhase->SetLineColor(2);
    tgCBOPhase->SetMarkerStyle(2);
    tgCBOPhase->GetYaxis()->SetRangeUser(0.8,1.2);
    
    //convert from relative to absolute frequency
    for (int i(0); i < tgCBOPhase->GetN(); i++){
      double x, y;
      tgCBOPhase->GetPoint(i,x,y);
      double errorY = tgCBOPhase->GetErrorY(i);
      double fracError = errorY/y;
      double freq = 1000.0/(y * CBOPeriod); // in kHz
      tgCBOPhase->SetPoint(i,x,freq);
      tgCBOPhase->SetPointError(i,0,freq*fracError);
    }
    
    tgCBOPhase->SetTitle(";Time [us];Radial CBO Frequency [kHz]");
    tgCBOPhase->Draw("AP");
    c3->SaveAs((outFolder+"/"+sample+"/CBOPhaseS"+station+ext).c_str()); 
    
    // Calculate rms vs time 
    TH1F* xRmsVsTime = new TH1F("xRmsVsTime","",h11->GetNbinsX(),avgRadPosVsTime->GetXaxis()->GetXmin(),avgRadPosVsTime->GetXaxis()->GetXmax());
    for (int ibin(1); ibin < h11->GetNbinsX(); ibin++) {
      TH1F* projectY = (TH1F*)h11->ProjectionY("",ibin,ibin+1);
      double prms = projectY->GetRMS();
      if (prms > 0) xRmsVsTime->SetBinContent(ibin,prms);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    //now perform the final fit, where we use the output of the old fit, but allow the parameters to float
    //////////////////////////////////////////////////////////////////////////////////////////////

    TF1* CBOFinal = new TF1("CBOFinal", CBOFinalFunc, tMin, tMax, 9);
    //linear
    double linearOffset = CBOAmp->GetParameter(4);
    double linearGrad = 0.0;
    //exp
    double expA = 1.0;
    double expOffset = exp1->GetParameter(1);
    double tCBO = exp1->GetParameter(2);
    //sin
    double sinA = exp1->GetParameter(0);
    double phase = CBOAmp->GetParameter(3);
    //double freq = 1000.0/(y * CBOPeriod); // in kHz
    double omega = CBOAmp->GetParameter(6) * CBOPeriod; //initial frequency for now
    double freqPower = 1.0;

    CBOFinal->SetParName(0, "linear offs");
    //CBOFinal->FixParameter(0, linearOffset);
    CBOFinal->SetParameter(0, linearOffset);
    CBOFinal->SetParLimits(0, 0.99*linearOffset, 1.01*linearOffset);
    
    CBOFinal->SetParName(1, "linear grad");
    CBOFinal->FixParameter(1, linearGrad);

    CBOFinal->SetParName(2, "exp amp    ");
    CBOFinal->FixParameter(2, expA);
    
    CBOFinal->SetParName(3, "exp offset ");
    //CBOFinal->FixParameter(3, expOffset);
    CBOFinal->SetParameter(3, expOffset);
    CBOFinal->SetParLimits(3, 0.5*expOffset, 1.5*expOffset);
    
    CBOFinal->SetParName(4, "tau CBO    ");
    //CBOFinal->FixParameter(4, tCBO);
    CBOFinal->SetParameter(4, tCBO);
    CBOFinal->SetParLimits(4, 0.5*tCBO, 1.5*tCBO);

    CBOFinal->SetParName(5, "sin amp    ");
    //CBOFinal->FixParameter(5, sinA);
    CBOFinal->SetParameter(5, sinA);
    CBOFinal->SetParLimits(5, 0.5*sinA, 1.5*sinA);
    
    CBOFinal->SetParName(6, "sin phase  ");
    //CBOFinal->FixParameter(6, phase);
    CBOFinal->SetParameter(6, phase);
    CBOFinal->SetParLimits(6, 0.5*phase, 1.5*phase);
    
    CBOFinal->SetParName(7, "omega      ");
    //CBOFinal->FixParameter(7, omega);
    CBOFinal->SetParameter(7, omega);
    CBOFinal->SetParLimits(7, 0.99*omega, 1.01*omega);

    CBOFinal->SetParName(8, "omega power");
    //CBOFinal->FixParameter(8, freqPower);
    CBOFinal->SetParameter(8, freqPower);
    CBOFinal->SetParLimits(8, 0.95, 1.05);

    CBOFinal->SetNpx(10000);

    TCanvas* c4 = new TCanvas("c4", "for final CBO fit", 1600, 1200);
    c4->SetFillColor(10);
    c4->cd();
    TPad *pad1 = new TPad("pad1", "The pad 70% of the height",    0.0,  0.3, 1.0,  1.0, 10);
    TPad *pad2 = new TPad("pad2", "The pad 30% of the height, BL",0.0,  0.0, 0.33, 0.3, 10);
    TPad *pad3 = new TPad("pad3", "The pad 30% of the height, BM",0.33, 0.0, 0.66, 0.3, 10);
    TPad *pad4 = new TPad("pad4", "The pad 30% of the height, BR",0.66, 0.0, 1.0, 0.3, 10);
    //pad1->SetFillColor(10);
    //pad2->SetFillColor(10);

    pad1->Draw();
    pad2->Draw();
    pad3->Draw();
    pad4->Draw();

    pad1->cd();
    pad1->SetTopMargin(0.02);
    pad1->SetLeftMargin(0.08);
    pad1->SetRightMargin(0.02);
    avgRadPosVsTimeClone->GetXaxis()->SetRangeUser(tMin, tMax);
    avgRadPosVsTimeClone->GetYaxis()->SetRangeUser(-1.0, 23.0);
    avgRadPosVsTimeClone->SetTitle("");
    double range = tMax - tMin;
    TPaveText* pt = new TPaveText(tMin + (range/2.0) - 10.0, 20.0, tMin + (range/2.0) + 10.0, 22.5);
    pt->AddText(("Station " + station).c_str());
    pt->SetTextFont(gStyle->GetTextFont());
    pt->SetFillColor(10);
    avgRadPosVsTimeClone->GetYaxis()->SetTitle("Average Radial position [mm]");
    avgRadPosVsTimeClone->GetXaxis()->SetTitleSize(0.05);
    avgRadPosVsTimeClone->GetXaxis()->SetTitleOffset(0.8);
    avgRadPosVsTimeClone->GetYaxis()->SetTitleSize(0.05);
    avgRadPosVsTimeClone->GetYaxis()->SetTitleOffset(0.7);
    avgRadPosVsTimeClone->GetXaxis()->CenterTitle();
    avgRadPosVsTimeClone->GetYaxis()->CenterTitle();
    avgRadPosVsTimeClone->Fit(CBOFinal,"R");
    avgRadPosVsTimeClone->Draw();
    pt->Draw("SAME");

    pad2->cd();
    pad2->SetLeftMargin(0.04);
    pad2->SetRightMargin(0.0);
    pad2->SetTopMargin(0.01);
    TProfile* avgRadPosVsTimeClone2 = (TProfile*)avgRadPosVsTimeClone->Clone();
    avgRadPosVsTimeClone2->GetXaxis()->SetRangeUser(tMin, tMin+20.0);
    stringstream ss2;
    ss2 << std::setprecision(4) << tMin << " - " << tMin+20.0 << " [#mus]";
    TPaveText* pt2 = new TPaveText(tMin+2., 20.0, tMin+18., 22.5);
    pt2->AddText(ss2.str().c_str());
    pt2->SetFillColor(10);
    avgRadPosVsTimeClone2->GetXaxis()->SetTitle("");
    avgRadPosVsTimeClone2->Draw();
    pt2->Draw("SAME");

    pad3->cd();
    pad3->SetLeftMargin(0.0);
    pad3->SetRightMargin(0.0);
    pad3->SetTopMargin(0.01);
    TProfile* avgRadPosVsTimeClone3 = (TProfile*)avgRadPosVsTimeClone->Clone();
    avgRadPosVsTimeClone3->GetXaxis()->SetRangeUser((tMin + tMax)/2.0 , (tMin + tMax)/2.0 + 20.0);
    stringstream ss3;
    ss3 << std::setprecision(4) << (tMin + tMax)/2.0 << " - " << (tMin + tMax) /2.0 +20.0 << " [#mus]";
    TPaveText* pt3 = new TPaveText( (tMin+tMax)/2. +2., 20.0, (tMin+tMax)/2. +18., 22.5);
    pt3->AddText(ss3.str().c_str());
    pt3->SetFillColor(10);
    avgRadPosVsTimeClone3->GetXaxis()->SetTitle("");
    avgRadPosVsTimeClone3->Draw();
    pt3->Draw("SAME");

    pad4->cd();
    pad4->SetLeftMargin(0.0);
    pad4->SetRightMargin(0.02);
    pad4->SetTopMargin(0.01);
    TProfile* avgRadPosVsTimeClone4 = (TProfile*)avgRadPosVsTimeClone->Clone();
    avgRadPosVsTimeClone4->GetXaxis()->SetRangeUser( tMax - 20.0 , tMax);
    stringstream ss4;
    ss4 << std::setprecision(4) << tMax - 20.0 << " - " << tMax << " [#mus]";
    TPaveText* pt4 = new TPaveText(tMax-18., 20.0, tMax-2., 22.5);
    pt4->AddText(ss4.str().c_str());
    pt4->SetFillColor(10);
    avgRadPosVsTimeClone4->GetXaxis()->SetTitle("");
    avgRadPosVsTimeClone4->Draw();
    pt4->Draw("SAME");

    c4->SaveAs((outFolder+"/"+sample+"/FinalCBOAmpFit_"+station+ext).c_str());
  }
  return;

}
