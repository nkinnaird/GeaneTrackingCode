double CBOWidthFinalFunc(double* x, double* p){

  double t = x[0];
  //the offset, allowing it to vary with time
  double offset = p[0];
  double offsetGradient = p[1];
  //the exponetial part
  double ampExp = p[2];
  double offsetExp = p[3];
  double tcbo = p[4];

  //the main sinusoidal wave
  double amp = p[5];
  double phase = p[6];
  double omega = p[7];
  double freqPower = p[8];
  
  //the second sinusoidal wave amplitude
  double amplitudeRatio = p[9];
  
  //additional offset for the second sin wave
  double offset2 = p[10];

  double linear = offset + offsetGradient * t;
  double exponential = ampExp * exp( (offsetExp - t) / tcbo);
  double sinusoid  =                  amp * sin(    TMath::TwoPi() * (t-phase) / pow(omega,freqPower)); 
  double sinusoid2 = offset2 + amplitudeRatio * amp * sin(2 * TMath::TwoPi() * (t-phase) / pow(omega,freqPower)); 
  
  //return linear;
  return linear + (exponential * (sinusoid - sinusoid2));
  
}


//for the width fix the phase and let the offset float
double CBOWidthFunc(double* x, double* p){

  double t = x[0];
  double tMin = p[0];
  double tMax = p[1];
  double CBOPeriod = p[2]; // here the CBO period has already been multiplied by nOsc
  double phase = p[3];
  double nOsc = p[4];
  double amplitudeRatio = p[5];

  if(t < tMin or t > tMax){
    return 0;
  }
  
  //TO DO HERE - maybe nOScToFit is required?
  int segment = (t-tMin)/CBOPeriod;
  //std::cout << "segment = " << segment << "\tp[2*segment+5] = p[" << 2*segment+5 << "] = " << p[2*segment+5] << endl;
  double A =                p[3*segment+6]*sin(    TMath::TwoPi()*(t-phase)/( p[3*segment+7] * (CBOPeriod / nOsc)));
  double B = amplitudeRatio*p[3*segment+6]*sin(2 * TMath::TwoPi()*(t-phase)/( p[3*segment+7] * (CBOPeriod / nOsc)));
  return A - B + p[3*segment+8];

}

void plotsCBOWidthFit() {

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

  string outFolder = "CBOFits/Width/";
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
    string stationName = "station" + station;

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
    
    ///////////////////////////////////////
    //calculate the width
    //////////////////////////////////////

    double timeBinWidth = avgRadPosVsTime->GetBinWidth(1); // make it the same as the profile, in micro seconds
    double tOrig = timeBinWidth * 100.0; //start time, make it start exactly on a bin edge for exact matching (bins ~150ns wide)
    double tW = tOrig; //for modifying in loop
    double tEnd = 100.0;
    int nbins = int((tEnd - tW) / timeBinWidth); 
    tEnd = tOrig + (nbins * timeBinWidth); //modify tEnd so it - tOrig is a multiple of the bin width
    //double tEnd = tW + (width * 300.0);
    TH2F* h1clone2 = (TH2F*)h11_shifted->Clone();
    
    //can rebin the x axis (radial position) to remove statistical bumps
    int rebinFactor = 4;
    cout << "original profile bin width: " << timeBinWidth << "\n";

    //histograms to keep the calculated averages and widths 
    TH1D* hMean  = new TH1D("hmean", stationName.c_str(), nbins, tOrig, tEnd);
    TH1D* hWidth = new TH1D("hwidth", stationName.c_str(),nbins, tOrig, tEnd);

    hMean->GetYaxis()->SetTitle("Average Radial Position [mm]");
    hWidth->GetYaxis()->SetTitle("Standard Deviation [mm]");
    hMean->GetXaxis()->SetTitle("Time After Injection [us]");
    hWidth->GetXaxis()->SetTitle("Time After Injection [us]");

    //for gaussian fits
    TF1* f1 = new TF1("f1", "[0]*exp( -0.5 * (([1] - x) / [2])^2 )", -50, 50);
    f1->SetParameter(0,30);
    f1->SetParameter(1,15);
    f1->SetParameter(2,15);
    f1->SetParLimits(2,0,50);

    //we have 3 options: 
    //the first controls whether to use the fitted gaussian mean as the mean, or whether to use the TProfile mean
    //the second controls whether to use the fitted gaussian width as the width, or whether to use the stddeviation
    //the final one controls whether to recalculate the stddev and mean for the bins > 10% of the maximum
    bool useFittedMean = false;
    bool useFittedWidth = false;
    bool useReFit = true;
    double fraction = 0.0;

    cout<< "hmean bin width: " << hMean->GetBinWidth(1) << " tprof bin width: " << avgRadPosVsTime->GetBinWidth(1) << "\n";

    int iFit = 0;
    while (tW < tEnd){
      
      // get the bins corresponding to the desired width and project the TH2F
      int Bin1 = h1clone2->GetXaxis()->FindBin(tW);
      int Bin2 = h1clone2->GetXaxis()->FindBin(tW+timeBinWidth);
      TH1D* h1projY = h1clone2->ProjectionY("tmp", Bin1, Bin2);
      h1projY->Rebin(rebinFactor);

      h1projY->Draw();
      //c1->SaveAs(Form("%s/prefit_%s%s",outFolder.c_str(), stationName.c_str(), ext.c_str()));
      
      double maxY = h1projY->GetMaximum();
      int xOfMaxY = h1projY->GetXaxis()->GetBinCenter(h1projY->GetMaximumBin());

      double rangeMin = xOfMaxY - 40.0;
      double rangeMax = xOfMaxY + 40.0;
      f1->SetRange(rangeMin, rangeMax);
      h1projY->Fit("f1", "QR");
      
      //fill mean and stddev hists
      //TODO: get mean from TProfile, not fits
      int tBin = hMean->FindBin(tW);// + (timeBinWidth/2.0));
      int profBin = avgRadPosVsTime->FindBin(tW);// + (timeBinWidth/2.0));
      hMean->SetBinContent(tBin, avgRadPosVsTime->GetBinContent(profBin));
      hMean->SetBinError(tBin, avgRadPosVsTime->GetBinError(profBin));
      if (iFit < 5){
	cout << "for time: " << tW 
	     << " prof bin: " << profBin << " at: " << avgRadPosVsTime->GetBinCenter(profBin) << " has width: " << avgRadPosVsTime->GetBinWidth(profBin)
	     << " tbin: " << tBin << " at: " << hMean->GetBinCenter(tBin) << " has width " << hMean->GetBinWidth(tBin) << " profile has: " << hMean->GetBinContent(tBin) << "\n";
      }
      if (useFittedMean){
	hMean->SetBinContent(tBin, f1->GetParameter(1));
	hMean->SetBinError(tBin, f1->GetParError(1));
      }
      //cout << "for time: " << tW << " fitted mean: " << hMean->GetBinContent(tBin) << "\n";
      hWidth->SetBinContent(tBin, f1->GetParameter(2));
      hWidth->SetBinError(tBin, f1->GetParError(2));
      
      if (!useFittedWidth){
	// using the mean from the fit, calculate the stddev in all events over 10% max
	double fitLowR = 9999.0;
	double fitHighR = 0.0;
	double stddev = 0.0;
	int ntot = 0;
	for (int ibin(1); ibin < h1projY->GetXaxis()->GetNbins(); ibin++){
	  double mean = (useFittedMean)? f1->GetParameter(1) : hMean->GetBinContent(tBin);
	  double x = h1projY->GetBinCenter(ibin);
	  int n = h1projY->GetBinContent(ibin);
	  if (n > maxY*fraction){
	    stddev += n * pow(mean - x,2);
	    ntot+= n;
	    if (x < fitLowR) fitLowR = x;
	    if (x > fitHighR) fitHighR = x;
	  }
	}
	
        //don't know how to calculate error in std deviation, so cheat for now by at least putting some dependence on the number of entries
	double cheatError = 20.0 / sqrt(double(ntot)); 
	
	stddev = sqrt(stddev / (double)(ntot-1));      
	hWidth->SetBinContent(tBin, stddev);
	hWidth->SetBinError(tBin, cheatError);
      
	//refit and look at difference...
	f1->SetRange(fitLowR, fitHighR);
	h1projY->Fit("f1", "QR");
	double mean2 = (useFittedMean)? f1->GetParameter(1) : hMean->GetBinContent(tBin);
	
	fitLowR = 9999.9;
	fitHighR = 0.0;
	double stddev2 = 0.0;
	ntot = 0;
	for (int ibin(1); ibin < h1projY->GetXaxis()->GetNbins(); ibin++){
	  double x = h1projY->GetBinCenter(ibin);
	  int n = h1projY->GetBinContent(ibin);
	  if (n > maxY*fraction){
	    stddev2 += n * pow(mean2 - x,2);
	    ntot+= n;
	    if (x < fitLowR) fitLowR = x;
	    if (x > fitHighR) fitHighR = x;
	  }
	}
	stddev2 = sqrt(stddev2 / (double)(ntot-1));
	//cout << "old mean: " << mean1 << " new mean: " << mean2 << " ""\n";
	//cout << "old stddev: " << stddev << " new stddev: " << stddev2 << " ""\n";
	
	if (useReFit){
	  if (useFittedMean){
	    hMean->SetBinContent(tBin, f1->GetParameter(1));
	    hMean->SetBinError(tBin, f1->GetParError(1));
	  }
	  if (useFittedWidth){
	    hWidth->SetBinContent(tBin, f1->GetParameter(2));
	    hWidth->SetBinError(tBin, f1->GetParError(2));
	  }
	  else{
	    hWidth->SetBinContent(tBin, stddev2);
	    hWidth->SetBinError(tBin, cheatError);
	  }
	}
      }

      //// to save output into dedicated Fit folder to make video...
      //std::stringstream ss;
      //ss << std::setw(4) << std::setfill('0') << iFit;
      //std::string s = ss.str();
      //
      //std::stringstream ss2;
      //ss2 << std::setprecision(4) << "time: " << tW << "-" << tW+width << " us";
      //std::string s2 = ss2.str();
      //
      //TPaveText *pt = new TPaveText(0.15, 0.75, 0.30, 0.85, "NDC");
      //pt->AddText(s2.c_str());
      //pt->SetBorderSize(0);
      //pt->SetFillColor(0);
      //
      //h1projY->SetTitle(stationName.c_str());
      //h1projY->GetYaxis()->SetRangeUser(0,200);
      //h1projY->GetXaxis()->SetTitleOffset(0.9);
      //h1projY->Draw("");
      //pt->Draw("SAME");
      //c1->SaveAs(Form("%s/Fit/fit%s_%s%s",outFolder.c_str(), s.c_str(), stationName.c_str(),ext.c_str()));

      tW += timeBinWidth;
      iFit++;
    }

    ////as a sanity check, when useFitedMean = false, hMean should be exactly the profile plot
    //avgRadPosVsTime->GetXaxis()->SetRangeUser(tOrig, tEnd);
    //avgRadPosVsTime->Draw("Hist P");
    //hMean->Draw("SAME");
    //c1->SaveAs(Form("%s/tmp_%s%s",outFolder.c_str(), stationName.c_str(), ext.c_str()));
    //return;

    hMean->Draw();
    c1->SaveAs(Form("%s/mean_%s%s",outFolder.c_str(), stationName.c_str(), ext.c_str()));
    
    hMean->SetMarkerStyle(33);
    hMean->Draw("HIST P");
    c1->SaveAs(Form("%s/mean_noError_%s%s",outFolder.c_str(), stationName.c_str(), ext.c_str()));

    hWidth->Draw();
    c1->SaveAs(Form("%s/stddev_%s%s",outFolder.c_str(), stationName.c_str(), ext.c_str()));

    hWidth->SetMarkerStyle(33);
    hWidth->Draw("HIST P");
    c1->SaveAs(Form("%s/stddev_noError_%s%s",outFolder.c_str(), stationName.c_str(), ext.c_str()));

    ////////////////////////////////
    // Amplitude of CBO Width - intial fit
    ////////////////////////////////
    double nOscToFit = 4.;
    double CBOPeriod =  nOscToFit * 2.43168; //us - the factor of 4 determines how many oscillation to treat as constant amplitude for the initial fit
    double tMin = 20.0; //65.0 - shiftTime;
    double tMax = tEnd; // 250 - shiftTime;
    int nSegments = (tMax-tMin)/CBOPeriod;
    tMax = nSegments*CBOPeriod + tMin;

    TF1* CBOWidth = new TF1("CBOWidth", CBOWidthFunc, tMin, tMax, 3*nSegments + 6); // 0: Ti, 1: Tend, 2: fit period, 3: global phase , 4: nOsc 5: ampRatio, 6: amp, 7:local freq 8:offset

    double phaseW = (stationName == "station12")? 1.3 : 2.0;
    CBOWidth->FixParameter(0,tMin);
    CBOWidth->FixParameter(1,tMax);
    CBOWidth->FixParameter(2,CBOPeriod);
    CBOWidth->SetParameter(3,phaseW);//global phase
    CBOWidth->FixParameter(4,nOscToFit); 
    CBOWidth->SetParameter(5,0.5);//amplitude ratio
    CBOWidth->SetParLimits(5,0.0, 2.0);//amplitude ratio
    
    for(int par = 6; par < CBOWidth->GetNpar(); par+=3){
      CBOWidth->SetParameter(par,5.0); //amplitude
      CBOWidth->SetParLimits(par,0.0, 8.0); //amplitude
      CBOWidth->SetParameter(par+1,1.0); // local freq
      CBOWidth->SetParLimits(par+1,0.9, 1.1); // local freq
      CBOWidth->SetParameter(par+2,15.0); // floating offset
      CBOWidth->SetParLimits(par+2,10.0, 30.0); // floating offset
    }
    CBOWidth->SetNpx(2000);

    //Plot the width at different time intervals
    hWidth->GetYaxis()->SetRangeUser(hWidth->GetMinimum(), hWidth->GetMaximum());
    hWidth->SetMarkerStyle(33);
    hWidth->Draw();
    hWidth->Fit(CBOWidth, "RL");
    cout << "chi square: " << CBOWidth->GetChisquare() << " ndf: " << CBOWidth->GetNDF() << " giving: " << CBOWidth->GetChisquare() / double(CBOWidth->GetNDF()) << "\n";
    CBOWidth->Draw("SAME");
    c1->SaveAs(Form("%s/CBOWidthFit_%s%s",outFolder.c_str(), stationName.c_str(), ext.c_str()));

    TH1* widthClone = (TH1*) hWidth->Clone();
    widthClone->GetXaxis()->SetRangeUser(tOrig, tOrig+15.0);
    widthClone->Draw(); //HIST P for no errors
    CBOWidth->Draw("SAME");
    c1->SaveAs(Form("%s/CBOWidthFit_early_%s%s",outFolder.c_str(), stationName.c_str(), ext.c_str()));

    widthClone->GetXaxis()->SetRangeUser(60.,75.);
    widthClone->Draw();
    CBOWidth->Draw("SAME");
    c1->SaveAs(Form("%s/CBOWidthFit_mid_%s%s",outFolder.c_str(), stationName.c_str(), ext.c_str()));

    widthClone->GetXaxis()->SetRangeUser(90.,105.);
    widthClone->Draw();
    CBOWidth->Draw("SAME");
    c1->SaveAs(Form("%s/CBOWidthFit_midlate_%s%s",outFolder.c_str(), stationName.c_str(), ext.c_str()));

    widthClone->GetXaxis()->SetRangeUser(120.,135.);
    widthClone->Draw();
    CBOWidth->Draw("SAME");
    c1->SaveAs(Form("%s/CBOWidthFit_late_%s%s",outFolder.c_str(), stationName.c_str(), ext.c_str()));

    // plot the parameters as a function of time
    TGraphErrors* tgCBOWidth = new TGraphErrors();
    TGraphErrors* tgCBOWFreq = new TGraphErrors();
    TGraphErrors* tgCBOOffset = new TGraphErrors();
    double averageOffset = 0.0;
    double tgCBOWidthInitial = 0.0;
    int counter = 0;
    for(int i = 0; i < 3*nSegments; i+=3){
      int nPts = tgCBOWidth->GetN();
      tgCBOWidth->SetPoint(nPts, tMin+CBOPeriod*nPts, CBOWidth->GetParameter(6+i));
      tgCBOWidth->SetPointError(nPts, 0, CBOWidth->GetParError(6+i));
      tgCBOWFreq->SetPoint(nPts, tMin+CBOPeriod*nPts, CBOWidth->GetParameter(7+i));
      tgCBOWFreq->SetPointError(nPts, 0, CBOWidth->GetParError(7+i));
      tgCBOOffset->SetPoint(nPts, tMin+CBOPeriod*nPts, CBOWidth->GetParameter(8+i));
      tgCBOOffset->SetPointError(nPts, 0, CBOWidth->GetParError(8+i));
      
      //for the fits later
      averageOffset += CBOWidth->GetParameter(8+i);
      if (i == 0) tgCBOWidthInitial = CBOWidth->GetParameter(6+i);

      counter++;
      cout << nPts << " SETTING Width vs TIME: " << " p[ " << 6+i << "]: " << CBOWidth->GetParameter(6+i) << " CBO T: " << CBOPeriod  
	   << " t: " << tMin+CBOPeriod*nPts << " from i: " << i << " npar: " << CBOWidth->GetNpar() <<"\n";
    }

    averageOffset /= double(counter);

    tgCBOWidth->SetMarkerColor(2);
    tgCBOWidth->SetLineColor(2);
    tgCBOWidth->SetMarkerStyle(2);
    //tgCBOWidth->GetXaxis()->SetRangeUser(,100.);
    tgCBOWidth->GetYaxis()->SetRangeUser(0.,2.);
    tgCBOWidth->SetTitle(";Time [#mus];Radial CBO width amplitude [mm]");
    tgCBOWidth->Draw("AP");
    c1->SaveAs(Form("%s/CBOWidth_time_%s%s",outFolder.c_str(), stationName.c_str(), ext.c_str()));

    tgCBOWFreq->SetMarkerColor(2);
    tgCBOWFreq->SetLineColor(2);
    tgCBOWFreq->SetMarkerStyle(2);
    tgCBOWFreq->GetYaxis()->SetRangeUser(0.97,1.03);
    tgCBOWFreq->SetTitle(";Time [#mus];local frequency for width fits");
    tgCBOWFreq->Draw("AP");
    c1->SaveAs(Form("%s/CBOWFreq_time_%s%s",outFolder.c_str(), stationName.c_str(), ext.c_str()));

    tgCBOOffset->SetMarkerColor(2);
    tgCBOOffset->SetLineColor(2);
    tgCBOOffset->SetMarkerStyle(2);
    tgCBOOffset->GetYaxis()->SetRangeUser(14.,25.);
    tgCBOOffset->SetTitle(";Time [us]; Radial CBO width [mm]");
    tgCBOOffset->Draw("AP");
    c1->SaveAs(Form("%s/CBOWidthOffset_time_%s%s",outFolder.c_str(), stationName.c_str(), ext.c_str()));

    cout << stationName << " original phase for width: " << phaseW << " after fit: " << CBOWidth->GetParameter(3) << "\n";
    //

    //plot mean and stddev on teh same canvas
    TCanvas* c2 = new TCanvas("c2","",800,600);
    c2->SetBottomMargin(0.13);
    //gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    Float_t small = 1e-5;
    
    c2->Divide(1,2,small,small);
    c2->cd(1);
    gPad->SetBottomMargin(small);
    //gPad->SetLeftMargin(1.3);
    
    //tEnd = 40.0;
    double ymin_width = hWidth->GetMinimum();
    double ymax_width = hWidth->GetMaximum();

    hWidth->SetLineWidth(2);
    hWidth->SetLineColor(kBlue);
    hWidth->SetMarkerColor(kBlack);
    hWidth->GetXaxis()->SetRangeUser(tMin, tMax);
    hWidth->GetYaxis()->SetTitle("#sigma(x) [mm]");
    hWidth->GetYaxis()->CenterTitle();
    hWidth->GetYaxis()->SetRangeUser(ymin_width,ymax_width);
    hWidth->GetYaxis()->SetTitleOffset(0.6);
    hWidth->GetYaxis()->SetTitleSize(0.08);
    hWidth->GetYaxis()->SetLabelSize(0.06);
    hWidth->Draw("HIST P");
    CBOWidth->Draw("SAME");

    vector<double> CBOtimes;
    double startTime = 0.0;
    //for (int i(0); i < 10; i++){
    while (startTime < tMax){
      if (startTime < tMin) {
	startTime += (CBOPeriod / nOscToFit);
	continue;
      }
      CBOtimes.push_back(startTime);
      startTime += (CBOPeriod / nOscToFit);
    }
    
    for (auto t : CBOtimes){
      TLine* tl = new TLine(t, ymin_width, t, ymax_width);
      tl->SetLineStyle(2);
      tl->SetLineColor(3);
      tl->Draw("SAME");
    }

    c2->cd(2);
    gPad->SetTopMargin(small);
    gPad->SetBottomMargin(0.15); // change this number
    //gPad->SetLeftMargin(1.3);
    gPad->SetTickx();
    
    double ymin_mean = hMean->GetMinimum();
    double ymax_mean = hMean->GetMaximum();

    hMean->SetLineWidth(2);
    hMean->SetLineColor(kRed);
    hMean->SetMarkerColor(kBlack);
    hMean->SetTitle("");
    hMean->GetYaxis()->SetTitle("Average radius [mm]");
    hMean->GetXaxis()->SetTitle("Time since injection[#mus]");
    hMean->GetXaxis()->CenterTitle();
    hMean->GetXaxis()->SetRangeUser(tMin,tMax);
    hMean->GetYaxis()->SetRangeUser(ymin_mean,ymax_mean);
    hMean->GetYaxis()->CenterTitle();
    hMean->GetXaxis()->SetTitleSize(0.08);
    hMean->GetYaxis()->SetTitleSize(0.08);
    hMean->GetXaxis()->SetTitleOffset(0.7);
    hMean->GetYaxis()->SetTitleOffset(0.6);
    hMean->GetXaxis()->SetLabelSize(0.06);
    hMean->GetYaxis()->SetLabelSize(0.06);
    hMean->Draw();

    for (auto t : CBOtimes){
      TLine* tl = new TLine(t, ymin_mean, t, ymax_mean);
      tl->SetLineStyle(2);
      tl->SetLineColor(3);
      tl->Draw("SAME");
    }

    c2->SaveAs(Form("%s/meanStddev_%s%s",outFolder.c_str(), stationName.c_str(), ext.c_str()));
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    //now perform the final fit, where we use the output of the old fit, but allow the parameters to float
    //////////////////////////////////////////////////////////////////////////////////////////////

    TF1* CBOWidthFinal = new TF1("CBOWidthFinal", CBOWidthFinalFunc, tMin, tMax, 11);
    //linear
    double linearOffset = averageOffset; //CBOWidth->GetParameter(4);
    double linearGrad = 0.0;
    //exp
    double expA = 1.0;
    double expOffset = 1.0;//exp1->GetParameter(1);
    double tCBO = 80.0;//exp1->GetParameter(2);
    //sin
    double sinA = tgCBOWidthInitial;;
    double phase = CBOWidth->GetParameter(3);
    //double freq = 1000.0/(y * CBOPeriod); // in kHz
    double omega = CBOWidth->GetParameter(2)/ 4.0; //initial frequency for now
    double freqPower = 1.0;
    double ampRatio = CBOWidth->GetParameter(5); 

    double offset2 = 0.0;

    bool fixall = false;

    CBOWidthFinal->SetParName(0, "linear offs");
    if (fixall) CBOWidthFinal->FixParameter(0, linearOffset);
    else {
      CBOWidthFinal->SetParameter(0, linearOffset);
      CBOWidthFinal->SetParLimits(0, 0.95*linearOffset, 1.05*linearOffset);
    }
    
    CBOWidthFinal->SetParName(1, "linear grad");
    if (fixall) CBOWidthFinal->FixParameter(1, linearGrad);
    else {
      CBOWidthFinal->SetParameter(1, linearGrad);
      CBOWidthFinal->SetParLimits(1, -0.9, 1.1);
    }

    CBOWidthFinal->SetParName(2, "exp amp    ");
    CBOWidthFinal->FixParameter(2, expA);
    
    CBOWidthFinal->SetParName(3, "exp offset ");
    if (fixall) CBOWidthFinal->FixParameter(3, expOffset);
    else {
      CBOWidthFinal->SetParameter(3, expOffset);
      CBOWidthFinal->SetParLimits(3, 0.5*expOffset, 1.5*expOffset);
    }
    
    CBOWidthFinal->SetParName(4, "tau CBO    ");
    if (fixall) CBOWidthFinal->FixParameter(4, tCBO);
    else {
      CBOWidthFinal->SetParameter(4, tCBO);
      CBOWidthFinal->SetParLimits(4, 0.5*tCBO, 2.5*tCBO);
    }

    CBOWidthFinal->SetParName(5, "sin amp    ");
    if (fixall) CBOWidthFinal->FixParameter(5, sinA);
    else {
      CBOWidthFinal->SetParameter(5, sinA);
      CBOWidthFinal->SetParLimits(5, 0.5*sinA, 1.5*sinA);
    }
    
    CBOWidthFinal->SetParName(6, "sin phase  ");
    if (fixall) CBOWidthFinal->FixParameter(6, phase);
    else {
      CBOWidthFinal->SetParameter(6, phase);
      CBOWidthFinal->SetParLimits(6, 0.1*phase, 3.0*phase);
    }
    
    CBOWidthFinal->SetParName(7, "omega      ");
    if (fixall) CBOWidthFinal->FixParameter(7, omega);
    else {
      CBOWidthFinal->SetParameter(7, omega);
      CBOWidthFinal->SetParLimits(7, 0.99*omega, 1.01*omega);
    }

    CBOWidthFinal->SetParName(8, "omega power");
    if (fixall)CBOWidthFinal->FixParameter(8, freqPower);
    else {
      CBOWidthFinal->SetParameter(8, freqPower);
      CBOWidthFinal->SetParLimits(8, 0.95, 1.05);
    }

    CBOWidthFinal->SetParName(9, "amp ratio ");
    if (fixall)CBOWidthFinal->FixParameter(9, ampRatio);
    else {
      CBOWidthFinal->SetParameter(9, ampRatio);
      CBOWidthFinal->SetParLimits(9, 0.01, 2.0);
    }

    CBOWidthFinal->SetParName(10, "offsetAlt");
    if (fixall)CBOWidthFinal->FixParameter(10, offset2);
    else {
      CBOWidthFinal->SetParameter(10, offset2);
      CBOWidthFinal->SetParLimits(10, -1.0, 1.0);
    }

    CBOWidthFinal->SetNpx(10000);

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

    double offsetFactor = 1.03;

    pad1->cd();
    pad1->SetTopMargin(0.02);
    pad1->SetLeftMargin(0.08);
    pad1->SetRightMargin(0.02);
    hWidth->GetXaxis()->SetRangeUser(tMin, tMax);
    hWidth->GetYaxis()->SetRangeUser(18.0, hWidth->GetMaximum()*1.08);
    hWidth->GetYaxis()->SetLabelSize(0.04);
    hWidth->SetTitle("");
    double range = tMax - tMin;
    TPaveText* pt = new TPaveText(tMin + (range/2.0) - 10.0, hWidth->GetMaximum()*offsetFactor-2.0, tMin + (range/2.0) + 10.0, hWidth->GetMaximum()*offsetFactor-1.2);
    pt->AddText(("Station " + station).c_str());
    pt->SetTextFont(gStyle->GetTextFont());
    pt->SetFillColor(10);
    hWidth->GetYaxis()->SetTitle("Radial Beam Width [mm]");
    hWidth->GetXaxis()->SetTitleSize(0.05);
    hWidth->GetXaxis()->SetTitleOffset(0.8);
    hWidth->GetYaxis()->SetTitleSize(0.05);
    hWidth->GetYaxis()->SetTitleOffset(0.7);
    hWidth->GetXaxis()->CenterTitle();
    hWidth->GetYaxis()->CenterTitle();
    hWidth->Fit(CBOWidthFinal,"R");
    hWidth->Draw();
    pt->Draw("SAME");

    pad2->cd();
    pad2->SetLeftMargin(0.04);
    pad2->SetRightMargin(0.0);
    pad2->SetTopMargin(0.01);
    TProfile* hWidth2 = (TProfile*)hWidth->Clone();
    hWidth2->GetXaxis()->SetRangeUser(tMin, tMin+20.0);
    stringstream ss2;
    ss2 << std::setprecision(4) << tMin << " - " << tMin+20.0 << " [#mus]";
    TPaveText* pt2 = new TPaveText(tMin+2., hWidth->GetMaximum()*offsetFactor-2.0, tMin+18., hWidth->GetMaximum()*offsetFactor-1.0);
    pt2->AddText(ss2.str().c_str());
    pt2->SetFillColor(10);
    hWidth2->GetXaxis()->SetTitle("");
    hWidth2->Draw();
    pt2->Draw("SAME");

    pad3->cd();
    pad3->SetLeftMargin(0.0);
    pad3->SetRightMargin(0.0);
    pad3->SetTopMargin(0.01);
    TProfile* hWidth3 = (TProfile*)hWidth->Clone();
    hWidth3->GetXaxis()->SetRangeUser((tMin + tMax)/2.0 , (tMin + tMax)/2.0 + 20.0);
    stringstream ss3;
    ss3 << std::setprecision(4) << (tMin + tMax)/2.0 << " - " << (tMin + tMax) /2.0 +20.0 << " [#mus]";
    TPaveText* pt3 = new TPaveText( (tMin+tMax)/2. +2., hWidth->GetMaximum()*offsetFactor-2.0, (tMin+tMax)/2. +18., hWidth->GetMaximum()*offsetFactor-1.0);
    pt3->AddText(ss3.str().c_str());
    pt3->SetFillColor(10);
    hWidth3->GetXaxis()->SetTitle("");
    hWidth3->Draw();
    pt3->Draw("SAME");

    pad4->cd();
    pad4->SetLeftMargin(0.0);
    pad4->SetRightMargin(0.02);
    pad4->SetTopMargin(0.01);
    TProfile* hWidth4 = (TProfile*)hWidth->Clone();
    hWidth4->GetXaxis()->SetRangeUser( tMax - 20.0 , tMax);
    stringstream ss4;
    ss4 << std::setprecision(4) << tMax - 20.0 << " - " << tMax << " [#mus]";
    TPaveText* pt4 = new TPaveText(tMax-18., hWidth->GetMaximum()*offsetFactor-2.0, tMax-2., hWidth->GetMaximum()*offsetFactor-1.0);
    pt4->AddText(ss4.str().c_str());
    pt4->SetFillColor(10);
    hWidth4->GetXaxis()->SetTitle("");
    hWidth4->Draw();
    pt4->Draw("SAME");

    c4->SaveAs(Form("%s/FinalCBOWidthFit_%s%s",outFolder.c_str(), stationName.c_str(), ext.c_str()));

  }
  return;

}
