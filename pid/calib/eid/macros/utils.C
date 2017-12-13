void loadSlices(TFile *inputFile, TH1D *slices[numberSectors][numberSlices], const int numberSectors, const int numberSlices, string histoTitle){
  for(int sector=0; sector<numberSectors; sector++){
    for (int slice=0; slice<numberSlices; slice++){
      string title = Form("%s_%d_sector_%d",histoTitle.c_str(),slice,sector);
      slices[sector][slice] = (TH1D*) inputFile->Get(title.c_str()); 

      cout << " >>> Loaded " << title << " with N=" << slices[sector][slice]->GetEntries() << endl;
    }
  }


}


void printSlices(TH1D *slices[numberSectors][numberSlices], const int numberSectors, const int numberSlices, int histoType){

  // Histogram Type 
  // [0] -> EC Sampling Fraction Slices
  // [1] -> CC Segment Slices 

  TCanvas *sliceCanvas = new TCanvas("sliceCanvas","",800,800);
  TLatex title, xTitle, yTitle; 
  title.SetNDC();
  title.SetTextSize(0.04);
  title.SetTextFont(42);

  xTitle.SetNDC();
  xTitle.SetTextSize(0.04);
  xTitle.SetTextFont(42);

  yTitle.SetNDC();
  yTitle.SetTextSize(0.04);
  yTitle.SetTextAngle(90);
  yTitle.SetTextFont(42);

  string xTitleString, yTitleString, titleString; 
  double xTitleXPos, xTitleYPos, yTitleXPos, yTitleYPos, 
    titleXPos, titleYPos;

  xTitleXPos = 0.68; xTitleYPos = 0.05; 
  yTitleXPos = 0.05; yTitleYPos = 0.65; 
  titleXPos  = 0.48;  titleYPos = 0.92;

  for(int sector=0; sector<numberSectors; sector++){
    for (int slice=0; slice<numberSlices; slice++){
      slices[sector][slice]->SetFillColorAlpha(62, 0.4);
      slices[sector][slice]->Draw();      

      gPad->SetGridx();
      gPad->SetGridy();
      gPad->SetMargin(0.15, 0.15, 0.15, 0.15);

      if (histoType == 0){
	xTitleString = "etot/p";
	yTitleString = "Counts";
	titleString  = Form("Sector:%d Slice:%d",sector,slice);
      }

      if (histoType == 1){
	xTitleString = "#theta_{CC}";
	yTitleString = "Counts";
	titleString  = Form("Sector:%d Slice:%d",sector,slice);
      }

      title .DrawLatex(titleXPos,   titleYPos,  titleString.c_str());
      xTitle.DrawLatex(xTitleXPos, xTitleYPos, xTitleString.c_str());
      yTitle.DrawLatex(yTitleXPos, yTitleYPos, yTitleString.c_str());

      if (histoType == 0){ sliceCanvas->SaveAs(Form("../img/ecSamplingSlice%dSector%d.png",slice,sector)); }
      if (histoType == 1){ sliceCanvas->SaveAs(Form("../img/ccThetaCorrelationSlice%dSector%d.png",slice,sector)); }
    }
  }


}
