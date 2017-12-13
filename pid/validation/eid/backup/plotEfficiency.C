{
  TFile *inputFile      = TFile::Open("Efficiency.root");
  TFile *inputFileLoose = TFile::Open("EfficiencyLoose.root");
  TFile *inputFileTight = TFile::Open("EfficiencyTight.root");

  vector<TH1F*> histo;
  vector<TH1F*> histoLoose;
  vector<TH1F*> histoTight;

  // Load from file 
  TList *tableOfContents      = inputFile->GetListOfKeys();
  TList *tableOfContentsLoose = inputFile->GetListOfKeys();
  TList *tableOfContentsTight = inputFile->GetListOfKeys();

  TIter next(tableOfContents);
  TIter nextLoose(tableOfContentsLoose);
  TIter nextTight(tableOfContentsTight);

  while(TObject *objectFromFile = next()){    
    TString currentObjectName  = objectFromFile->GetName();
    histo.push_back((TH1F*) inputFile->Get(currentObjectName));
    cout << "[LoadingHistograms] Found " << currentObjectName << endl;
  }
  while(TObject *objectFromFile = nextLoose()){    
    TString currentObjectName  = objectFromFile->GetName();
    histoLoose.push_back((TH1F*) inputFileLoose->Get(currentObjectName));
    cout << "[LoadingHistograms] Found " << currentObjectName << endl;
  }
  while(TObject *objectFromFile = nextTight()){    
    TString currentObjectName  = objectFromFile->GetName();
    histoTight.push_back((TH1F*) inputFileTight->Get(currentObjectName));
    cout << "[LoadingHistograms] Found " << currentObjectName << endl;
  }


  double xBins = (double) histo[0]->GetXaxis()->GetNbins();
  TCanvas *effCanvas = new TCanvas("effCanvas","",800,800); 
  TLatex *title = new TLatex(); 
  TLatex *xLabel = new TLatex();
  TLatex *yLabel = new TLatex();
  TLatex *stats  = new TLatex();

  title->SetNDC();
  title->SetTextSize(0.04);
  title->SetTextFont(42);

  stats->SetNDC();
  stats->SetTextSize(0.03);
  stats->SetTextFont(42);

  xLabel->SetNDC();
  xLabel->SetTextSize(0.03);
  xLabel->SetTextFont(42);

  yLabel->SetNDC();
  yLabel->SetTextSize(0.03);
  yLabel->SetTextFont(42);
  yLabel->SetTextAngle(90.0);
 
  gPad->SetGridx();
  gPad->SetGridy();

  for (int i=0; i<histo.size(); ++i){
    //    histo[i]->SetFillColorAlpha(65, 0.3);
    histo[i]     ->SetMinimum(0.0); 
    histo[i]     ->SetMaximum(1.0); 
    histo[i]     ->SetLineColor(55); 
    histoLoose[i]->SetLineColor(75); 
    histoTight[i]->SetLineColor(95); 
    histo[i]     ->Draw("l"); 
    histoLoose[i]->Draw("lsame"); 
    histoTight[i]->Draw("lsame"); 
    //    histo[i]->GetXaxis()->SetLabelSize(0.025);
    //    histo[i]->GetYaxis()->SetLabelSize(0.025);
    title   ->DrawLatex(0.28, 0.87, histo[i]->GetTitle()); 
    xLabel  ->DrawLatex(0.74, 0.07, "P [GeV/c]"); 
    yLabel  ->DrawLatex(0.05, 0.63, "Signal Efficiency"); 
    stats   ->SetTextColor(75);
    stats   ->DrawLatex(0.02, 0.08, Form("Loose  : %.2f",histoLoose[i]->GetEntries()/xBins));
    stats   ->SetTextColor(55);
    stats   ->DrawLatex(0.02, 0.05, Form("Nominal: %.2f",histo[i]->GetEntries()/xBins));
    stats   ->SetTextColor(95);
    stats   ->DrawLatex(0.02, 0.02, Form("Tight  : %.2f",histoTight[i]->GetEntries()/xBins));
    effCanvas->Print(Form("img/eff/cut%d.png",i)); 

  }

}
