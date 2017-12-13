{
  gStyle->SetPalette(62);

  TFile *inputFile      = TFile::Open("KinematicImpactMC.root");
  //  TFile *inputFile      = TFile::Open("KinematicImpact.root");
  string imagePath      = "/volatile/clas12/dmriser/plots/inclusive/eidCuts/";
  string configName     = "KeppelRad";

  TH1F *all = (TH1F*) inputFile->Get("all");
  vector<TH1F*> pass;
  vector<TH1F*> fail;

  // Load from file 
  TList *tableOfContents = inputFile->GetListOfKeys();

  TIter next(tableOfContents);

  while(TObject *objectFromFile = next()){    
    TString currentObjectName  = objectFromFile->GetName();

    if (currentObjectName.Contains("pCut"))      { pass.push_back((TH1F*) inputFile->Get(currentObjectName)); }
    else if (currentObjectName.Contains("fCut")) { fail.push_back((TH1F*) inputFile->Get(currentObjectName)); }


    cout << "[LoadingHistograms] Found " << currentObjectName << endl;
  }

  double xBins = (double) pass[0]->GetXaxis()->GetNbins();
  TCanvas *effCanvas = new TCanvas("effCanvas","",1200,400); 
  TLatex *title      = new TLatex(); 
  TLatex *xLabel     = new TLatex();
  TLatex *yLabel     = new TLatex();
  TLatex *stats      = new TLatex();

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
  gPad->SetLogz(); 

  for (int i=0; i<pass.size(); ++i){
    effCanvas->Divide(3,1);

    effCanvas->cd(1);
    all->Draw("colz");
    title   ->DrawLatex(0.28, 0.87," All Negatives "); 
    yLabel  ->DrawLatex(0.05, 0.63, "Q^{2} [GeV^{2}/c^{2}]"); 

    effCanvas->cd(2);
    pass[i] ->Draw("colz"); 
    title   ->DrawLatex(0.28, 0.87, Form("Pass %s",pass[i]->GetTitle())); 


    effCanvas->cd(3);
    fail[i] ->Draw("colz"); 
    title   ->DrawLatex(0.28, 0.87, Form("Fail %s",fail[i]->GetTitle())); 
    xLabel  ->DrawLatex(0.74, 0.07, "W [GeV/c^{2}]"); 

    effCanvas->cd();
    stats->DrawLatex(0.05, 0.02,Form("Left: Pass %.2f, Right: Fail %.2f",pass[i]->GetEntries()/(pass[i]->GetEntries() + fail[i]->GetEntries()),fail[i]->GetEntries()/(pass[i]->GetEntries() + fail[i]->GetEntries())));
  
    effCanvas->Print(Form("%s%sCut%d.png",imagePath.c_str(),configName.c_str(),i)); 
    effCanvas->Clear();
  }
}
