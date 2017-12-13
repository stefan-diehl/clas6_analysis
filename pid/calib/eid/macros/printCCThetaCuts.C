{

  TFile *inputFile = TFile::Open("../e16.root");

  TH2I         *cc[6];
  TGraphErrors *mu[6];
  TGraphErrors *sigma[6];

  double NSIGMA = 3.00; 

  TCanvas *compareCanvas = new TCanvas("compareCanvas","",800,1200);  
  compareCanvas->Divide(2,3);

  for(int sector=0; sector<6; sector++){
    cc[sector]    = (TH2I*) inputFile->Get(Form("ccThetaCorrelation_%d",sector));
    mu[sector]    = (TGraphErrors*) inputFile->Get(Form("ccMuGraph%d",sector));
    sigma[sector] = (TGraphErrors*) inputFile->Get(Form("ccSigmaGraph%d",sector));

    compareCanvas->cd(sector+1);
    gPad->SetLogz();
    cc[sector]   ->Draw("colz"); 
    mu[sector]   ->Draw("same");
  }

  

}
