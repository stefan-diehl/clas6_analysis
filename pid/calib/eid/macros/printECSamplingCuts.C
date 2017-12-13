{

  TFile *inputFile = TFile::Open("data.root");

  TH2I *ecSampling[6];
  TF1  *mu[6];
  TF1  *sigma[6];
  TF1  *upper[6];
  TF1  *lower[6];

  double NSIGMA = 1.50; 

  TCanvas *compareCanvas = new TCanvas("compareCanvas","",800,1200);  
  compareCanvas->Divide(2,3);

  for(int sector=0; sector<6; sector++){
    ecSampling[sector] = (TH2I*) inputFile->Get(Form("ecSamplingFraction_%d",sector));
    mu[sector]         = (TF1*)  inputFile->Get(Form("ecPol3Mu%d",   sector));
    sigma[sector]      = (TF1*)  inputFile->Get(Form("ecPol3Sigma%d",sector));


    double xMin = mu[sector]->GetXmin();
    double xMax = mu[sector]->GetXmax();
    upper[sector]      = new TF1(Form("upper%d",sector),"pol3",xMin,xMax);
    lower[sector]      = new TF1(Form("lower%d",sector),"pol3",xMin,xMax);

    upper[sector]->SetParameter(0, mu[sector]->GetParameter(0)+NSIGMA*sigma[sector]->GetParameter(0));
    upper[sector]->SetParameter(1, mu[sector]->GetParameter(1)+NSIGMA*sigma[sector]->GetParameter(1));
    upper[sector]->SetParameter(2, mu[sector]->GetParameter(2)+NSIGMA*sigma[sector]->GetParameter(2));
    upper[sector]->SetParameter(3, mu[sector]->GetParameter(3)+NSIGMA*sigma[sector]->GetParameter(3));

    lower[sector]->SetParameter(0, mu[sector]->GetParameter(0)-NSIGMA*sigma[sector]->GetParameter(0));
    lower[sector]->SetParameter(1, mu[sector]->GetParameter(1)-NSIGMA*sigma[sector]->GetParameter(1));
    lower[sector]->SetParameter(2, mu[sector]->GetParameter(2)-NSIGMA*sigma[sector]->GetParameter(2));
    lower[sector]->SetParameter(3, mu[sector]->GetParameter(3)-NSIGMA*sigma[sector]->GetParameter(3));

    cout << "[Pol3 Params] UPPER, Sector=" << sector << " A=" << upper[sector]->GetParameter(3) << 
      " B=" << upper[sector]->GetParameter(2) << 
      " C=" << upper[sector]->GetParameter(1) << 
      " D=" << upper[sector]->GetParameter(0) << endl;


    compareCanvas      ->cd(sector+1);
    ecSampling[sector] ->Draw("colz"); 
    mu[sector]         ->Draw("same");
    upper[sector]      ->Draw("same");
    lower[sector]      ->Draw("same");
  }


}
