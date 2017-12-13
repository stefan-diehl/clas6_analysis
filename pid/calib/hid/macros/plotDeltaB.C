{

  std::string inputFile("/volatile/clas12/dmriser/rootFiles/pid/meson/return_meson1.root");
  TFile *file = TFile::Open(inputFile.c_str());

  std::string mesonTitle = "kp";
  std::string particleLatex = "K^{+}";
  TH2F *deltaBeta[6];

  for(int i=1; i<7; i++){
    deltaBeta[i-1] = (TH2F*) file->Get(Form("MesonHistograms/h2_p_dbeta_%d_%s",i,mesonTitle.c_str()));
  }

  TLine *upper = new TLine(0.0,0.02,4.5,0.02); 
  upper->SetLineStyle(8); 
  upper->SetLineWidth(2); 

  TLine *lower = new TLine(0.0,-0.02,4.5,-0.02); 
  lower->SetLineStyle(8); 
  lower->SetLineWidth(2); 

  TCanvas *can = new TCanvas("can","",800,1100);
  can->Divide(2,3);

  TLatex tit, xtit, ytit;
  tit.SetNDC();
  tit.SetTextFont(102);
  tit.SetTextSize(0.03);
  xtit.SetNDC();
  xtit.SetTextFont(102);
  xtit.SetTextSize(0.05);
  ytit.SetNDC();
  ytit.SetTextFont(102);
  ytit.SetTextSize(0.05);
  ytit.SetTextAngle(90.0);

  double xtitx = 0.0;
  double xtity = 0.0;

  for(int s=0; s<6; s++){
    can->cd(s+1);
 
    //    TLine *protonMass = new TLine(0.938,0.0,0.938,wAfter[s]->GetMaximum());
    //    protonMass->SetLineStyle(8);
    //    protonMass->SetLineWidth(2);

    gPad->SetLogz(); 
    deltaBeta[s]->Draw("colz");
    upper->Draw(); 
    lower->Draw(); 

    if(s==0 || s==2 || s==4){
      gPad->SetLeftMargin(0.2);
      gPad->SetRightMargin(0.05);
      ytit.DrawLatex(0.1,0.48,"#Delta #beta");

      xtitx = 0.66;
    }

    else {
      gPad->SetRightMargin(0.2);
      gPad->SetLeftMargin(0.05);

      xtitx = 0.52;
    }

    if(s==0 || s==1){
      gPad->SetTopMargin(0.2);
      gPad->SetBottomMargin(0.05);
      xtity = 0.68;
    }
    else if (s==4 || s==5){
      gPad->SetTopMargin(0.05);
      gPad->SetBottomMargin(0.2);
      xtit.DrawLatex(0.46,0.1,"P (GeV/c)");
      xtity = 0.81;
    }
    else{
      gPad->SetTopMargin(0.05);
      gPad->SetBottomMargin(0.05);
      xtity = 0.81;
    }
  }

  can->cd();
  tit.DrawLatex(0.1, 0.95,"#Delta #beta histograms");


}
