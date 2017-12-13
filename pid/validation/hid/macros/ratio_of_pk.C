{
  
  // ------------------------------------------------
  //     User Stuff 
  // ------------------------------------------------
  TFile *f = TFile::Open("/volatile/clas12/dmriser/rootFiles/pid/meson/validation.root"); 
  std::string outputPath("/volatile/clas12/dmriser/plots/pid/kp/"); 

  // ------------------------------------------------
  //     Other Stuff 
  // ------------------------------------------------
  TH1D *xPip = (TH1D*) f->Get("StandardHistograms/h1_xbj_sect0_pip"); 
  TH1D *xKp  = (TH1D*) f->Get("StandardHistograms/h1_xbj_sect0_kp"); 

  xPip->Scale(1/xPip->Integral());
  xKp->Scale(1/xKp->Integral());
  
  TH1D *ratio = xKp->Clone(); 
  ratio->Divide(xPip); 
  ratio->SetMaximum(1.0);
  ratio->SetMarkerStyle(8); 
  ratio->Draw("p"); 
}
