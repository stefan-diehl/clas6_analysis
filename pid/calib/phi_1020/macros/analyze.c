{

  TFile *inputFile = TFile::Open("/volatile/clas12/dmriser/farm_out/phi_1020_pass8/pid.root");
  //  TFile *inputFile = TFile::Open("../out/big_dingus.root");
  TTree *events    = (TTree*) inputFile->Get("events");

  TH1F *mm_epkx_nocut = new TH1F("mm_epkx_nocut", "", 200, 0.3, 0.7);
  TH1F *mm_epkx_boxcut = new TH1F("mm_epkx_boxcut", "", 200, 0.3, 0.7);
  TH1F *mm_epkx_allcut = new TH1F("mm_epkx_allcut", "", 200, 0.3, 0.7);

  TH1F *mm_ekx_nocut = new TH1F("mm_ekx_nocut", "", 200, 0.0, 2.0);
  TH1F *mm_ekx_boxcut = new TH1F("mm_ekx_boxcut", "", 200, 0.0, 2.0);

  TH1F *mm_epx_nocut  = new TH1F("mm_epx_nocut", "", 200, 0.5, 2.1);
  TH1F *mm_epx_allcut = new TH1F("mm_epx_allcut", "", 200, 0.5, 2.1);
  TH1F *mm_epx_allcut = new TH1F("mm_epx_anglecut", "", 200, 0.5, 2.1);

  TH1F *im_kk_nocut = new TH1F("im_kk_nocut", "",     60, 0.95, 1.15); 
  TH1F *im_kk_epkxcut = new TH1F("im_kk_epkxcut", "", 60, 0.95, 1.15); 
  TH1F *im_kk_allcut = new TH1F("im_kk_allcut", "",   60, 0.95, 1.15); 
  TH1F *im_kk_lambdacut = new TH1F("im_kk_lambdacut", "",   60, 0.95, 1.15); 

  TH2F *mm_compare_mass = new TH2F("mm_compare_mass", "mm_compare_mass", 200, 0.0, 1.4, 200, 0.0, 1.4);
  TH2F *mm_epx_epkx = new TH2F("mm_epx_epkx", "", 200, 0.8, 2.0, 200, 0.0, 1.2); 
  TH2F *mm_epkx_im_kk = new TH2F("mm_epkx_im_kk", "", 200, 0.8, 1.3, 200, 0.3, 0.7); 

  TH2F *im_kk_theta_kk = new TH2F("im_kk_theta_kk", "",200, 0.90, 1.6, 200, 0.0, 75.0);   
  TH2F *mm_epx_theta_kk = new TH2F("mm_epx_theta_kk", "",200, 0.50, 1.6, 200, 0.0, 75.0); 
  TH2F *im_kk_im_pk = new TH2F("im_kk_im_pk", "", 200, 0.9, 1.3, 200, 0.9, 2.0); 
  
  // real right 
  std::string pidcut("alpha_kp > 0.32 && alpha_prot > 0.32");
  std::string boxcut("missing_mass_epkx_mpi > 0.4 && missing_mass_epkx_mk > 0.2");
  std::string epkxcut("missing_mass_epkx_mk > 0.475 && missing_mass_epkx_mk < 0.518");
  std::string epxcut("missing_mass_epx > 0.980 && missing_mass_epx < 1.060");
  std::string anglecut("theta_kk_lab > 5.0 && theta_kk_lab < 12.0");
  std::string lambdacut("im_pk < 1.48 || im_pk > 1.56");  

  // do projections 
  events->Draw("missing_mass_epkx_mk >> mm_epkx_nocut", pidcut.c_str(), "");
  events->Draw("missing_mass_epkx_mk >> mm_epkx_boxcut", Form("%s && %s", boxcut.c_str(), pidcut.c_str()), "");
  events->Draw("missing_mass_epkx_mk >> mm_epkx_allcut", Form("%s && %s && %s && %s", pidcut.c_str(), boxcut.c_str(), epkxcut.c_str(), epxcut.c_str()), "");

  events->Draw("missing_mass_ekx_mk >> mm_ekx_nocut", pidcut.c_str(), "");
  events->Draw("missing_mass_ekx_mk >> mm_ekx_boxcut", Form("%s && %s", pidcut.c_str(), boxcut.c_str()), "");

  events->Draw("missing_mass_epx >> mm_epx_nocut", pidcut.c_str(), "");
  events->Draw("missing_mass_epx >> mm_epx_allcut", Form("%s && %s && %s", pidcut.c_str(), boxcut.c_str(), epkxcut.c_str()) , "");
  events->Draw("missing_mass_epx >> mm_epx_anglecut", Form("%s && %s && %s && %s", pidcut.c_str(), anglecut.c_str(), boxcut.c_str(), epkxcut.c_str()) , "");

  events->Draw("im_kk >> im_kk_nocut", pidcut.c_str(), ""); 
  events->Draw("im_kk >> im_kk_allcut", Form("%s && %s && %s", pidcut.c_str(), boxcut.c_str(), epkxcut.c_str()), "");
  events->Draw("im_kk >> im_kk_lambdacut", Form("%s && %s && %s && %s", pidcut.c_str(), lambdacut.c_str(), boxcut.c_str(), epkxcut.c_str()), "");
  //  events->Draw("im_kk >> im_kk_allcut", Form("%s && %s && %s", boxcut.c_str(), epkxcut.c_str(), epxcut.c_str()), "");
  events->Draw("im_kk >> im_kk_epkxcut", Form("%s && %s", pidcut.c_str(), epkxcut.c_str()), "");

  events->Draw("missing_mass_epkx_mpi:missing_mass_epkx_mk >> mm_compare_mass", pidcut.c_str(), "colz");
  events->Draw("missing_mass_epkx_mk:missing_mass_epx >> mm_epx_epkx", Form("%s && %s", pidcut.c_str(), boxcut.c_str()), "colz");

  //  events->Draw("missing_mass_epkx_mk:im_kk >> mm_epkx_im_kk", Form("%s && %s", boxcut.c_str(), epxcut.c_str()), "colz"); 
  events->Draw("missing_mass_epkx_mk:im_kk >> mm_epkx_im_kk", boxcut.c_str(), "colz"); 

  events->Draw("theta_kk_lab:im_kk >> im_kk_theta_kk", Form("%s && %s && %s", lambdacut.c_str(), boxcut.c_str(), epkxcut.c_str()));
  //  events->Draw("theta_kk_lab:missing_mass_epx >> mm_epx_theta_kk", Form("%s && %s", boxcut.c_str(), epkxcut.c_str()));
  events->Draw("theta_kk_lab:missing_mass_epx >> mm_epx_theta_kk", Form("%s && %s", pidcut.c_str(),boxcut.c_str()));
  events->Draw("im_pk:im_kk >> im_kk_im_pk", Form("%s && %s && %s", pidcut.c_str(), boxcut.c_str(), epkxcut.c_str()), "colz");


  //drawing shits 
  TLatex *lab = new TLatex(); 
  lab->SetNDC();
  lab->SetTextSize(0.04); 
  lab->SetTextFont(42); 

  TCanvas *can = new TCanvas("can", "", 1200, 500); 
  can->Divide(2,1); 

  can->cd(1); 
  mm_epkx_nocut->SetLineColor(kBlack);
  mm_epkx_nocut->SetFillColor(kBlack);
  mm_epkx_nocut->SetFillStyle(3004);
  mm_epkx_nocut->Draw(); 
  lab->DrawLatex(0.4, 0.05, "M_{X} (epkx)"); 

  can->cd(2); 
  mm_epkx_im_kk->Draw("colz"); 
  lab->DrawLatex(0.4, 0.05, "IM_{KK}"); 

  lab->SetTextAngle(90.0); 
  lab->DrawLatex(0.05, 0.45, "M_{X} (epkx)"); 
  lab->SetTextAngle(0.0); 

  can->Print("mm_epkx.pdf"); 
  can->Clear();

  // 0-------------------------------0 
  //             next plot 
  // 0-------------------------------0 
  can->Divide(2,1);
  can->cd(1);
  mm_compare_mass->Draw("colz"); 
  
  TLine *hline = new TLine(0.2, 0.0, 0.2, 1.4);
  TLine *vline = new TLine(0.0, 0.4, 1.4, 0.4);
  
  hline->SetLineColor(kWhite); 
  vline->SetLineColor(kWhite); 

  hline->Draw(); 
  vline->Draw(); 
  lab->DrawLatex(0.3, 0.05, "M_{X} (epkX) assume kaon mass"); 

  lab->SetTextAngle(90.0); 
  lab->DrawLatex(0.025, 0.2, "M_{X} (epkX) assume pion mass"); 
  lab->SetTextAngle(0.0); 

  can->cd(2); 
  mm_ekx_nocut->SetLineColor(kBlack);   
  mm_ekx_nocut->Draw(); 

  mm_ekx_boxcut->SetLineColor(kBlack); 
  mm_ekx_boxcut->SetFillColor(kBlack); 
  mm_ekx_boxcut->SetFillStyle(3004); 
  mm_ekx_boxcut->Draw("same"); 
  lab->DrawLatex(0.2, 0.05, "M_{X} (ekx) with/without cut on left"); 

  can->Print("mm_ekx.pdf"); 

  // 0-------------------------------0 
  //             next plot 
  // 0-------------------------------0 

  can->Clear(); 
  can->Divide(2,1); 
  
  can->cd(1);
  mm_epx_epkx->Draw("colz");
  lab->DrawLatex(0.3, 0.05, "M_{X} (epx)"); 

  lab->SetTextAngle(90.0); 
  lab->DrawLatex(0.05, 0.3, "M_{X} (epkx)"); 
  lab->SetTextAngle(0.0); 

  TLine *hcut_top   = new TLine(0.8, 0.47, 2.0, 0.47); 
  TLine *hcut_bot   = new TLine(0.8, 0.53, 2.0, 0.53); 
  TLine *vcut_left  = new TLine(0.98, 0.0, 0.98, 1.2); 
  TLine *vcut_right = new TLine(1.06, 0.0, 1.06, 1.2); 

  hcut_top  ->SetLineColor(kWhite); 
  hcut_bot  ->SetLineColor(kWhite); 
  vcut_left ->SetLineColor(kWhite); 
  vcut_right->SetLineColor(kWhite); 

  hcut_bot  ->Draw(); 
  hcut_top  ->Draw(); 
  vcut_left ->Draw(); 
  vcut_right->Draw(); 


  can->cd(2); 

  im_kk_lambdacut->SetLineColor(kBlack); 
  im_kk_lambdacut->SetFillColor(kBlack); 
  im_kk_lambdacut->SetFillStyle(3004); 
  im_kk_lambdacut->SetMinimum(0.0); 
  im_kk_lambdacut->Draw();   

  im_kk_allcut->SetLineColor(99); 
  //  im_kk_allcut->SetFillColor(kBlack); 
  //  im_kk_allcut->SetFillStyle(3004); 
  im_kk_allcut->SetMinimum(0.0); 
  im_kk_allcut->Draw("same");   

  lab->DrawLatex(0.3, 0.05, "IM_{KK} all other cuts applied"); 

  can->Print("im_kk.pdf"); 

  // 0-------------------------------0 
  //  \          next plot          \
  // 0-------------------------------0 
  TCanvas *can2 = new TCanvas("can2", "", 1200, 700); 

  mm_epx_nocut->SetMinimum(0.0);
  mm_epx_nocut->SetLineColor(kBlack); 
  mm_epx_nocut->Draw();
  
  mm_epx_allcut->SetLineColor(kBlack);
  mm_epx_allcut->SetFillColor(kBlack);
  mm_epx_allcut->SetFillStyle(3004);
  mm_epx_allcut->Draw("same");

  //  mm_epx_anglecut->SetLineColor(99); 
  //  mm_epx_anglecut->SetFillColor(99); 
  //  mm_epx_anglecut->SetFillStyle(3004); 
  //  mm_epx_anglecut->Draw("same"); 

  lab->DrawLatex(0.3, 0.05, "M_{X} (epx) w + w/o cut on epkx (k-)");

  can2->Print("mm_epx.pdf");

  // 0-------------------------------0 
  //  \          next plot          \
  // 0-------------------------------0 
  can2->Clear(); 
  im_kk_theta_kk->Draw("colz"); 

  lab->DrawLatex(0.45, 0.05, "IM_{KK}");
  lab->SetTextAngle(90.0); 
  lab->DrawLatex(0.05, 0.45, "#theta_{KK}");
  lab->SetTextAngle(0.0); 

  can2->Print("im_kk_theta_kk.pdf");

  // 0-------------------------------0 
  //  \          next plot          \
  // 0-------------------------------0 
  can2->Clear(); 
  mm_epx_theta_kk->Draw("colz"); 

  lab->DrawLatex(0.45, 0.05, "M_{X} (epx)");
  lab->SetTextAngle(90.0); 
  lab->DrawLatex(0.05, 0.45, "#theta_{KK}");
  lab->SetTextAngle(0.0); 

  can2->Print("mm_epx_theta_kk.pdf");

  // 0-------------------------------0 
  //  \          next plot          \
  // 0-------------------------------0 
  can2->Clear(); 
  im_kk_im_pk->Draw("colz"); 

  lab->DrawLatex(0.45, 0.05, "IM_{KK}");
  lab->SetTextAngle(90.0); 
  lab->DrawLatex(0.05, 0.45, "IM_{PK^{-}}");
  lab->SetTextAngle(0.0); 

  can2->Print("im_kk_im_pk.pdf");
}
