{

  TFile *input = TFile::Open("/volatile/clas12/dmriser/farm_out/phi_1020_pass0/phi.root");

  TH1F *mm_epx[4];
  mm_epx[0] = (TH1F*) input->Get("mm_epx_no_cuts");
  mm_epx[1] = (TH1F*) input->Get("mm_epx_tag_kp");
  mm_epx[2] = (TH1F*) input->Get("mm_epx_tag_kp_cut_km_mm");
  mm_epx[3] = (TH1F*) input->Get("mm_epx_tag_kp_cut_im_kk");

  TH1F *mm_epkx[2];
  mm_epkx[0] = (TH1F*) input->Get("mm_epkx_no_cut");
  mm_epkx[1] = (TH1F*) input->Get("mm_epkx_cut_mm_phi");

  TH1I *im_kk = input->Get("im_kk"); 

  TLatex *label = new TLatex();
  label->SetNDC(); 
  label->SetTextFont(42); 
  label->SetTextSize(0.04);

  TCanvas *can = new TCanvas("can", "", 1600, 1200);
  can->Divide(2, 2);

  can->cd(1);
  mm_epx[0]->SetFillStyle(3004); 
  mm_epx[0]->SetFillColorAlpha(kBlack, 1.0); 
  mm_epx[0]->SetLineColor(kBlack);
  mm_epx[0]->Draw();
  label->DrawLatex(0.4, 0.05, "M_{X} (ep #rightarrow epX)");

  can->cd(2);
  mm_epx[1]->SetFillStyle(3004); 
  mm_epx[1]->SetFillColorAlpha(kBlack, 1.0); 
  mm_epx[1]->SetLineColor(kBlack);
  mm_epx[1]->Draw();
  label->DrawLatex(0.3, 0.05, "M_{X} (ep #rightarrow epX) +Tagged K^{+}");

  can->cd(3);
  mm_epx[2]->SetFillStyle(3004); 
  mm_epx[2]->SetFillColorAlpha(kBlack, 1.0); 
  mm_epx[2]->SetLineColor(kBlack);
  mm_epx[2]->Draw();
  label->DrawLatex(0.15, 0.05, "M_{X} (ep #rightarrow epX) +Tagged K^{+} and M_{X} (ep #rightarrow epkX) = K^{-}");

  can->cd(4);
  mm_epx[3]->SetFillStyle(3004); 
  mm_epx[3]->SetFillColorAlpha(kBlack, 1.0); 
  mm_epx[3]->SetLineColor(kBlack);
  mm_epx[3]->Draw();
  label->DrawLatex(0.15, 0.05, "M_{X} (ep #rightarrow epX) (3) + cut on IM_{KK} = M_{\phi}");

  // print this mother 
  can->Print("mm_epx.pdf");

  TCanvas *can2 = new TCanvas("can2", "", 1600, 600);
  can2->Divide(2, 1);

  can2->cd(1); 
  mm_epkx[0]->SetFillStyle(3004);
  mm_epkx[0]->SetFillColorAlpha(kBlack, 1.0);
  mm_epkx[0]->SetLineColor(kBlack);
  mm_epkx[0]->Draw();   
  label->DrawLatex(0.35, 0.05, "M_{X} (ep #rightarrow epK^{+}X)"); 

  can2->cd(2); 
  mm_epkx[1]->SetFillStyle(3004);
  mm_epkx[1]->SetFillColorAlpha(kBlack, 1.0);
  mm_epkx[1]->SetLineColor(kBlack);
  mm_epkx[1]->Draw();   
  label->DrawLatex(0.18, 0.05, "M_{X} (ep #rightarrow epK^{+}X) + cut on M_{X} (ep #rightarrow epX) = M_{\phi}"); 

  can2->Print("mm_epkx.pdf");

  
  TCanvas *can3 = new TCanvas("can3", "", 800, 600);
  im_kk->SetFillStyle(3004); 
  im_kk->SetFillColorAlpha(kBlack, 1.0); 
  im_kk->SetLineColor(kBlack); 
  im_kk->Draw(); 
  label->DrawLatex(0.15, 0.05, "IM_{K^{+} K^{-}} where K^{-} is inferred from M_{X} (ep #rightarrow epkx)"); 

  can3->Print("im_kk.pdf");
}

