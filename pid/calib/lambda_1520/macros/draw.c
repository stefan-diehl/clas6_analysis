{

  TFile *input = TFile::Open("../out/out.root");

  TH1I *mm_ek[4];
  mm_ek[0] = (TH1I*) input->Get("mm_ek_no_cut");
  mm_ek[1] = (TH1I*) input->Get("mm_ek_prot");
  mm_ek[2] = (TH1I*) input->Get("mm_ek_pip");
  mm_ek[3] = (TH1I*) input->Get("mm_ek_pim");


  TH1I *mm_ekpion[2];
  mm_ekpion[0] = (TH1I*) input->Get("mm_ekpion_plus");
  mm_ekpion[1] = (TH1I*) input->Get("mm_ekpion_minus");

 TH1I *mm_ekproton = (TH1I*) input->Get("mm_ekproton");

  mm_ek[0]->SetLineColor(kBlack);
  mm_ek[0]->SetMinimum(0.0);

  for(int i=1; i<4; i++){
    mm_ek[i]->SetLineColor(kBlack);
    mm_ek[i]->SetFillColor(kBlack);
    mm_ek[i]->SetFillStyle(3004);
  }

  TLatex *tit = new TLatex();
  tit->SetNDC();
  tit->SetTextFont(102);
  tit->SetTextSize(0.05);

  TCanvas *can = new TCanvas("can","",1600,1200); 
  can->Divide(2,2);

  can->cd(1);
  mm_ek[0]->Draw();
  tit->DrawLatex(0.45, 0.05, "M_{eK^{+}}");

  can->cd(2);
  //  mm_ek[0]->Draw();
  mm_ek[1]->Draw();
  tit->DrawLatex(0.35, 0.05, "M_{eK^{+}} (with p^{+} K^{-})");

  can->cd(3);
  //  mm_ek[0]->Draw();
  mm_ek[2]->Draw();
  tit->DrawLatex(0.35, 0.05, "M_{eK^{+}} (with #pi^{+} #Sigma^{-})");

  can->cd(4);
  //  mm_ek[0]->Draw();
  mm_ek[3]->Draw();
  tit->DrawLatex(0.35, 0.05, "M_{eK^{+}} (with #pi^{-} #Sigma^{+})");

  can->Print("mm_ek.pdf");

  TCanvas *can2 = new TCanvas("can2","",1600,1200);
  can2->cd();
  
  can2->Divide(2,2);
  
  can2->cd(1);
  mm_ekpion[0]->Draw(); 
  tit->DrawLatex(0.35, 0.05, "M_{eK^{+}#pi^{+}}");

  can2->cd(2);
  mm_ekpion[1]->Draw(); 
  tit->DrawLatex(0.35, 0.05, "M_{eK^{+}#pi^{-}}");

  can2->cd(3);
  mm_ekproton->Draw(); 
  tit->DrawLatex(0.35, 0.05, "M_{eK^{+}p^{+}}");

  can2->Print("mm_other_channels.pdf");
}
