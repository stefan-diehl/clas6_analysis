{

  TFile *f = TFile::Open("out.root");
  TTree *t = f->Get("events");

  kp_km = new TH1I("kp_km","",100,-10,70);
  kp_lambda = new TH1I("kp_lambda","",100,-10,70);

  t->Draw("angle_kp_km >> kp_km");
  t->Draw("angle_kp_lambda >> kp_lambda");

  kp_km->SetFillColorAlpha(99, 0.2);
  kp_lambda->SetFillColorAlpha(55, 0.2);

  TCanvas *c = new TCanvas("c","",1200,500);
  c->Divide(2, 1);

  c->cd(1); 
  kp_km->Draw();

  c->cd(2); 
  kp_lambda->Draw();

}
