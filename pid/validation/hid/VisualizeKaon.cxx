
/*

  David Riser 
  VisualizeKaon.cxx

  Draw histograms for 
  Kaon ID. 

  April 16, 2017 

 */

#include <iostream>
#include <vector>

#include "CommonTools.h"
#include "h22Option.h"
#include "MesonHistograms.h"
#include "Parameters.h"
#include "StandardHistograms.h"

#include "TCanvas.h"
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"

int PrintUsage(); 

class KaonVisualization {
public:
  KaonVisualization(){
  }

  ~KaonVisualization(){
  }

  void Configure(int argc, char *argv[]){
    opts.args["INPUT"].args = "UNSET";
    opts.args["INPUT"].type = 1;
    opts.args["INPUT"].name = "Input file";

    opts.args["PARS"].args = "/u/home/dmriser/Analysis_v2/lists/data.pars";
    opts.args["PARS"].type = 1;
    opts.args["PARS"].name = "Parameter file";

    opts.args["OUT"].args = "/volatile/clas12/dmriser/plots/pid/";
    opts.args["OUT"].type = 1;
    opts.args["OUT"].name = "Output folder";

    opts.set(argc, argv);
    
    kpHistos = new MesonHistograms("kp",321); 
    posHistos = new MesonHistograms("pos",321); 

    kmHistos = new MesonHistograms("km",-321); 
    negHistos = new MesonHistograms("neg",-321); 
    
    posCutDVzHistos = new MesonHistograms("posCutDVz", 321); 
    negCutDVzHistos = new MesonHistograms("negCutDVz", -321); 
    
    posCutDISHistos = new MesonHistograms("posCutDIS", 321); 
    negCutDISHistos = new MesonHistograms("negCutDIS", -321); 
    
    posCutDISMMHistos = new MesonHistograms("posCutDISMM", 321); 
    negCutDISMMHistos = new MesonHistograms("negCutDISMM", -321); 
    
    posCutDCHistos = new MesonHistograms("posCutDC", 321); 

    kpEvents = new StandardHistograms("kp",1); 
    kmEvents = new StandardHistograms("km",1); 

    std::string input = opts.args["INPUT"].args;
    if(input != "UNSET"){
      kpHistos->Load(input);
      kmHistos->Load(input);

      posHistos->Load(input);
      negHistos->Load(input);

      posCutDVzHistos->Load(input);
      posCutDCHistos->Load(input);
      posCutDISHistos->Load(input);
      posCutDISMMHistos->Load(input);

      negCutDVzHistos->Load(input);
      negCutDISHistos->Load(input);
      negCutDISMMHistos->Load(input);

      kpEvents->Load(input); 
      kmEvents->Load(input); 

      for(int s=0; s<7; s++){
	std::cout << "[KaonVisualization] Loaded Sector " << s << " with Entries = " << kpHistos->h1_p[s]->GetEntries() << std::endl;  
      }

    }
    else {
      std::cerr << "[KaonVisualization::Configure] Input file not set." << std::endl; 
      exit(0);
    }

    pars.loadParameters(opts.args["PARS"].args); 

    std::cout << "[KaonVisualization::Configure] Finished configuration successfully. " << std::endl;  
  }

  void Process(){
    std::string outputPath = opts.args["OUT"].args; 

    gStyle->SetOptStat(0); 
    gStyle->SetOptTitle(0); 
    //    Global::Visualization::SetCustomPalette();
    //    Global::Visualization::SetExpPalette();
    //    Global::Visualization::SetRandPalette(18081988);
    Global::Visualization::SetBentCoolWarmPalette();

    int colorIndex = 98; 

    TLatex title, xTitle, yTitle; 
    title.SetNDC(); 
    title.SetTextFont(102);
    title.SetTextSize(0.05);

    xTitle.SetNDC(); 
    xTitle.SetTextFont(102);
    xTitle.SetTextSize(0.04);

    yTitle.SetNDC(); 
    yTitle.SetTextFont(102);
    yTitle.SetTextSize(0.04);
    yTitle.SetTextAngle(90.0);

    TCanvas *can = new TCanvas("can","",1200,1000);
    //    TCanvas *can = new TCanvas("can","",800,500);
    TCanvas *singleCanvas = new TCanvas("singleCanvas","",800,500); 


    // -------------------------------------------------------
    //                  beta vs. p before any cuts 
    // -------------------------------------------------------
    for (int s=0; s<7; s++){
      singleCanvas->cd(); 

      posHistos->h2_p_beta[s]->Draw("colz"); 
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Positive Tracks, Cut Applied: None"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "#beta"); 

      singleCanvas->Print(Form("%skp/BetaPNoCutsSector%d.png",outputPath.c_str(), s));
      singleCanvas->Clear(); 
    }

    for (int s=0; s<7; s++){
      singleCanvas->cd(); 

      negHistos->h2_p_beta[s]->Draw("colz"); 
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Negative Tracks, Cut Applied: None"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "#beta"); 

      singleCanvas->Print(Form("%skm/BetaPNoCutsSector%d.png",outputPath.c_str(), s));
      singleCanvas->Clear(); 
    }

    // -------------------------------------------------------
    //                  beta vs. p after dvz cut  
    // -------------------------------------------------------
    for (int s=0; s<7; s++){
      singleCanvas->cd(); 

      posCutDVzHistos->h2_p_beta[s]->Draw("colz"); 
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Positive Tracks, Cut Applied: #DeltaV_{z}"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "#beta"); 

      singleCanvas->Print(Form("%skp/BetaPDVzCutSector%d.png",outputPath.c_str(), s));
      singleCanvas->Clear(); 
    }

    for (int s=0; s<7; s++){
      singleCanvas->cd(); 

      negCutDVzHistos->h2_p_beta[s]->Draw("colz"); 
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Negative Tracks, Cut Applied:  #DeltaV_{z}"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "#beta"); 

      singleCanvas->Print(Form("%skm/BetaPDVzCutSector%d.png",outputPath.c_str(), s));
      singleCanvas->Clear(); 
    }


    // -------------------------------------------------------
    //                  summary of cuts on beta vs. p
    // -------------------------------------------------------

    for(int s=0; s<7; s++){
      can->Divide(2,2); 

      can->cd(1);
      posHistos->h2_p_beta[s]->Draw("colz"); 
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: None"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "#beta"); 
      
      can->cd(2); 
      posCutDVzHistos->h2_p_beta[s]->Draw("colz"); 
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: #Deltav_{z}"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "#beta"); 

      can->cd(3); 
      posCutDCHistos->h2_p_beta[s]->Draw("colz"); 
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: DC R_{1} Fiducial"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "#beta"); 

      can->cd(4); 
      kpHistos->h2_p_beta[s]->Draw("colz");
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: All"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "#beta"); 

      can->Print(Form("%skp/BetaPSummarySector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }

    for(int s=0; s<7; s++){
      can->Divide(2,2); 

      can->cd(1);
      negHistos->h2_p_beta[s]->Draw("colz"); 
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: None"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "#beta"); 
      
      can->cd(2); 
      negCutDVzHistos->h2_p_beta[s]->Draw("colz"); 
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: #Deltav_{z}"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "#beta"); 

      can->cd(3); 
      negCutDISMMHistos->h2_p_beta[s]->Draw("colz"); 
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: DIS, M_{X}"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "#beta"); 

      can->cd(4); 
      kmHistos->h2_p_beta[s]->Draw("colz");
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: All"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "#beta"); 

      can->Print(Form("%skm/BetaPSummarySector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }

    // -------------------------------------------------------
    //                  summary of cuts on dc r1
    // -------------------------------------------------------
    for(int s=0; s<7; s++){
      can->Divide(2,2); 

      can->cd(1);
      posHistos->h2_dcx_dcy[s]->Draw("colz");
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: None"); 
      xTitle.DrawLatex(0.45, 0.02, "x (cm)"); 
      yTitle.DrawLatex(0.03, 0.48, "y (cm)"); 
      
      can->cd(2); 
      posCutDVzHistos->h2_dcx_dcy[s]->Draw("colz");
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: #Deltav_{z}"); 
      xTitle.DrawLatex(0.45, 0.02, "x (cm)"); 
      yTitle.DrawLatex(0.03, 0.48, "y (cm)"); 

      can->cd(3); 
      posCutDCHistos->h2_dcx_dcy[s]->Draw("colz");
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: DC R_{1} Fiducial"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "#beta"); 

      can->cd(4); 
      kpHistos->h2_dcx_dcy[s]->Draw("colz");
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: All"); 
      xTitle.DrawLatex(0.45, 0.02, "x (cm)"); 
      yTitle.DrawLatex(0.03, 0.48, "y (cm)"); 

      can->Print(Form("%skp/DCR1SummarySector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }

    // -------------------------------------------------------
    //                  summary of cuts on tofmass
    // -------------------------------------------------------
    for(int s=0; s<7; s++){
      can->Divide(2,2); 

      can->cd(1);
      posHistos->h2_p_tofmass[s]->Draw("colz");
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: None"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "M_{ToF}"); 
      
      can->cd(2); 
      posCutDVzHistos->h2_p_tofmass[s]->Draw("colz");
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: #Deltav_{z}"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "M_{ToF}"); 

      can->cd(3); 
      posCutDCHistos->h2_p_tofmass[s]->Draw("colz");
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: DC R_{1} Fiducial"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "M_{ToF}"); 

      can->cd(4); 
      kpHistos->h2_p_tofmass[s]->Draw("colz");
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: All"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "M_{ToF}"); 

      can->Print(Form("%skp/TofMassSummarySector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }

    for(int s=0; s<7; s++){
      can->Divide(2,2); 

      can->cd(1);
      negHistos->h2_p_tofmass[s]->Draw("colz");
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: None"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "M_{ToF}"); 
      
      can->cd(2); 
      negCutDVzHistos->h2_p_tofmass[s]->Draw("colz");
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: #Deltav_{z}"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "M_{ToF}"); 

      can->cd(3); 
      negCutDISMMHistos->h2_p_tofmass[s]->Draw("colz");
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: DIS, M_{X}"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "M_{ToF}"); 

      can->cd(4); 
      kmHistos->h2_p_tofmass[s]->Draw("colz");
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      title.DrawLatex(0.2, 0.93, "Cut Applied: All"); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "M_{ToF}"); 

      can->Print(Form("%skm/TofMassSummarySector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }
    

    // -------------------------------------------------------
    //                     p tofmass
    // -------------------------------------------------------
    for(int s=0; s<7; s++){
      negCutDVzHistos->h2_p_tofmass[s]->Draw("colz");
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      double cutUpper = pars.getParameter("KP_TOFMASS_MAX").getValue(0);
      double cutLower = pars.getParameter("KP_TOFMASS_MIN").getValue(0);

      TLine upper(0.5, cutUpper, 3.5, cutUpper); 
      TLine lower(0.5, cutLower, 3.5, cutLower); 

      upper.SetLineColor(99); 
      lower.SetLineColor(99); 

      upper.SetLineWidth(2); 
      lower.SetLineWidth(2); 

      upper.Draw(); 
      lower.Draw(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "M_{ToF}"); 

      can->Print(Form("%skm/TofMassSector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }

    // -------------------------------------------------------
    //                     p tofmass
    // -------------------------------------------------------
    for(int s=0; s<7; s++){
      posCutDVzHistos->h2_p_tofmass[s]->Draw("colz");
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      double cutUpper = pars.getParameter("KP_TOFMASS_MAX").getValue(0);
      double cutLower = pars.getParameter("KP_TOFMASS_MIN").getValue(0);

      TLine upper(0.5, cutUpper, 3.5, cutUpper); 
      TLine lower(0.5, cutLower, 3.5, cutLower); 

      upper.SetLineColor(99); 
      lower.SetLineColor(99); 

      upper.SetLineWidth(2); 
      lower.SetLineWidth(2); 

      upper.Draw(); 
      lower.Draw(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "M_{ToF}"); 

      can->Print(Form("%skp/TofMassSector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }

    // -------------------------------------------------------
    //                     p beta 
    // -------------------------------------------------------
    for(int s=0; s<6; s++){
      posCutDVzHistos->h2_p_beta[s+1]->Draw("colz");
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      TF1 *upper = new TF1("upper", "pol3", 0.5, 3.0); 
      TF1 *lower = new TF1("lower", "pol3", 0.5, 3.0); 
      TF1 *mid = new TF1("mid", "pol3", 0.5, 3.0); 

      mid->SetParameter(3, pars.getParameter("KP_BVP_MU_A").getValue(s)); 
      mid->SetParameter(2, pars.getParameter("KP_BVP_MU_B").getValue(s));
      mid->SetParameter(1, pars.getParameter("KP_BVP_MU_C").getValue(s));
      mid->SetParameter(0, pars.getParameter("KP_BVP_MU_D").getValue(s));
      
      upper->SetParameter(3, pars.getParameter("KP_BVP_MU_A").getValue(s) 
			  +pars.getParameter("KP_BVP_NSIGMA").getValue(0) 
			  *pars.getParameter("KP_BVP_SIGMA_A").getValue(s));
      upper->SetParameter(2, pars.getParameter("KP_BVP_MU_B").getValue(s) 
			  +pars.getParameter("KP_BVP_NSIGMA").getValue(0) 
			  *pars.getParameter("KP_BVP_SIGMA_B").getValue(s));
      upper->SetParameter(1, pars.getParameter("KP_BVP_MU_C").getValue(s) 
			  +pars.getParameter("KP_BVP_NSIGMA").getValue(0) 
			  *pars.getParameter("KP_BVP_SIGMA_C").getValue(s));
      upper->SetParameter(0, pars.getParameter("KP_BVP_MU_D").getValue(s) 
			  +pars.getParameter("KP_BVP_NSIGMA").getValue(0) 
			  *pars.getParameter("KP_BVP_SIGMA_D").getValue(s));

      lower->SetParameter(3, pars.getParameter("KP_BVP_MU_A").getValue(s) 
			  -pars.getParameter("KP_BVP_NSIGMA").getValue(0) 
			  *pars.getParameter("KP_BVP_SIGMA_A").getValue(s));
      lower->SetParameter(2, pars.getParameter("KP_BVP_MU_B").getValue(s) 
			  -pars.getParameter("KP_BVP_NSIGMA").getValue(0) 
			  *pars.getParameter("KP_BVP_SIGMA_B").getValue(s));
      lower->SetParameter(1, pars.getParameter("KP_BVP_MU_C").getValue(s) 
			  -pars.getParameter("KP_BVP_NSIGMA").getValue(0) 
			  *pars.getParameter("KP_BVP_SIGMA_C").getValue(s));
      lower->SetParameter(0, pars.getParameter("KP_BVP_MU_D").getValue(s) 
			  -pars.getParameter("KP_BVP_NSIGMA").getValue(0) 
			  *pars.getParameter("KP_BVP_SIGMA_D").getValue(s));

      upper->SetLineColor(99); 
      lower->SetLineColor(99); 
      mid->SetLineColor(kBlack); 

      upper->SetLineWidth(2); 
      lower->SetLineWidth(2); 
      mid->SetLineWidth(2); 
      mid->SetLineStyle(7); 

      upper->Draw("same"); 
      lower->Draw("same"); 
      mid->Draw("same");

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s+1)); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "#beta"); 

      can->Print(Form("%skp/BetaPSector%d.png",outputPath.c_str(), s+1));
      can->Clear(); 
    }

    for(int s=0; s<6; s++){
      negCutDVzHistos->h2_p_beta[s+1]->Draw("colz");
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      TF1 *upper = new TF1("upper", "pol3", 0.5, 3.0); 
      TF1 *lower = new TF1("lower", "pol3", 0.5, 3.0); 
      TF1 *mid = new TF1("mid", "pol3", 0.5, 3.0); 

      mid->SetParameter(3, pars.getParameter("KM_BVP_MU_A").getValue(s)); 
      mid->SetParameter(2, pars.getParameter("KM_BVP_MU_B").getValue(s));
      mid->SetParameter(1, pars.getParameter("KM_BVP_MU_C").getValue(s));
      mid->SetParameter(0, pars.getParameter("KM_BVP_MU_D").getValue(s));
      
      upper->SetParameter(3, pars.getParameter("KM_BVP_MU_A").getValue(s) 
			  +pars.getParameter("KM_BVP_NSIGMA").getValue(0) 
			  *pars.getParameter("KM_BVP_SIGMA_A").getValue(s));
      upper->SetParameter(2, pars.getParameter("KM_BVP_MU_B").getValue(s) 
			  +pars.getParameter("KM_BVP_NSIGMA").getValue(0) 
			  *pars.getParameter("KM_BVP_SIGMA_B").getValue(s));
      upper->SetParameter(1, pars.getParameter("KM_BVP_MU_C").getValue(s) 
			  +pars.getParameter("KM_BVP_NSIGMA").getValue(0) 
			  *pars.getParameter("KM_BVP_SIGMA_C").getValue(s));
      upper->SetParameter(0, pars.getParameter("KM_BVP_MU_D").getValue(s) 
			  +pars.getParameter("KM_BVP_NSIGMA").getValue(0) 
			  *pars.getParameter("KM_BVP_SIGMA_D").getValue(s));

      lower->SetParameter(3, pars.getParameter("KM_BVP_MU_A").getValue(s) 
			  -pars.getParameter("KM_BVP_NSIGMA").getValue(0) 
			  *pars.getParameter("KM_BVP_SIGMA_A").getValue(s));
      lower->SetParameter(2, pars.getParameter("KM_BVP_MU_B").getValue(s) 
			  -pars.getParameter("KM_BVP_NSIGMA").getValue(0) 
			  *pars.getParameter("KM_BVP_SIGMA_B").getValue(s));
      lower->SetParameter(1, pars.getParameter("KM_BVP_MU_C").getValue(s) 
			  -pars.getParameter("KM_BVP_NSIGMA").getValue(0) 
			  *pars.getParameter("KM_BVP_SIGMA_C").getValue(s));
      lower->SetParameter(0, pars.getParameter("KM_BVP_MU_D").getValue(s) 
			  -pars.getParameter("KM_BVP_NSIGMA").getValue(0) 
			  *pars.getParameter("KM_BVP_SIGMA_D").getValue(s));

      upper->SetLineColor(99); 
      lower->SetLineColor(99); 
      mid->SetLineColor(kBlack); 

      upper->SetLineWidth(2); 
      lower->SetLineWidth(2); 
      mid->SetLineWidth(2); 
      mid->SetLineStyle(7); 

      upper->Draw("same"); 
      lower->Draw("same"); 
      mid->Draw("same");

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s+1)); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "#beta"); 

      can->Print(Form("%skm/BetaPSector%d.png",outputPath.c_str(), s+1));
      can->Clear(); 
    }

    // -------------------------------------------------------
    //                      delta vz 
    // -------------------------------------------------------
    for(int s=0; s<7; s++){
      posHistos->h1_dvz[s]->SetFillColorAlpha(colorIndex, 0.3);
      posHistos->h1_dvz[s]->SetLineColor(colorIndex);
      posHistos->h1_dvz[s]->Draw();

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      
      double limit = pars.getParameter("KP_DVZ").getValue(0);
      TLine left(-1*limit,0,-1*limit,0.85*posHistos->h1_dvz[s]->GetMaximum()); 
      TLine right(limit,0,limit,0.85*posHistos->h1_dvz[s]->GetMaximum()); 

      left.SetLineWidth(2); 
      right.SetLineWidth(2); 

      left.SetLineColor(99); 
      right.SetLineColor(99); 

      left.Draw(); 
      right.Draw(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "#Deltav_{z} (cm)"); 
      yTitle.DrawLatex(0.03, 0.48, "Counts"); 

      can->Print(Form("%skp/DeltaVzSector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }
 
   for(int s=0; s<7; s++){
      negHistos->h1_dvz[s]->SetFillColorAlpha(colorIndex, 0.3);
      negHistos->h1_dvz[s]->SetLineColor(colorIndex);
      negHistos->h1_dvz[s]->Draw();

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 
      
      double limit = pars.getParameter("KP_DVZ").getValue(0);
      TLine left(-1*limit,0,-1*limit,0.85*posHistos->h1_dvz[s]->GetMaximum()); 
      TLine right(limit,0,limit,0.85*posHistos->h1_dvz[s]->GetMaximum()); 

      left.SetLineWidth(2); 
      right.SetLineWidth(2); 

      left.SetLineColor(99); 
      right.SetLineColor(99); 

      left.Draw(); 
      right.Draw(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "#Deltav_{z} (cm)"); 
      yTitle.DrawLatex(0.03, 0.48, "Counts"); 

      can->Print(Form("%skm/DeltaVzSector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }
  
    // -------------------------------------------------------
    //                      dc r1 
    // -------------------------------------------------------
    for(int s=0; s<7; s++){
      posHistos->h2_dcx_dcy[s]->Draw("colz");

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      double height = pars.getParameter("KP_DCR1_HEIGHT").getValue(0);
      double angle = pars.getParameter("KP_DCR1_ANGLE").getValue(0);
      double slope = 1/tan(0.5*to_radians*angle);

      TF1 upper("upper","(x-[0])/[1]",0,60); 
      upper.SetParameter(0, height); 
      upper.SetParameter(1, slope); 

      TF1 lower("lower","-1*(x-[0])/[1]",0,60); 
      lower.SetParameter(0, height); 
      lower.SetParameter(1, slope); 

      upper.SetLineWidth(2); 
      lower.SetLineWidth(2); 
      upper.SetLineColor(99); 
      lower.SetLineColor(99); 
      upper.Draw("same"); 
      lower.Draw("same"); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "x (cm)"); 
      yTitle.DrawLatex(0.03, 0.48, "y (cm)"); 

      can->Print(Form("%skp/DCR1Sector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }

    // -------------------------------------------------------
    //                     missing mass p 
    // -------------------------------------------------------
    for(int s=0; s<7; s++){
      kpHistos->h2_p_mm2[s]->Draw("colz");

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "M_{X}^{2}"); 

      can->Print(Form("%skp/MissingMassPSector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }

    // -------------------------------------------------------
    //                   prelim kinematics 
    // -------------------------------------------------------
    for(int s=0; s<7; s++){
      kpEvents->h1_phiH[s]->SetFillColorAlpha(colorIndex, 0.3);
      kpEvents->h1_phiH[s]->SetLineColor(colorIndex);
      kpEvents->h1_phiH[s]->Draw();

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "#phi_{h} (deg)"); 
      yTitle.DrawLatex(0.03, 0.48, "Counts"); 

      can->Print(Form("%skp/PhiHadronSector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }

    for(int s=0; s<7; s++){
      kmEvents->h1_phiH[s]->SetFillColorAlpha(colorIndex, 0.3);
      kmEvents->h1_phiH[s]->SetLineColor(colorIndex);
      kmEvents->h1_phiH[s]->Draw();

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "#phi_{h} (deg)"); 
      yTitle.DrawLatex(0.03, 0.48, "Counts"); 

      can->Print(Form("%skm/PhiHadronSector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }

    for(int s=0; s<7; s++){
      singleCanvas->cd(); 

      kpEvents->h1_mm2[s]->SetFillColorAlpha(colorIndex, 0.3);
      kpEvents->h1_mm2[s]->SetLineColor(colorIndex);
      kpEvents->h1_mm2[s]->Draw();

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "M_{X}^{2}"); 
      yTitle.DrawLatex(0.03, 0.48, "Counts"); 

      singleCanvas->Print(Form("%skp/MM2Sector%d.png",outputPath.c_str(), s));
      singleCanvas->Clear(); 
    }

    for(int s=0; s<7; s++){
      singleCanvas->cd(); 

      kmEvents->h1_mm2[s]->SetFillColorAlpha(colorIndex, 0.3);
      kmEvents->h1_mm2[s]->SetLineColor(colorIndex);
      kmEvents->h1_mm2[s]->Draw();

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "M_{X}^{2}"); 
      yTitle.DrawLatex(0.03, 0.48, "Counts"); 

      singleCanvas->Print(Form("%skm/MM2Sector%d.png",outputPath.c_str(), s));
      singleCanvas->Clear(); 
    }


    for(int s=0; s<7; s++){
      kpEvents->h1_pt2[s]->SetFillColorAlpha(colorIndex, 0.3);
      kpEvents->h1_pt2[s]->SetLineColor(colorIndex);
      kpEvents->h1_pt2[s]->Draw();

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "P_{T}^{2} (GeV^{2}/c^{2})"); 
      yTitle.DrawLatex(0.03, 0.48, "Counts"); 

      can->Print(Form("%skp/Pt2Sector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }


    for(int s=0; s<7; s++){
      kmEvents->h1_pt2[s]->SetFillColorAlpha(colorIndex, 0.3);
      kmEvents->h1_pt2[s]->SetLineColor(colorIndex);
      kmEvents->h1_pt2[s]->Draw();

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "P_{T}^{2} (GeV^{2}/c^{2})"); 
      yTitle.DrawLatex(0.03, 0.48, "Counts"); 

      can->Print(Form("%skm/Pt2Sector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }

    for(int s=0; s<7; s++){
      kpEvents->h1_z[s]->SetFillColorAlpha(colorIndex, 0.3);
      kpEvents->h1_z[s]->SetLineColor(colorIndex);
      kpEvents->h1_z[s]->Draw();

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "z"); 
      yTitle.DrawLatex(0.03, 0.48, "Counts"); 

      can->Print(Form("%skp/ZSector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }


    for(int s=0; s<7; s++){
      kmEvents->h1_z[s]->SetFillColorAlpha(colorIndex, 0.3);
      kmEvents->h1_z[s]->SetLineColor(colorIndex);
      kmEvents->h1_z[s]->Draw();

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "z"); 
      yTitle.DrawLatex(0.03, 0.48, "Counts"); 

      can->Print(Form("%skm/ZSector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }

    for(int s=0; s<7; s++){
      kpEvents->h1_q2[s]->SetFillColorAlpha(colorIndex, 0.3);
      kpEvents->h1_q2[s]->Draw();

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "Q^{2} (GeV^{2}/c^{2})"); 
      yTitle.DrawLatex(0.03, 0.48, "Counts"); 

      can->Print(Form("%skp/Q2Sector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }


    for(int s=0; s<7; s++){
      kmEvents->h1_q2[s]->SetFillColorAlpha(colorIndex, 0.3);
      kmEvents->h1_q2[s]->Draw();

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "Q^{2} (GeV^{2}/c^{2})"); 
      yTitle.DrawLatex(0.03, 0.48, "Counts"); 

      can->Print(Form("%skm/Q2Sector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }

    for(int s=0; s<7; s++){
      kpEvents->h2_z_pt2[s]->Draw("colz");

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "z"); 
      yTitle.DrawLatex(0.03, 0.48, "P_{T}^{2} (GeV^{2}/c^{2})"); 

      can->Print(Form("%skp/ZPt2Sector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }
 

    for(int s=0; s<7; s++){
      kmEvents->h2_z_pt2[s]->Draw("colz");

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "z"); 
      yTitle.DrawLatex(0.03, 0.48, "P_{T}^{2} (GeV^{2}/c^{2})"); 

      can->Print(Form("%skm/ZPt2Sector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }
 
    for(int s=0; s<7; s++){
      kpEvents->h2_part1_p_samplingFraction[s]->Draw("colz");

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "E_{dep}/p"); 

      can->Print(Form("%skp/KaonSFSector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }
 
 
    for(int s=0; s<7; s++){
      kmEvents->h2_part1_p_samplingFraction[s]->Draw("colz");

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetLogz(); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "p (GeV/c)"); 
      yTitle.DrawLatex(0.03, 0.48, "E_{dep}/p"); 

      can->Print(Form("%skm/KaonSFSector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }
 
    for(int s=0; s<7; s++){
      kmEvents->h1_ele_samplingFraction[s]->SetLineColor(99);
      kmEvents->h1_ele_samplingFraction[s]->SetFillColorAlpha(99, 0.2);
      kmEvents->h1_ele_samplingFraction[s]->Draw();

      kmEvents->h1_part1_samplingFraction[s]->SetLineColor(60);
      kmEvents->h1_part1_samplingFraction[s]->SetFillColorAlpha(60, 0.2);
      kmEvents->h1_part1_samplingFraction[s]->Draw("same");

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "E_{dep}/p"); 
      yTitle.DrawLatex(0.03, 0.48, "Counts"); 

      can->Print(Form("%skm/ElectronKaonSFSector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }
 
 
    for(int s=0; s<7; s++){
      kpEvents->h1_ele_samplingFraction[s]->SetLineColor(99);
      kpEvents->h1_ele_samplingFraction[s]->SetFillColorAlpha(99, 0.2);
      kpEvents->h1_ele_samplingFraction[s]->Draw();

      kpEvents->h1_part1_samplingFraction[s]->SetLineColor(60);
      kpEvents->h1_part1_samplingFraction[s]->SetFillColorAlpha(60, 0.2);
      kpEvents->h1_part1_samplingFraction[s]->Draw("same");

      gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
      gPad->SetGridx(); 
      gPad->SetGridy(); 

      title.DrawLatex(0.4, 0.93, Form("Sector %d",s)); 
      xTitle.DrawLatex(0.45, 0.02, "E_{dep}/p"); 
      yTitle.DrawLatex(0.03, 0.48, "Counts"); 

      can->Print(Form("%skp/ElectronKaonSFSector%d.png",outputPath.c_str(), s));
      can->Clear(); 
    }
 
  }

private:
  h22Options          opts; 
  Parameters          pars; 
  MesonHistograms    *kpHistos; 
  MesonHistograms    *kmHistos; 

  MesonHistograms    *posHistos; 
  MesonHistograms    *negHistos; 

  MesonHistograms    *posCutDVzHistos; 
  MesonHistograms    *negCutDVzHistos; 

  MesonHistograms    *posCutDISHistos; 
  MesonHistograms    *negCutDISHistos; 

  MesonHistograms    *posCutDISMMHistos; 
  MesonHistograms    *negCutDISMMHistos; 

  MesonHistograms    *posCutDCHistos; 

  StandardHistograms *kpEvents; 
  StandardHistograms *kmEvents; 
};

int main(int argc, char *argv[]){

  if (argc < 2){
    return PrintUsage();
  }
  
  KaonVisualization kViz; 
  kViz.Configure(argc, argv); 
  kViz.Process();

  return 0;
}

int PrintUsage(){

  std::cout << "Please pass in the file containing MesonHistograms created by kaon.cxx" << std::endl; 
  std::cout << "with the flag ./VisualizeKaon -INPUT=kaon.root -PARS=kaon.pars -OUT=images/" << std::endl; 

  return -1;
}
