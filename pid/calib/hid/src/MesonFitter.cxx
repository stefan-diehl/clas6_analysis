////////////////////////////////////////
/*

 David Riser, University of Connecticut
 
 July 13, 2016
 
 Modified:
 March 8, 2017
 April 14, 2017 

 Make slices 
 
 */
////////////////////////////////////////

#include <iostream>

#include "CommonTools.h"
#include "h22Option.h"
#include "MesonHistograms.h"
#include "Parameters.h"
#include "ParameterSet.h"
#include "SliceFitter.h"
#include "TF1Integrator.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"


class MesonBetaSliceFitter {
public: 
  MesonBetaSliceFitter(MesonHistograms *h, int n, double min, double max) : histos(h) {
    fitter = new SliceFitter(n, min, max); 
  } 

  void Fit(std::string name){

    std::string meanString(Form("x/sqrt(x^2 + %f)",pow(pid_to_mass(histos->GetPid()),2))); 
    fitter->SetExpectedMean(meanString); 
    fitter->SetLowerTolerance(0.96);
    fitter->SetUpperTolerance(1.02); 
    
    for(int s=1; s<7; s++){
      std::string title(Form("%s_sect%d",name.c_str(),s)); 
      fitter->Fit(histos->h2_p_beta[s], title); 

      std::cout << "Fit slice " << s << std::endl; 
      
      muGraph[s-1] = fitter->GetGraphMu(title);
      sigmaGraph[s-1] = fitter->GetGraphSigma(title); 

      muFit[s-1] = fitter->GetFitToMu("pol3",title); 
      sigmaFit[s-1] = fitter->GetFitToSigma("pol3",title); 

      std::vector<TF1> theseFits = fitter->GetFits(); 
      std::vector<TH1D> theseSlices = fitter->GetSlices(); 

      std::cout << "[MesonBetaSliceFitter::Fit] For sector " << s << " we have slices = " << theseSlices.size() << std::endl;
      std::cout << "[MesonBetaSliceFitter::Fit] For sector " << s << " we have fits = "   << theseFits.size() << std::endl;

      fits.push_back(theseFits);
      slices.push_back(theseSlices); 
    }
    
  }

  std::string GetParticleName(){
    std::string particle; 

    int pid = histos->GetPid(); 

    if (pid == 211){
      particle = "PIP"; 
    } else if (pid == -211){
      particle = "PIM"; 
    } else if (pid == 321){
      particle = "KP"; 
    } else if (pid == -321){
      particle = "KM"; 
    } else {
      particle = "unknown"; 
    }

    return particle; 
  }

  void WriteParameters(Parameters *pars){
    std::string particle = GetParticleName(); 

    ParameterSet muA, muB, muC, muD;
    ParameterSet sigmaA, sigmaB, sigmaC, sigmaD;
    ParameterSet nSigma; 

    muA.setName(Form("%s_BVP_MU_A",particle.c_str())); 
    muB.setName(Form("%s_BVP_MU_B",particle.c_str())); 
    muC.setName(Form("%s_BVP_MU_C",particle.c_str())); 
    muD.setName(Form("%s_BVP_MU_D",particle.c_str())); 
    sigmaA.setName(Form("%s_BVP_SIGMA_A",particle.c_str())); 
    sigmaB.setName(Form("%s_BVP_SIGMA_B",particle.c_str())); 
    sigmaC.setName(Form("%s_BVP_SIGMA_C",particle.c_str())); 
    sigmaD.setName(Form("%s_BVP_SIGMA_D",particle.c_str())); 
    nSigma.setName(Form("%s_BVP_NSIGMA",particle.c_str())); 

    nSigma.addValueAndError(3.0, 0.0); 

    for (int s=0; s<6; s++){
      muA.addValueAndError(muFit[s].GetParameter(3), muFit[s].GetParError(3)); 
      muB.addValueAndError(muFit[s].GetParameter(2), muFit[s].GetParError(2)); 
      muC.addValueAndError(muFit[s].GetParameter(1), muFit[s].GetParError(1)); 
      muD.addValueAndError(muFit[s].GetParameter(0), muFit[s].GetParError(0)); 

      sigmaA.addValueAndError(sigmaFit[s].GetParameter(3), sigmaFit[s].GetParError(3)); 
      sigmaB.addValueAndError(sigmaFit[s].GetParameter(2), sigmaFit[s].GetParError(2)); 
      sigmaC.addValueAndError(sigmaFit[s].GetParameter(1), sigmaFit[s].GetParError(1)); 
      sigmaD.addValueAndError(sigmaFit[s].GetParameter(0), sigmaFit[s].GetParError(0)); 
    }

    pars->addParameterSet(muA); 
    pars->addParameterSet(muB); 
    pars->addParameterSet(muC); 
    pars->addParameterSet(muD); 
    pars->addParameterSet(sigmaA); 
    pars->addParameterSet(sigmaB); 
    pars->addParameterSet(sigmaC); 
    pars->addParameterSet(sigmaD); 
    pars->addParameterSet(nSigma); 
  }

  void Save(TFile *out){
    if(out->IsOpen()){
      for(int s=0; s<6; s++){
	muGraph[s].Write();
	sigmaGraph[s].Write(); 

	muFit[s].Write();
	sigmaFit[s].Write(); 

	for(int b=0; b<fits[s].size(); b++){
	  slices[s][b].Write(); 
	  fits[s][b].Write(); 
	}
      }

    } else {
      std::cerr << "[MesonBetaSliceFitter::Save] Output file is not open. " << std::endl; 
    }
    
  } 
  
protected:
  MesonHistograms *histos; 
  SliceFitter     *fitter; 

  TGraphErrors muGraph[6];
  TGraphErrors sigmaGraph[6]; 

  TF1 muFit[6]; 
  TF1 sigmaFit[6]; 
  
  std::vector<std::vector<TF1> >  fits;
  std::vector<std::vector<TH1D> > slices; 
  
};

class MesonSlices { 
public:

  MesonSlices(){

  }

  ~MesonSlices(){
  }

  std::vector<std::vector<TH1D*> > beta; 
  std::vector<std::vector<TH1D*> > deltaBeta; 
  std::vector<std::vector<TH1D*> > mass; 
  std::vector<TH1D*>               mass_projection;     // Doesn't show everlap in signal.
  std::vector<TF1*>                mass_projection_fit; // Doesn't really work very well. 
  std::vector<std::vector<TF1*> >  mass_fit;
  std::vector<std::vector<TF1*> >  beta_fit;
  std::vector<double> bins; 

  void SetBins(std::vector<double> b){
    bins = b; 
  }

  void Save(TFile *outputFile){
    if (outputFile->IsOpen()){
      outputFile->mkdir("slices_beta");
      outputFile->mkdir("slices_dbeta");
      outputFile->mkdir("slices_mass");

      outputFile->cd();
      outputFile->cd("slices_beta");
      for(int s=0; s<beta.size(); s++){
	for(int b=0; b<beta[s].size(); b++){
	  beta[s][b]->Write();
	  beta_fit[s][b]->Write();
	}
      }

      outputFile->cd();
      outputFile->cd("slices_dbeta");
      for(int s=0; s<deltaBeta.size(); s++){
	for(int b=0; b<deltaBeta[s].size(); b++){
	  deltaBeta[s][b]->Write();
	}
      }

      outputFile->cd();
      outputFile->cd("slices_mass");
      for(int s=0; s<mass.size(); s++){
	for(int b=0; b<mass[s].size(); b++){
	  mass[s][b]->Write();
	  mass_fit[s][b]->Write();
	}
      }

      for(int s=0; s<mass_projection.size(); s++){
	mass_projection[s]->Write(); 
	mass_projection_fit[s]->Write(); 
      }

      outputFile->cd();
    } else {
      std::cerr << "[MesonFittingService::Save] The output file was not opened successfully. " << std::endl; 
    }
  }

};

class MesonFittingService {
  
public:

  MesonFittingService(MesonHistograms *h) : histos(h) {
    slices = new MesonSlices(); 

    //    fPMin = 0.5; 
    //    fPMax = 3.5; 
    fPMin = 0.5; 
    fPMax = 1.8; 
    fPWidth = (fPMax-fPMin)/fNSlices; 

    for(int b=0; b<fNSlices; b++){
      fBins.push_back(fPMin +fPWidth*b);
    }

    slices->SetBins(fBins); 
  }

  ~MesonFittingService(){
  }

  MesonHistograms   *histos;   
  MesonSlices       *slices; 

  void SetPMin(double min) { fPMin = min; }
  void SetPMax(double max) { fPMax = max; }

  const int GetNumberOfSlices() const {
    return fNSlices; 
  }

  double GetBin(int index) {
    if(index < fNSlices){
      return fBins[index];
    } else {
      return 0.00; 
    }
  }

  double GetPMin() const { 
    return fPMin; 
  }
 
  double GetPMax() const { 
    return fPMax; 
  }
  
  void Execute(){
    Slice();
  }

  void Save(TFile *out){
    slices->Save(out);
  }

protected:
  const static int fNSlices = 24; 
  std::vector<double> fBins; 
  double fPMin, fPMax, fPWidth; 

  void Slice(){
    for(int s=0; s<6; s++){
      std::vector<TH1D*> temp_beta; 
      std::vector<TH1D*> temp_deltaBeta; 
      std::vector<TH1D*> temp_mass; 
      std::vector<TF1*> temp_mass_fit; 
      std::vector<TF1*> temp_beta_fit; 

      for(int b=0; b<fNSlices; b++){
	double p = b*fPWidth + fPMin; 
	int startBin = histos->h2_p_beta[s+1]->GetXaxis()->FindBin(p); 
	int stopBin  = histos->h2_p_beta[s+1]->GetXaxis()->FindBin(p+fPWidth); 
	temp_beta     .push_back(histos->h2_p_beta[s+1]->ProjectionY(Form("h1_beta_pid%d_slice%d_sect%d",histos->GetPid(),b,s), startBin, stopBin));
	temp_deltaBeta.push_back(histos->h2_p_dbeta[s+1]->ProjectionY(Form("h1_dbeta_pid%d_slice%d_sect%d",histos->GetPid(),b,s), startBin, stopBin));
	temp_mass     .push_back(histos->h2_p_tofmass[s+1]->ProjectionY(Form("h1_mass_pid%d_slice%d_sect%d",histos->GetPid(),b,s), startBin, stopBin));
	
	double pionWidth = 0.1 + 0.05*p; 
	double kaonWidth = 0.025;
	double shiftedPionMass = 0.05*p + pid_to_mass(211); 

	// Fitting TOF Mass in momentum bin
	// Cheating to help I do seperate fits first. 
	// I also add a momentum dependent term that 
	// shifts the pion mass up in momentum 
	// based on a simple linear reg. to p=0 
	// and p=1.5 mass values. 
	TF1 *pionMassFit = new TF1("pionFit","gaus", 0.08, 0.17); 
	pionMassFit->SetParameter(1, shiftedPionMass); 
	pionMassFit->SetParameter(2, pionWidth); 
	pionMassFit->SetRange(shiftedPionMass *0.5, shiftedPionMass *1.1); 

	TF1 *kaonMassFit = new TF1("kaonFit","gaus", 0.4, 0.55); 
	kaonMassFit->SetParameter(1, pid_to_mass(321)); 
	kaonMassFit->SetParameter(2, kaonWidth); 

	temp_mass[b]->Fit("pionFit","RQ"); 
	temp_mass[b]->Fit("kaonFit","RQ"); 

	TF1 *temp_mfit = new TF1(Form("f_mass_pid%d_slice%d_sect%d",  histos->GetPid(), b, s),"gaus(0)+gaus(3)", 0.0, 0.75); 
	temp_mfit->SetParameter(0, pionMassFit->GetParameter(0)); 
	temp_mfit->SetParameter(1, pionMassFit->GetParameter(1)); 
	temp_mfit->SetParameter(2, pionMassFit->GetParameter(2)); 
	temp_mfit->SetParameter(3, kaonMassFit->GetParameter(0)); 
	temp_mfit->SetParameter(4, kaonMassFit->GetParameter(1)); 
	temp_mfit->SetParameter(5, kaonMassFit->GetParameter(2)); 
	temp_mass[b]->Fit(Form("f_mass_pid%d_slice%d_sect%d",  histos->GetPid(), b, s),"RQ"); 
	temp_mass_fit.push_back(temp_mfit); 

	// Now doing beta fits 
	// First we need to guess where the peak 
	// should be so that we can have a nice fit 
	double mom      = slices->bins[b]; 
	double pionPeak = mom/sqrt(pow(mom,2)+pow(pid_to_mass(211),2)); 
	double kaonPeak = mom/sqrt(pow(mom,2)+pow(pid_to_mass(321),2)); 

	TF1 *pionBetaFit = new TF1("pionBetaFit","gaus", pionPeak *0.95, pionPeak*1.05); 
	pionMassFit->SetParameter(1, pionPeak); 

	TF1 *kaonBetaFit = new TF1("kaonBetaFit","gaus", kaonPeak *0.96, kaonPeak*1.04); 
	kaonMassFit->SetParameter(1, kaonPeak); 

	temp_beta[b]->Fit("pionBetaFit","RQ"); 
	temp_beta[b]->Fit("kaonBetaFit","RQ"); 

	TF1 *temp_bfit = new TF1(Form("f_beta_pid%d_slice%d_sect%d",  histos->GetPid(), b, s),"gaus(0)+gaus(3)", 0.2, 1.2); 
	temp_bfit->SetParameter(0, pionBetaFit->GetParameter(0)); 
	temp_bfit->SetParameter(1, pionBetaFit->GetParameter(1)); 
	temp_bfit->SetParameter(2, pionBetaFit->GetParameter(2)); 
	temp_bfit->SetParameter(3, kaonBetaFit->GetParameter(0)); 
	temp_bfit->SetParameter(4, kaonBetaFit->GetParameter(1)); 
	temp_bfit->SetParameter(5, kaonBetaFit->GetParameter(2)); 
	//	temp_mass[b]->Fit(Form("f_mass_pid%d_slice%d_sect%d",  histos->GetPid(), b, s),"RQ"); 
	temp_beta_fit.push_back(temp_bfit); 
	
	std::cout << "[MesonFittingService::Slice] Finished momentum slice " << b << " for sector " << s << std::endl; 
      }

      // do the projection over all p 
      TH1D *thisMassProjection = histos->h2_p_tofmass[s+1]->ProjectionY(Form("h1_mass_pid%d_sect%d",histos->GetPid(),s));
      slices->mass_projection.push_back(thisMassProjection); 

      TF1 *pionMassFit = new TF1("pionFit","gaus", 0.08, 0.17); 
      TF1 *kaonMassFit = new TF1("kaonFit","gaus", 0.4, 0.55); 

      pionMassFit->SetParameter(1, pid_to_mass(211)); 
      kaonMassFit->SetParameter(1, pid_to_mass(321)); 

      thisMassProjection->Fit("pionFit","RQ"); 
      thisMassProjection->Fit("kaonFit","RQ"); 

      TF1 *temp_mfit = new TF1(Form("f_mass_pid%d_sect%d",  histos->GetPid(), s),"gaus(0)+gaus(3)", 0.0, 0.75); 
      temp_mfit->SetParameter(0, pionMassFit->GetParameter(0)); 
      temp_mfit->SetParameter(1, pionMassFit->GetParameter(1)); 
      temp_mfit->SetParameter(2, pionMassFit->GetParameter(2)); 
      temp_mfit->SetParameter(3, kaonMassFit->GetParameter(0)); 
      temp_mfit->SetParameter(4, kaonMassFit->GetParameter(1)); 
      temp_mfit->SetParameter(5, kaonMassFit->GetParameter(2)); 
      
      
      slices->mass_projection_fit.push_back(temp_mfit); 
      slices->beta.push_back(temp_beta); 
      slices->deltaBeta.push_back(temp_deltaBeta); 
      slices->mass.push_back(temp_mass); 
      slices->mass_fit.push_back(temp_mass_fit); 
      slices->beta_fit.push_back(temp_beta_fit); 
    }
  }

};

class MesonCutEfficiencyService {
public:
  MesonCutEfficiencyService(MesonSlices *s) : slices(s){
    fIntegrator.SetNumberSteps(10000);

    // integrate up to upper mass cut 
    fIntegrator.SetUpperLimit(0.75); 
  }

  ~MesonCutEfficiencyService(){
  }
  

  void Execute(){
    // If we haven't already 
    SetupCutValues(); 

    // setup the functions 
    for(int s=0; s<slices->mass.size(); s++){

      std::vector<TH1D*> temp_eff; 
      std::vector<TH1D*> temp_cont; 
      std::vector<TH1D*> temp_stat; 

      // keep track of best cut as funct of momentum 
      std::vector<double> bestCutValue; 

      for(int b=0; b<slices->mass[s].size(); b++){
	TF1 *pionFit = new TF1("pionFit", "gaus", slices->mass[s][b]->GetXaxis()->GetBinLowEdge(1), 
			       slices->mass[s][b]->GetXaxis()->GetBinUpEdge( slices->mass[s][b]->GetXaxis()->GetNbins() )); 
	TF1 *kaonFit = new TF1("kaonFit", "gaus", slices->mass[s][b]->GetXaxis()->GetBinLowEdge(1), 
			       slices->mass[s][b]->GetXaxis()->GetBinUpEdge( slices->mass[s][b]->GetXaxis()->GetNbins() )); 
	
	pionFit->SetParameter(0, slices->mass_fit[s][b]->GetParameter(0)); 
	pionFit->SetParameter(1, slices->mass_fit[s][b]->GetParameter(1)); 
	pionFit->SetParameter(2, slices->mass_fit[s][b]->GetParameter(2)); 
	kaonFit->SetParameter(0, slices->mass_fit[s][b]->GetParameter(3)); 
	kaonFit->SetParameter(1, slices->mass_fit[s][b]->GetParameter(4)); 
	kaonFit->SetParameter(2, slices->mass_fit[s][b]->GetParameter(5)); 

	TH1D *h1_eff = new TH1D(Form("h1_eff_slice%d_sect%d",b,s), Form("h1_eff_slice%d_sect%d",b,s), fNumberCutValues, fCutMin, fCutMax); 
	TH1D *h1_cont = new TH1D(Form("h1_cont_slice%d_sect%d",b,s), Form("h1_cont_slice%d_sect%d",b,s), fNumberCutValues, fCutMin, fCutMax); 
	TH1D *h1_stat = new TH1D(Form("h1_stat_slice%d_sect%d",b,s), Form("h1_stat_slice%d_sect%d",b,s), fNumberCutValues, fCutMin, fCutMax); 

	// total kaon count 
	fIntegrator.SetLowerLimit(-0.2);
	fIntegrator.SetFunction(kaonFit);
	double kaonTotal = fIntegrator.Integrate(); 
	
	for(int c=0; c<fCutValues.size(); c++){
	  fIntegrator.SetFunction(pionFit);
	  fIntegrator.SetLowerLimit(fCutValues[c]);
	  double pionVal = fIntegrator.Integrate(); 

	  fIntegrator.SetFunction(kaonFit);
	  double kaonVal = fIntegrator.Integrate(); 

	  double eff = kaonVal/kaonTotal; 
	  double cont = 1-kaonVal/(kaonVal + pionVal); 
	  double stat = eff-cont; 

	  h1_eff->SetBinContent(c+1, eff); 
	  h1_cont->SetBinContent(c+1, cont); 
	  h1_stat->SetBinContent(c+1, stat); 
	}

	temp_eff.push_back(h1_eff); 	
	temp_cont.push_back(h1_cont); 	
	temp_stat.push_back(h1_stat); 	
	
	bestCutValue.push_back(h1_stat->GetBinCenter( h1_stat->GetMaximumBin() )); 

	std::cout << "[MesonCutEfficiencyService] Integrating momentum bin = " << b << " for sector = " << s << std::endl;
      }

      efficiency.push_back(temp_eff); 
      contamination.push_back(temp_cont); 
      test_statistic.push_back(temp_stat); 

      fBestCutGraph[s] = new TGraph(slices->bins.size(), slices->bins.data(), bestCutValue.data());
      fBestCutGraph[s]->SetName(Form("g_bestcut_sect%d",s)); 

      std::string bestCutFitName(Form("f_bestcut_sect%d",s));
      fBestCut[s] = new TF1(bestCutFitName.c_str(),"[0]",0.9,1.5); 
      fBestCutGraph[s]->Fit(bestCutFitName.c_str(),"RQ");
    }
  }

  void DumpPars(std::string conf){
    for(int s=0; s<6; s++){
      std::cout << " conf=" << conf << " sect=" << s << " cut=" << fBestCut[s]->GetParameter(0) << std::endl; 
    }
  }

  void SetCutMin(double min){
    fCutMin = min; 

    fIntegrator.SetLowerLimit(fCutMin); 
  }

  void SetCutMax(double max){
    fCutMax = max; 
  }

  void SetupCutValues(){
    fCutValues.clear();

    for(int i=0; i<fNumberCutValues; i++){
      fCutValues.push_back(i*(fCutMax-fCutMin)/fNumberCutValues +fCutMin);
    }
  }

  void Save(TFile *out){
    out->cd();

    out->mkdir("efficiency"); 
    out->cd("efficiency/"); 

    for(int s=0; s<efficiency.size(); s++){
      for(int b=0; b<efficiency[s].size(); b++){
	efficiency[s][b]->Write(); 
      }
    }

    out->mkdir("contamination"); 
    out->cd("contamination/"); 

    for(int s=0; s<contamination.size(); s++){
      for(int b=0; b<contamination[s].size(); b++){
	contamination[s][b]->Write(); 
      }
    }

    out->mkdir("test_statistic"); 
    out->cd("test_statistic/"); 

    for(int s=0; s<test_statistic.size(); s++){
      for(int b=0; b<test_statistic[s].size(); b++){
	test_statistic[s][b]->Write(); 
      }
    }

    for(int s=0; s<6; s++){
      fBestCutGraph[s]->Write(); 
      fBestCut[s]->Write(); 
    }

    out->cd(); 
  } 

  MesonSlices *slices; 
  std::vector<std::vector<TH1D*> > efficiency; 
  std::vector<std::vector<TH1D*> > contamination; 
  std::vector<std::vector<TH1D*> > test_statistic; 
  TGraph                          *fBestCutGraph[6]; 
  TF1                             *fBestCut[6];

protected:
  const static int fNumberCutValues = 100;
  double fCutMin, fCutMax; 
  std::vector<double>              fCutValues; 
  TF1Integrator_Simpsons           fIntegrator;

};

class MesonCutEfficiencyPlottingService {
public:
  MesonCutEfficiencyPlottingService(MesonHistograms *h, MesonCutEfficiencyService *s) : histos(h), efficiencyService(s) {
    fOutputPath = "unset"; 
  }
  
  void SetOutputPath(std::string p){
    fOutputPath = p; 
  }

  std::string GetOutputPath() const {
    return fOutputPath; 
  }

  void Execute(){

    if(fOutputPath != "unset"){

      gStyle->SetOptTitle(0);
      gStyle->SetOptStat(0);
      
      TCanvas *can = new TCanvas("can","",800, 500); 
      TLatex title; 
      title.SetNDC(); 
      title.SetTextFont(102); 
      title.SetTextSize(0.05); 

      int fillColor = 13; 
      double fillAlpha = 0.4;

      //      Global::Visualization::SetCustomPalette(); 
      //      Global::Visualization::SetBentCoolWarmPalette(); 
      Global::Visualization::SetPurpleYellowPalette(); 

      // draw best cut on 2-d distributions 
      for (int s=0; s<6; s++){
	histos->h2_p_tofmass[s+1]->Draw("colz"); 
	gPad->SetLogz(); 
	gPad->SetGridx(); 
	gPad->SetGridy(); 
	gPad->SetMargin(0.15,0.15, 0.15, 0.15);

	efficiencyService->fBestCut[s]->SetLineStyle(8);
	efficiencyService->fBestCut[s]->SetLineWidth(2);
	efficiencyService->fBestCut[s]->SetRange(0.5, 2.5);
	efficiencyService->fBestCut[s]->Draw("lsame");

	title.DrawLatex(0.44, 0.92, Form("sector %d", s+1)); 
	title.DrawLatex(0.47, 0.08, "p (GeV/c)"); 

	can->Print(Form("%sBestTOFCutSector%d.png",fOutputPath.c_str(),s)); 
      }

      // draw the plots of eff/ cont/ stat 
      for(int s=0; s<efficiencyService->efficiency.size(); s++){
	for(int b=0; b<efficiencyService->efficiency[s].size(); b++){
	  efficiencyService->efficiency[s][b]    ->SetMaximum(1.0);
	  efficiencyService->efficiency[s][b]    ->SetMinimum(0.0);

	  efficiencyService->efficiency[s][b]    ->SetLineColor(99);
	  efficiencyService->contamination[s][b] ->SetLineColor(77);
	  efficiencyService->test_statistic[s][b]->SetLineColor(55);
	  efficiencyService->efficiency[s][b]    ->SetLineWidth(2);
	  efficiencyService->contamination[s][b] ->SetLineWidth(2);
	  efficiencyService->test_statistic[s][b]->SetLineWidth(2);

	  efficiencyService->efficiency[s][b]    ->Draw("l");
	  efficiencyService->contamination[s][b] ->Draw("lsame");
	  efficiencyService->test_statistic[s][b]->Draw("lsame");

	  gPad->SetGridx();
	  gPad->SetGridy();
	  gPad->SetMargin(0.15,0.15, 0.15, 0.15);

	  title.DrawLatex(0.21, 0.92, Form("Efficiency K^{+}, Sector %d, p = %.3f (Gev/c)",s,efficiencyService->slices->bins[b])); 
	  title.DrawLatex(0.21, 0.81, "#color[99]{Efficiency}, #color[77]{Contamination}, #color[55]{Difference}"); 
	  title.DrawLatex(0.45, 0.08, "M_{cut}"); 

	  can->Print(Form("%sEffSlice%dSector%d.png",fOutputPath.c_str(),b,s)); 
	}
      }

    } else {
      std::cerr << "[MesonCutEfficiencyPlottingService] Output path not valid! " << std::endl;
    }
    
    
  }
  
  MesonCutEfficiencyService *efficiencyService;
  MesonHistograms           *histos; 
protected:
  std::string fOutputPath; 

};

class MesonPlottingService {
public:
  MesonPlottingService(MesonHistograms *h, MesonSlices *s) : histos(h), slices(s) { 
    fOutputPath = "unset"; 
  }
  ~MesonPlottingService(){

  }
  
  void SetOutputPath(std::string p){
    fOutputPath = p; 
  }

  std::string GetOutputPath() const {
    return fOutputPath; 
  }


  MesonSlices     *slices; 
  MesonHistograms *histos; 

  void Execute(){
    if(fOutputPath != "unset"){

      gStyle->SetOptTitle(0);
      gStyle->SetOptStat(0);
      
      TCanvas *can = new TCanvas("can","",800, 800); 
      TLatex title; 
      title.SetNDC(); 
      title.SetTextFont(102); 
      title.SetTextSize(0.05); 

      int fillColor = 13; 
      double fillAlpha = 0.4;

      int numberBins = slices->bins.size() +3;
      int numberPages = ceil(numberBins/4); 

      Global::Visualization::SetBentCoolWarmPalette();
      //      Global::Visualization::SetCustomPalette(); 


      for(int s=0; s<6; s++){
	can->Print(Form("%sCompareSlicesSector%d.pdf[", fOutputPath.c_str(), s));

	int bin=0; 
	for(int p=0; p<numberPages; p++){
	  can->Divide(3, 4); 

	  for(int row=1; row<5; row++){
	    int index = 3*row -2; 

	    if (p == 0 && row == 1){
	      can->cd(1); 
	      gPad->SetLogz();
	      gPad->SetMargin(0.2, 0.01, 0.2, 0.01); 
	      histos->h2_p_beta[s+1]->Draw("colz");
	      title.DrawLatex(0.58, 0.08, "p (GeV/c)");

	      can->cd(2); 
	      gPad->SetLogz();
	      gPad->SetMargin(0.2, 0.01, 0.2, 0.01); 
	      histos->h2_p_dbeta[s+1]->Draw("colz");
	      title.DrawLatex(0.58, 0.08, "p (GeV/c)");

	      can->cd(3); 
	      gPad->SetLogz();
	      gPad->SetMargin(0.2, 0.01, 0.2, 0.01); 
	      histos->h2_p_tofmass[s+1]->Draw("colz");
	      title.DrawLatex(0.58, 0.08, "p (GeV/c)");
	    }

	    else {
	      
	      if (bin <numberBins){
		can->cd(index);
		gPad->SetMargin(0.2, 0.01, 0.2, 0.01); 
		gPad->SetLogy(); 
		slices->beta[s][bin]->SetFillColorAlpha(fillColor,fillAlpha);
		slices->beta[s][bin]->SetLineColor(fillColor);
		slices->beta[s][bin]->Draw();
		title.DrawLatex(0.24, 0.93, Form("p = %.3f (GeV/c)",slices->bins[bin])); 
		title.DrawLatex(0.58, 0.08, "#beta");
		
		can->cd(index+1); 
		gPad->SetMargin(0.2, 0.01, 0.2, 0.01); 
		gPad->SetLogy(); 
		slices->deltaBeta[s][bin]->SetFillColorAlpha(fillColor,fillAlpha);
		slices->deltaBeta[s][bin]->SetLineColor(fillColor);
		slices->deltaBeta[s][bin]->Draw();
		title.DrawLatex(0.24, 0.93, Form("p = %.3f (GeV/c)",slices->bins[bin])); 
		title.DrawLatex(0.58, 0.08, "#Delta#beta");
		
		can->cd(index+2); 
		gPad->SetMargin(0.2, 0.01, 0.2, 0.01); 
		//		gPad->SetLogy(); 
		slices->mass[s][bin]->SetFillColorAlpha(fillColor,fillAlpha);
		slices->mass[s][bin]->SetLineColor(fillColor);
		slices->mass[s][bin]->Draw();

		TF1 *pionFit = new TF1("pionFit", "gaus", slices->mass[s][bin]->GetXaxis()->GetBinLowEdge(1), 
				       slices->mass[s][bin]->GetXaxis()->GetBinUpEdge( slices->mass[s][bin]->GetXaxis()->GetNbins() )); 
 		TF1 *kaonFit = new TF1("kaonFit", "gaus", slices->mass[s][bin]->GetXaxis()->GetBinLowEdge(1), 
				       slices->mass[s][bin]->GetXaxis()->GetBinUpEdge( slices->mass[s][bin]->GetXaxis()->GetNbins() )); 
 
		pionFit->SetParameter(0, slices->mass_fit[s][bin]->GetParameter(0)); 
		pionFit->SetParameter(1, slices->mass_fit[s][bin]->GetParameter(1)); 
		pionFit->SetParameter(2, slices->mass_fit[s][bin]->GetParameter(2)); 
		kaonFit->SetParameter(0, slices->mass_fit[s][bin]->GetParameter(3)); 
		kaonFit->SetParameter(1, slices->mass_fit[s][bin]->GetParameter(4)); 
		kaonFit->SetParameter(2, slices->mass_fit[s][bin]->GetParameter(5)); 

		pionFit->SetLineColor(99); 
		kaonFit->SetLineColor(77);
		slices->mass_fit[s][bin]->SetLineColor(55); 

		pionFit->Draw("same"); 
		kaonFit->Draw("same");
		//		slices->mass_fit[s][bin]->Draw("same"); 

		title.DrawLatex(0.24, 0.93, Form("p = %.3f (GeV/c)",slices->bins[bin])); 
		title.DrawLatex(0.58, 0.08, "M_{TOF}");
		
		bin++; 
	      }
	    }
	  }
	  
	  can->Print(Form("%sCompareSlicesSector%d.pdf", fOutputPath.c_str(), s));
	  can->Clear();
	}
	can->Print(Form("%sCompareSlicesSector%d.pdf]", fOutputPath.c_str(), s));
      }

      // print the fit mass distr. 
      TCanvas *sliceCan = new TCanvas("sliceCan", "", 800, 500); 
      
      for(int s=0; s<slices->mass.size(); s++){
	for(int b=0; b<slices->mass[s].size(); b++){
	  sliceCan->Clear(); 

	  slices->mass[s][b]->Draw();
	  gPad->SetGridx(); 
	  gPad->SetGridy(); 
	  gPad->SetMargin(0.15, 0.15, 0.15, 0.15); 
	  //		gPad->SetLogy(); 
	  slices->mass[s][b]->SetFillColorAlpha(fillColor,fillAlpha);
	  slices->mass[s][b]->SetLineColor(fillColor);
	  slices->mass[s][b]->Draw();
	  
	  TF1 *pionFit = new TF1("pionFit", "gaus", slices->mass[s][b]->GetXaxis()->GetBinLowEdge(1), 
				 slices->mass[s][b]->GetXaxis()->GetBinUpEdge( slices->mass[s][b]->GetXaxis()->GetNbins() )); 
	  TF1 *kaonFit = new TF1("kaonFit", "gaus", slices->mass[s][b]->GetXaxis()->GetBinLowEdge(1), 
				 slices->mass[s][b]->GetXaxis()->GetBinUpEdge( slices->mass[s][b]->GetXaxis()->GetNbins() )); 
	  
	  pionFit->SetParameter(0, slices->mass_fit[s][b]->GetParameter(0)); 
	  pionFit->SetParameter(1, slices->mass_fit[s][b]->GetParameter(1)); 
	  pionFit->SetParameter(2, slices->mass_fit[s][b]->GetParameter(2)); 
	  kaonFit->SetParameter(0, slices->mass_fit[s][b]->GetParameter(3)); 
	  kaonFit->SetParameter(1, slices->mass_fit[s][b]->GetParameter(4)); 
	  kaonFit->SetParameter(2, slices->mass_fit[s][b]->GetParameter(5)); 
	  
	  pionFit->SetLineColor(99); 
	  kaonFit->SetLineColor(77);
	  slices->mass_fit[s][b]->SetLineColor(55); 
	  
	  pionFit->Draw("same"); 
	  kaonFit->Draw("same");
	  //		slices->mass_fit[s][bin]->Draw("same"); 
	  
	  title.DrawLatex(0.24, 0.93, Form("p = %.3f (GeV/c)",slices->bins[b])); 
	  title.DrawLatex(0.58, 0.08, "M_{TOF}");
	  
	  sliceCan->Print(Form("%sFitTofMassSlice%dSector%d.png",fOutputPath.c_str(), b, s)); 
	}
      }
      
      
    } else {
      std::cerr << "[MesonPlottingService] Output path is not set. " << std::endl; 
    }
  }

protected: 
  std::string fOutputPath; 

};

int main(int argc, char * argv[]){

    // Setup Options
    h22Options opts;
    opts.args["INPUT"].args = "";
    opts.args["INPUT"].type = 1;
    opts.args["INPUT"].name = "Input file";
    opts.args["PARS"].args = "/u/home/dmriser/Analysis_v2/lists/data.pars";
    opts.args["PARS"].type = 1;
    opts.args["PARS"].name = "Parameter file";
    opts.set(argc,argv);

   
    Parameters *pars = new Parameters(); 
    pars->loadParameters(opts.args["PARS"].args); 
 
    if (opts.args["INPUT"].args != ""){
      
      MesonHistograms pip("pip", 211);
      pip.Load(opts.args["INPUT"].args);
            
      MesonHistograms pim("pim", -211);
      pim.Load(opts.args["INPUT"].args);
      
      MesonHistograms kp("kp", 321);
      kp.Load(opts.args["INPUT"].args);
      
      MesonHistograms km("km", -321);
      km.Load(opts.args["INPUT"].args);

      MesonBetaSliceFitter kaonBetaSlices(&kp, 50, 0.7, 3.0); 
      kaonBetaSlices.Fit("kp"); 
      kaonBetaSlices.WriteParameters(pars); 

      MesonBetaSliceFitter kaonNegBetaSlices(&km, 50, 0.7, 3.0); 
      kaonNegBetaSlices.Fit("km"); 
      kaonNegBetaSlices.WriteParameters(pars); 

      MesonBetaSliceFitter pionBetaSlices(&pip, 50, 0.7, 3.0); 
      pionBetaSlices.Fit("pp"); 
      pionBetaSlices.WriteParameters(pars); 

      MesonBetaSliceFitter pionNegBetaSlices(&pim, 50, 0.7, 3.0); 
      pionNegBetaSlices.Fit("pm"); 
      pionNegBetaSlices.WriteParameters(pars); 

      TFile *out = new TFile(opts.args["OUT"].args.c_str(), "recreate");
      kaonBetaSlices.Save(out); 
      kaonNegBetaSlices.Save(out); 
      kp.Save(out);
      
      
      MesonFittingService fitPip(&pip);
      fitPip.Execute();
      fitPip.Save(out);

      MesonFittingService fitPim(&pim);
      fitPim.Execute();
      fitPim.Save(out);

      MesonFittingService fitKp(&kp);
      fitKp.Execute();
      fitKp.Save(out);

      MesonFittingService fitKm(&km);
      fitKm.Execute();
      fitKm.Save(out);

      MesonCutEfficiencyService effKp(fitKp.slices); 
      effKp.SetCutMin(0.2); 
      effKp.SetCutMax(0.6); 
      effKp.Execute(); 
      effKp.Save(out); 


      MesonCutEfficiencyService effKm(fitKm.slices); 
      effKm.SetCutMin(0.2); 
      effKm.SetCutMax(0.6); 
      effKm.Execute(); 
      effKm.Save(out); 


      MesonCutEfficiencyPlottingService plotKpEff(&kp, &effKp);
      plotKpEff.SetOutputPath("/volatile/clas12/dmriser/plots/pid/kp/");
      plotKpEff.Execute(); 

      MesonCutEfficiencyPlottingService plotKmEff(&km, &effKm);
      plotKmEff.SetOutputPath("/volatile/clas12/dmriser/plots/pid/km/");
      plotKmEff.Execute(); 

      MesonPlottingService plotPip(&pip, fitPip.slices);
      plotPip.SetOutputPath("/volatile/clas12/dmriser/plots/pid/pip/"); 
      plotPip.Execute();

      MesonPlottingService plotPim(&pim, fitPim.slices);
      plotPim.SetOutputPath("/volatile/clas12/dmriser/plots/pid/pim/"); 
      plotPim.Execute();

      MesonPlottingService plotKp(&kp, fitKp.slices);
      plotKp.SetOutputPath("/volatile/clas12/dmriser/plots/pid/kp/"); 
      plotKp.Execute();

      MesonPlottingService plotKm(&km, fitKm.slices);
      plotKm.SetOutputPath("/volatile/clas12/dmriser/plots/pid/km/"); 
      plotKm.Execute();
      
      effKm.DumpPars("km"); 
      effKp.DumpPars("kp"); 
      
      out->Close(); 

    } else {
      std::cerr << "[main] No input file provided with flag -INPUT=file.root" << std::endl;
    }

    pars->saveParameters("meson.pars"); 
    return 0;
}

