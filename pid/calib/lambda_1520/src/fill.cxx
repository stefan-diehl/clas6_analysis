#include <iostream>

// h22 libs 
#include "CommonTools.h"
#include "Corrections.h"
#include "CheckPoints.h"
#include "DBins.h"
#include "h22Event.h"
#include "h22Reader.h"
#include "MesonHistograms.h"
#include "MomCorr.h"
#include "ParticleFilter.h"
#include "Parameters.h"
#include "PhysicsEvent.h"
#include "PhysicsEventBuilder.h"
#include "StatusBar.h"
#include "StandardHistograms.h"

// root 
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TVector3.h"

// this project 
#include "common/TreeOutput.h"

const static float M_SIGMA_PLUS  = 1.189;
const static float M_SIGMA_MINUS = 1.197;

class Analysis : public h22Reader {
public:  
  Analysis() {

    std::string path        = Global::Environment::GetAnalysisPath(); 
    std::string momCorrPath = Form("%s/momCorr/",path.c_str());
    momCorr                 = new MomCorr_e1f(momCorrPath);

    // setup reader options 
    GSIM = false; 
    Init();

    // setup tree 
    treeOutput = new TreeOutput(); 
    exclusiveTree = new TreeOutput();

    // needs parameters 
    params = new Parameters(); 
    params->loadParameters(Form("%s/lists/data_tofmass.pars", path.c_str())); 
    filter = new ParticleFilter(params);

    // histos 
    std::string type[2] = {"signal", "background"};

    mm_ek[0] = new TH1I("mm_ek_no_cut", "", 100, 1.1, 1.8);
    mm_ek[1] = new TH1I("mm_ek_prot", "", 100, 1.1, 1.8);
    mm_ek[2] = new TH1I("mm_ek_pip", "", 100, 1.1, 1.8);
    mm_ek[3] = new TH1I("mm_ek_pim", "", 100, 1.1, 1.8);

    mm_ekpion[0] = new TH1I("mm_ekpion_plus","",100,0.9,1.6);
    mm_ekpion[1] = new TH1I("mm_ekpion_minus","",100,0.9,1.6);

    mm_ekproton  = new TH1I("mm_ekproton","",100,0.0,0.9);

    for (int t=0; t<2; t++){
      theta[t] = new TH1I(Form("theta_%s",type[t].c_str()),"",100,-180,180); 
      theta_pt[t] = new TH2I(Form("theta_pt_%s",type[t].c_str()),"",100,-180,180,100,0,1.5);
      
      theta_kk[t]  = new TH1I(Form("theta_kk_%s",type[t].c_str()), "", 100, -10, 70);
      theta_kpp[t] = new TH1I(Form("theta_kpp_%s",type[t].c_str()), "", 100, -10, 70);
      theta_kmp[t] = new TH1I(Form("theta_kmp_%s",type[t].c_str()), "", 100, -10, 70);
      theta_kpp_kmp[t] = new TH2I(Form("theta_kpp_kmp_%s",type[t].c_str()), "", 100, -10, 70, 100, -10, 70);
      
      p_prot[t] = new TH1I(Form("p_prot_%s",type[t].c_str()), "", 100, 0.3, 3.5);
      p_kp[t] = new TH1I(Form("p_kp_%s",type[t].c_str()), "", 100, 0.3, 3.5);
      p_km[t] = new TH1I(Form("p_km_%s",type[t].c_str()), "", 100, 0.3, 3.5);
      
      p_theta_prot[t] = new TH2I(Form("p_theta_prot_%s",type[t].c_str()), "", 100, 0.3, 3.5, 100, -10, 70); 
      p_theta_kp[t]   = new TH2I(Form("p_theta_kp_%s",type[t].c_str()), "", 100, 0.3, 3.5, 100, -10, 70); 
      p_theta_km[t]   = new TH2I(Form("p_theta_km_%s",type[t].c_str()), "", 100, 0.3, 3.5, 100, -10, 70); 
    }

  }

  ~Analysis(){
  }

  void Loop(){
    
    // setup particle filter 
    filter->set_info(GSIM, GetRunNumber());

    StatusBar stat; 

    TStopwatch timer; 
    timer.Reset();
    timer.Start();

    // for every event
    for(int ientry=0; ientry<GetEntries(); ientry++){
      GetEntry(ientry);
      
      std::vector<int> electronIndices = filter->getVectorOfParticleIndices(event, 11);
      if(!electronIndices.empty()){
	
	int electronIndex = electronIndices[0];
	event.SetElectronIndex(electronIndex);
	Corrections::correctEvent(&event, GetRunNumber(), GSIM);
	
	std::vector<int> mesonIndices = filter->getVectorOfParticleIndices(event, 321);
	if (!mesonIndices.empty()){
	  
	  // build an event 
	  int mesonIndex = mesonIndices[0];
	  TLorentzVector electron = event.GetTLorentzVector(electronIndex, 11);  
	  TLorentzVector meson    = event.GetTLorentzVector(mesonIndex, 321);
	  
	  electron = momCorr->PcorN(electron, -1, 11);
	  PhysicsEvent ev = builder.getPhysicsEvent(electron,meson); 
	  
	  
	  // Dan Carman method for rejecting some pi+ events in the 
	  // kaon sample. 
	  //
	  // Classify an event and if it's a kaon you pretend 
	  // it's a pion.  Then check the missing mass spectrum 
	  // and if the event is near the neutron peak it must 
	  // really be a pion, so reassign it to be.  
	  //
	  // Currently all it does it throw out the pion, 
	  // if doesn't actually add it to the pion sample 
	  // because the code never makes it that far in the 
	  // pion case it would just return 0 pions and not enter.
	  TLorentzVector pion = event.GetTLorentzVector(mesonIndex, 211); 
	  PhysicsEvent pionEvent = builder.getPhysicsEvent(electron, pion); 
	  
	  // reject and dont process event 
	  if (sqrt(pionEvent.mm2) < 1.07){
	    continue; 
	  }

	  // fill the normal one 
	  mm_ek[0]->Fill(sqrt(ev.mm2));

	  // look for things 
	  std::vector<int> protonIndices = filter->getVectorOfParticleIndices(event, 2212);
	  std::vector<int> pipIndices    = filter->getVectorOfParticleIndices(event, 211);
	  std::vector<int> pimIndices    = filter->getVectorOfParticleIndices(event,-211);

	  if (!protonIndices.empty()){
	    TLorentzVector proton    = event.GetTLorentzVector(protonIndices[0], 2212);
	    PhysicsEvent protonEvent = builder.getPhysicsEvent(electron,meson,proton); 

	    mm_ekproton->Fill(sqrt(protonEvent.mm2));

	    // fill events with km 
	    if(sqrt(protonEvent.mm2) > 0.40 && 
	       sqrt(protonEvent.mm2) < 0.60){
	      mm_ek[1]->Fill(sqrt(ev.mm2));


	      treeOutput->x  = ev.x; 
	      treeOutput->w  = ev.w; 
	      treeOutput->q2 = ev.qq; 
	      treeOutput->z  = ev.z; 
	      treeOutput->pt = ev.pT; 
	      treeOutput->missing_mass = sqrt(ev.mm2);
	      treeOutput->hel = event.corr_hel; 

	      treeOutput->p_prot = proton.P();
	      treeOutput->p_kp = meson.P();
	      treeOutput->p_km = protonEvent.finalState.P();

	      treeOutput->angle_kp_km = to_degrees*meson.Angle(protonEvent.finalState.Vect());
	      treeOutput->angle_kp_prot = to_degrees*proton.Angle(meson.Vect());
	      treeOutput->angle_km_prot = to_degrees*protonEvent.finalState.Angle(proton.Vect());

	      treeOutput->phi_kp = ev.phiHadron; 

	      // we need a copy of the event with only the 
	      // kaon 
	      PhysicsEvent kaonOnlyEvent = builder.getPhysicsEvent(electron, meson);
	      treeOutput->angle_kp_lambda = to_degrees*kaonOnlyEvent.finalState.Angle(meson.Vect());
	      treeOutput->angle_km_lambda = to_degrees*kaonOnlyEvent.finalState.Angle(protonEvent.finalState.Vect());
	      treeOutput->angle_prot_lambda = to_degrees*kaonOnlyEvent.finalState.Angle(proton.Vect());	      
	      treeOutput->p_lambda = kaonOnlyEvent.finalState.P();
	      
	      TVector3 boostToCM = -(ev.virtualPhoton+ev.targetParticle).BoostVector();
	      TLorentzVector boostedKp = meson; 
	      TLorentzVector boostedKm = ev.finalState; 
	      TLorentzVector boostedLam = kaonOnlyEvent.finalState; 
	      boostedKp.Boost(boostToCM);
	      boostedKm.Boost(boostToCM);
	      boostedLam.Boost(boostToCM);
	      
	      treeOutput->angle_kp_lambda_cm = to_degrees*boostedKp.Angle(boostedLam.Vect());

	      treeOutput->theta_kp = boostedKp.Theta()*to_degrees; 
	      treeOutput->theta_km = boostedKm.Theta()*to_degrees; 

	      // fill tree 
	      treeOutput->tree->Fill();


	      // look for the real thing 
	      std::vector<int> kms = filter->getVectorOfParticleIndices(event, -321);
	      if(!kms.empty()){
		int kmIndex = kms[0]; 
		TLorentzVector km = event.GetTLorentzVector(kms[0], -321);

		PhysicsEvent exEv = builder.getPhysicsEvent(electron, meson, proton, km);
		
		exclusiveTree->x  = exEv.x; 
		exclusiveTree->w  = exEv.w; 
		exclusiveTree->q2 = exEv.qq; 
		exclusiveTree->z  = ev.z; 
		exclusiveTree->pt = ev.pT; 
		exclusiveTree->missing_mass = sqrt(ev.mm2);
		exclusiveTree->hel = event.corr_hel; 

		exclusiveTree->phi_kp = ev.phiHadron; 		
	      
		exclusiveTree->p_prot = proton.P();
		exclusiveTree->p_kp = meson.P();
		exclusiveTree->p_km = km.P();
		
		exclusiveTree->angle_kp_km   = to_degrees*meson.Angle(km.Vect());
		exclusiveTree->angle_kp_prot = to_degrees*proton.Angle(meson.Vect());
		exclusiveTree->angle_km_prot = to_degrees*km.Angle(proton.Vect());

		// invariant masses 
		exclusiveTree->im_km_prot = sqrt((km+proton).M2());
		exclusiveTree->im_km_kp   = sqrt((km+meson).M2());
		
		// we need a copy of the event with only the 
		// kaon 
		TLorentzVector lambda = km+proton; 

		exclusiveTree->angle_kp_lambda = to_degrees*lambda.Angle(meson.Vect());
		exclusiveTree->angle_km_lambda = to_degrees*lambda.Angle(protonEvent.finalState.Vect());
		exclusiveTree->angle_prot_lambda = to_degrees*lambda.Angle(proton.Vect());	      
		exclusiveTree->p_lambda = lambda.P();

		TVector3 boost = -(exEv.virtualPhoton + exEv.targetParticle).BoostVector();
		lambda.Boost(boost);
		meson.Boost(boost);
		km.Boost(boost);
		exclusiveTree->angle_kp_lambda_cm = to_degrees*meson.Angle(lambda.Vect());
		
		// boost to vf frame and take angles 
		exclusiveTree->theta_kp = meson.Theta()*to_degrees; 
		exclusiveTree->theta_km = km.Theta()*to_degrees; 

		// fill tree 
		exclusiveTree->tree->Fill();		
	      }

	      if(sqrt(ev.mm2) > 1.50 && sqrt(ev.mm2) < 1.54){
		theta[0]->Fill(ev.thetaHadron);
		theta_pt[0]->Fill(ev.thetaHadron, ev.pT);
		theta_kpp[0]->Fill(to_degrees*proton.Angle(meson.Vect()));
		theta_kk[0]->Fill(to_degrees*meson.Angle(ev.finalState.Vect()));
		theta_kmp[0]->Fill(to_degrees*proton.Angle(ev.finalState.Vect()));
 		theta_kpp_kmp[0]->Fill(to_degrees*proton.Angle(meson.Vect()), to_degrees*proton.Angle(ev.finalState.Vect()));
		p_prot[0]->Fill(proton.P()); 
		p_kp[0]  ->Fill(meson.P()); 
		p_km[0]  ->Fill(ev.finalState.P());
		p_theta_prot[0]->Fill(proton.P(), proton.Theta()*to_degrees); 
		p_theta_kp[0]  ->Fill(meson.P(), meson.Theta()*to_degrees); 
		p_theta_km[0]  ->Fill(ev.finalState.P(), ev.finalState.Theta()*to_degrees);
	      }
	      else if(sqrt(ev.mm2) > 1.40 && sqrt(ev.mm2) < 1.50 ||
		      sqrt(ev.mm2) > 1.54 && sqrt(ev.mm2) < 1.70){
		theta[1]->Fill(ev.thetaHadron);
		theta_pt[1]->Fill(ev.thetaHadron, ev.pT);
		theta_kpp[1]->Fill(to_degrees*proton.Angle(meson.Vect()));
		theta_kk[1]->Fill(to_degrees*meson.Angle(ev.finalState.Vect()));
		theta_kmp[1]->Fill(to_degrees*proton.Angle(ev.finalState.Vect()));
		theta_kpp_kmp[1]->Fill(to_degrees*proton.Angle(meson.Vect()), to_degrees*proton.Angle(ev.finalState.Vect()));
		p_prot[1]->Fill(proton.P()); 
		p_kp[1]  ->Fill(meson.P()); 
		p_km[1]  ->Fill(ev.finalState.P());
		p_theta_prot[1]->Fill(proton.P(), proton.Theta()*to_degrees); 
		p_theta_kp[1]  ->Fill(meson.P(), meson.Theta()*to_degrees); 
		p_theta_km[1]  ->Fill(ev.finalState.P(), ev.finalState.Theta()*to_degrees);
	      }

	    }

	  }

	  if(!pipIndices.empty()){
	    TLorentzVector pion    = event.GetTLorentzVector(pipIndices[0], 211);
	    PhysicsEvent pionEvent = builder.getPhysicsEvent(electron,meson,pion); 

	    mm_ekpion[0]->Fill(sqrt(pionEvent.mm2));

	    if(sqrt(pionEvent.mm2) > 1.13 && 
	       sqrt(pionEvent.mm2) < 1.26){
	      mm_ek[2]->Fill(sqrt(ev.mm2));
	    }
	  }

	  if(!pimIndices.empty()){
	    TLorentzVector pion    = event.GetTLorentzVector(pimIndices[0], -211);
	    PhysicsEvent pionEvent = builder.getPhysicsEvent(electron,meson,pion); 

	    mm_ekpion[1]->Fill(sqrt(pionEvent.mm2));

	    if(sqrt(pionEvent.mm2) > 1.13 && 
	       sqrt(pionEvent.mm2) < 1.26){
	      mm_ek[3]->Fill(sqrt(ev.mm2));
	    }
	  }

	}
      }
      

      if (ientry%10000 == 0){
	stat.PrintStatus(ientry, GetEntries()); 
      }
    }

    std::cout << std::endl;
    double loopTime  = timer.RealTime(); 
    double eventRate = GetEntries()/loopTime; 
    std::cout << "[GenericAnalysis::Loop] Total time: " 
	      << loopTime << " seconds, Event rate: " 
	      << eventRate << " events/sec. " << std::endl;

  }

  void Save(std::string outFile){
    TFile *out = new TFile(outFile.c_str(), "recreate");

    mm_ek[0]->Write(); 
    mm_ek[1]->Write(); 
    mm_ek[2]->Write(); 
    mm_ek[3]->Write(); 
    mm_ekpion[0]->Write(); 
    mm_ekpion[1]->Write(); 
    mm_ekproton->Write(); 

    for(int t=0; t<2; t++){
      theta[t]    ->Write(); 
      theta_pt[t] ->Write();
      theta_kk[t] ->Write(); 
      theta_kpp[t]->Write(); 
      theta_kmp[t]->Write(); 
      theta_kpp_kmp[t]->Write();
      
      p_prot[t]->Write();
      p_kp[t]->Write();
      p_km[t]->Write();
      p_theta_kp[t]->Write();
      p_theta_km[t]->Write();
      p_theta_prot[t]->Write();

    }

    treeOutput->Save(out); 

    out->Close();

    TFile *exOut = new TFile(Form("ex_%s", outFile.c_str()), "recreate"); 
    exclusiveTree->Save(exOut);
    exOut->Close();
  }


protected:
  MomCorr_e1f             *momCorr; 
  PhysicsEventBuilder      builder; 
  ParticleFilter          *filter; 
  Parameters              *params; 

  // histograms 
  // 0 - no cuts 
  // 1 - proton in event and cut on k- mass
  // 2 - pi+ in event and cut on sigma mass 
  // 3 - pi- in event and cut on sigma mass 
  TH1I *mm_ek[4];

  // centered at sigma mass 
  // 0 - pi+ 
  // 1 - pi- 
  TH1I *mm_ekpion[2];

  // centered at kaon mass 
  // 0 - proton 
  TH1I *mm_ekproton;

  // plot theta angle for events 
  TH1I *theta[2];
  TH2I *theta_pt[2]; 

  // angle between particles 
  TH1I *theta_kk[2]; 
  TH1I *theta_kpp[2];  
  TH1I *theta_kmp[2];  
  TH2I *theta_kpp_kmp[2];   

  // momenta 
  TH1I *p_prot[2]; 
  TH1I *p_kp[2]; 
  TH1I *p_km[2]; 
  TH2I *p_theta_prot[2]; 
  TH2I *p_theta_kp[2]; 
  TH2I *p_theta_km[2]; 

  // tree 
  TreeOutput *treeOutput; 
  TreeOutput *exclusiveTree; 
};

int main(int argc, char *argv[]){

  if (argc > 1){
  // setup analysis 
  Analysis analysis; 

  for (int iarg=1; iarg<argc; iarg++){
    std::string f(argv[iarg]);
    analysis.AddFile(f); 
  }

  // run analysis loop
  analysis.Loop();
  analysis.Save("out.root");
  
  } else {
    std::cerr << "No files found." << std::endl; 
  }

  return 0; 
}
