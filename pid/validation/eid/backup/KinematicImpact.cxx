#include <iostream>
#include <vector>
using std::cout; 
using std::endl;
using std::vector;

// My libs
#include "CommonTools.h"
#include "DataEventCut.h"
#include "DataEventSelector.h"
#include "GenericAnalysis.h"
#include "h22Option.h"
#include "Parameters.h"
#include "ParticleFilter.h"
#include "PhysicsEvent.h"
#include "PhysicsEventBuilder.h"

// Root 
#include "TFile.h"
#include "TMath.h"
#include "TH2.h"

class Efficiency : public GenericAnalysis {
  
public:
  TH2F *all, *combined;
  vector<TH2F*> pass;
  vector<TH2F*> fail;

  Parameters          *pars;
  ParticleFilter      *filter;
  DataEventSelector   *selector;
  PhysicsEventBuilder *builder;

  void Initialize();
  void ProcessEvent();
  void Save();
};

void Efficiency::Initialize(){
  pars = new Parameters();
  pars->loadParameters("/u/home/dmriser/mydoc/analysis/root_scripts/Analysis_v2/lists/keppelRad.pars");

  filter = new ParticleFilter(pars);
  filter->set_info(0,1);

  selector = new DataEventSelector();
  selector = filter->GetSelector(11);
  selector->EnableAll();

  all      = new TH2F("all","all",200,0.8,3.1,200,0.5,5.0);
  combined = new TH2F("combined","combined",200,0.8,3.1,200,0.5,5.0);

  for(int icut=0; icut<selector->cuts.size(); ++icut){
    std::string cutName = selector->cuts[icut]->name();
    TH2F *currentPassHisto = new TH2F(Form("pCut%d",icut),cutName.c_str(),200,0.8,3.1,200,0.5,5.0);
    TH2F *currentFailHisto = new TH2F(Form("fCut%d",icut),cutName.c_str(),200,0.8,3.1,200,0.5,5.0);
    pass.push_back(currentPassHisto);
    fail.push_back(currentFailHisto);
  }

  builder = new PhysicsEventBuilder();
}

void Efficiency::ProcessEvent(){

  for (int ipart=0; ipart<event.gpart; ++ipart){  
    if (event.q[ipart] == -1){ 
      PhysicsEvent physicsEvent = builder->getPhysicsEvent(event.GetTLorentzVector(ipart, 11)); 
      all->Fill(physicsEvent.w, physicsEvent.qq); 
      
      // Check every cut possible 
      for(int icut=0; icut<selector->cuts.size(); ++icut){
	if (selector->cuts[icut]->applies(event, ipart) && selector->cuts[icut]->passes(event, ipart)){
	  pass[icut]->Fill(physicsEvent.w, physicsEvent.qq);
	} 
	else {
	  fail[icut]->Fill(physicsEvent.w, physicsEvent.qq);
	}
      }

      if (selector->IsPassed(event, ipart)){ combined->Fill(physicsEvent.w, physicsEvent.qq); }
    }
  }

}

void Efficiency::Save(){

  //  for (int ihist=0; ihist<pass.size(); ++ihist){
  //    pass[ihist]->Divide(all); 
  //  }

  TFile *outputFile = new TFile("KinematicImpactMC.root","recreate"); 
  all     ->Write();
  combined->Write();

  for (int ihist=0; ihist<pass.size(); ++ihist){
    pass[ihist]->Write();
    fail[ihist]->Write();
  }
  //  outputFile->Write(); 
  outputFile->Close(); 
}

int main(int argc, char *argv[]){


  h22Options opts;
  opts.set(argc, argv);

  Efficiency Eff;

  for (int ifile=0; ifile<opts.ifiles.size(); ifile++){ Eff.AddFile(opts.ifiles[ifile]); }
  Eff.RunAnalysis(opts.args["N"].arg); 

  return 0;
}
