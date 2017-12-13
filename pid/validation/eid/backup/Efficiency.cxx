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

// Root 
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"

class Efficiency : public GenericAnalysis {
  
public:
  TH1F *allP, *combinedP;
  vector<TH1F*> effP;

  Parameters        *pars;
  ParticleFilter    *filter;
  DataEventSelector *selector;

  void Initialize();
  void ProcessEvent();
  void Save();

};

void Efficiency::Initialize(){
  pars = new Parameters();
  pars->loadParameters("/u/home/dmriser/mydoc/analysis/root_scripts/Analysis_v2/lists/keppelRadLoose.pars");

  filter = new ParticleFilter(pars);
  filter->set_info(0,1);

  selector = new DataEventSelector();
  selector = filter->GetSelector(11);
  selector->EnableAll();

  int numberPBins = 80; 
  double pMin = 0.7; double pMax = 4.5;

  allP      = new TH1F("pAll","pAll",numberPBins,pMin,pMax);
  combinedP = new TH1F("pCombined","pCombined",numberPBins,pMin,pMax);

  for(int icut=0; icut<selector->cuts.size(); ++icut){
    std::string cutName = selector->cuts[icut]->name();
    TH1F *currentHistoP = new TH1F(Form("pCut%d",icut),cutName.c_str(),numberPBins,pMin,pMax);
    effP.push_back(currentHistoP);
  }

}

void Efficiency::ProcessEvent(){

  for (int ipart=0; ipart<event.gpart; ++ipart){  
    if (event.mcid[ipart] == 11){ allP->Fill(event.p[ipart]); }

    // Check every cut possible 
    for(int icut=0; icut<selector->cuts.size(); ++icut){
      if (selector->cuts[icut]->applies(event, ipart) && selector->cuts[icut]->passes(event, ipart)){
	effP[icut]->Fill(event.p[ipart]);
      }
    }
    
    //    if (selector->passes(event, ipart)){ combinedP->Fill(event.p[ipart]); }
  }

}

void Efficiency::Save(){

  for (int ihist=0; ihist<effP.size(); ++ihist){
    effP[ihist]->Divide(allP); 
  }

  TFile *outputFile = new TFile("EfficiencyLoose.root","recreate"); 
  allP->Write();

  for (int ihist=0; ihist<effP.size(); ++ihist){
    effP[ihist]->Write();
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
