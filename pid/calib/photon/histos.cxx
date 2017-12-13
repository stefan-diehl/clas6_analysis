#ifndef histo_cxx
#define histo_cxx

#include <iostream>
using std::string; 

#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"

#include "h22Event.h"
#include "histos.h"

void Histograms::Init(){ 
  for (int s=0; s<7; s++){
    betaVsP[s] = new TH2I(Form("betaVsPSector%d",s),"",100,0.05,4.0,100,0.5,1.1);
  }
}

void Histograms::Fill(h22Event event, int ipart){ 
  int sector = event.ec_sect[ipart];

  if (sector > 0 && sector < 7){
    betaVsP[0]->Fill(event.p[ipart],event.b[ipart]);
    betaVsP[sector]->Fill(event.p[ipart],event.b[ipart]);
  }
 
}

void Histograms::Save(string outputFilename){
  TFile *out = new TFile(outputFilename.c_str(),"recreate"); 

  for (int s=0; s<7; s++){
    betaVsP[s]->Write();
  }

  out->Write();
  out->Close();
}

#endif
