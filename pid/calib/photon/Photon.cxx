#ifndef Photon_cxx
#define Photon_cxx

#include "Photon.h"
#include "histos.h"
#include "histos.cxx"

 // Put your includes here 
#include "h22Event.h" 
#include "h22Reader.h" 
#include "Parameters.h"
#include "ParameterSet.h"

#include <iostream>
using std::cout; 
using std::endl; 
using std::flush; 

 // Class Constructor/Destructor 
Photon::Photon(){ 

}

Photon::~Photon(){ 
 
}

void Photon::Loop(){
  Init();
  histos.Init();

  int nev = GetEntries();
  for(int ievent=0; ievent<nev; ievent++){
    GetEntry(ievent); 
    ProcessEvent(); 

    if (ievent%1000 == 0) { cout << "\r done " << ievent << " of " << nev << flush; }
  } 
  cout << endl;

}

void Photon::ProcessEvent(){

  for (int ipart=1; ipart<event.gpart; ipart++){
    if (event.q[ipart] == 0) { histos.Fill(event, ipart); }
  }

}

void Photon::WriteParameters(){
  ParameterSet betaVsPMin;
  betaVsPMin.setName("photonBetaVsPMin");
  betaVsPMin.addValue(0.95);
  betaVsPMin.addError(0.00);
  
  Parameters photonIDParameters; 
  photonIDParameters.addParameterSet(betaVsPMin);
  photonIDParameters.saveParameters("params.dat");
}

#endif
