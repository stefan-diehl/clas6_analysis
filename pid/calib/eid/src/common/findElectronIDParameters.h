#ifndef findElectronIDParameters_h
#define findElectronIDParameters_h


// C++ Includes 
#include <iostream>

 // Put your includes here 
#include "h22Event.h" 
#include "h22Reader.h" 
#include "Parameters.h"
#include "ParameterSet.h"

// CERN Root Includes 
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"

class findElectronIDParameters : public h22Reader {
 public:
  findElectronIDParameters();
  ~findElectronIDParameters();
  
  
  TH2I *ecSamplingFraction[6];
  TH2I *ccThetaCorrelation[6];

  const static int numberSlices = 16;
  TH1D *ecSlices[6][numberSlices];
  TH1D *ccSlices[6][18];

  TGraphErrors *ecMuGraph[6];
  TGraphErrors *ecSigmaGraph[6];
  TGraphErrors *ccMuGraph[6];
  TGraphErrors *ccSigmaGraph[6];
  TF1          *ecPol3Mu[6];
  TF1          *ecPol3Sigma[6];

  Parameters eidParameters; 
  
  void Loop(int numberEvents);
  void ProcessEvent();
  void WriteHardCodedParameters();
  void CalculateSamplingFractionCut();
  void CalculateCCThetaCut();
  void SaveParameters(std::string outputFilename); 
  void SaveHistograms(std::string outputFilename); 
};
#endif
