#ifndef findElectronIDParameters_cxx
#define findElectronIDParameters_cxx

#include "findElectronIDParameters.h"

#include <iostream>
#include <vector>

// Put your includes here 
#include "h22Event.h" 
#include "h22Reader.h" 
#include "Parameters.h"
#include "ParameterSet.h"

// CERN Root Includes
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"

// Class Constructor/Destructor 
findElectronIDParameters::findElectronIDParameters(){   
  for (int sect=0; sect<6; sect++){
    ecSamplingFraction[sect] = new TH2I(Form("ecSamplingFraction_%d",sect),"",100,0.5,4.5,100,0.01,0.5);
    ccThetaCorrelation[sect] = new TH2I(Form("ccThetaCorrelation_%d",sect),"",18,0,17,100,0,80);
  }
}

findElectronIDParameters::~findElectronIDParameters(){ 

}

void findElectronIDParameters::Loop(int numberEvents){

  // Event loop below. 
  int nen = GetEntries(); 
  if (numberEvents < nen){ nen = numberEvents; }

  for(int ievent=0; ievent<nen; ievent++){
    GetEntry(ievent); 
    ProcessEvent(); 
    if (ievent%1000 == 0) { std::cout << "\r done " << ievent << " of " << nen << std::flush; }
  } 

  std::cout << std::endl;
}

void findElectronIDParameters::ProcessEvent(){
  
  for (int ipart =0; ipart<event.gpart; ipart++){
    if (event.q[ipart]<0 && event.ec_ei[ipart]>0.05 && event.nphe[ipart]>20) {
      int sectorIndex = event.dc_sect[ipart] -1;
      
      if(sectorIndex > -1 && sectorIndex < 6){
	ecSamplingFraction[sectorIndex]->Fill(event.p[ipart], event.etot[ipart]/event.p[ipart]);
	ccThetaCorrelation[sectorIndex]->Fill((event.cc_segm[ipart]%1000)/10, event.GetThetaCC(ipart));
      }
    }
  }

}

void findElectronIDParameters::CalculateCCThetaCut(){

  // Routine is supposed to cut slices out of the etot/p vs. p 
  // then fit with gaussians in each bin and further fit the mean 
  // with pol3. 
  int binsPerSlice = 1;

  TF1 *fitGauss = new TF1("fitGauss","gaus");  
  TF1 *fitPol2  = new TF1("fitPol2","pol2", 4.0,12.0);

  // Used for TGraphErrors constructor
  double dummyAxis[numberSlices];
  double dummyAxisError[numberSlices];
  for(int d=0; d<18; d++){ dummyAxis[d] = (double) d; dummyAxisError[d] = 0.00; } 

  double mu[numberSlices];
  double sigma[numberSlices];

  double muError[numberSlices];
  double sigmaError[numberSlices];

  ParameterSet muParA, muParB, muParC, muParD, 
    sigmaParA, sigmaParB, sigmaParC, sigmaParD, nSigma;

  nSigma   .setName("EL_CCT_NSIGMA"); 
  muParA   .setName("EL_CCT_MU_A");
  muParB   .setName("EL_CCT_MU_B");
  muParC   .setName("EL_CCT_MU_C");
  muParD   .setName("EL_CCT_MU_D");
  sigmaParA.setName("EL_CCT_SIGMA_A");
  sigmaParB.setName("EL_CCT_SIGMA_B");
  sigmaParC.setName("EL_CCT_SIGMA_C");
  sigmaParD.setName("EL_CCT_SIGMA_D");


  for (int sector=0; sector<6; sector++){
    for(int slice=0; slice<18; slice++){
      std::string name = Form("ccSlices_%d_sector_%d",slice,sector);
      
      fitGauss->SetRange(-2.0 + 2.5*slice, 2.0 + 2.5*slice); 

      ccSlices[sector][slice] = ccThetaCorrelation[sector]->ProjectionY(name.c_str(),slice,slice+1);
      ccSlices[sector][slice]->Fit("fitGauss","RQ");
 
      mu[slice]         = fitGauss->GetParameter(1);
      sigma[slice]      = fitGauss->GetParameter(2);
      muError[slice]    = fitGauss->GetParError(1);
      sigmaError[slice] = fitGauss->GetParError(2);
    }
    
    ccMuGraph[sector]    = new TGraphErrors(18, dummyAxis,    mu, dummyAxisError, muError);
    ccSigmaGraph[sector] = new TGraphErrors(18, dummyAxis, sigma, dummyAxisError, sigmaError);

    ccMuGraph[sector]   ->SetName(Form("ccMuGraph%d",sector));
    ccSigmaGraph[sector]->SetName(Form("ccSigmaGraph%d",sector));

    ccMuGraph[sector]->Fit("fitPol2","RQ");
    muParA.addValue(0.0);
    muParA.addError(0.0);
    //    muParA.addValue(fitPol3->GetParameter(3));
    //    muParA.addError(fitPol3->GetParError(3));
    muParB.addValue(fitPol2->GetParameter(2));
    muParB.addError(fitPol2->GetParError(2));    
    muParC.addValue(fitPol2->GetParameter(1));
    muParC.addError(fitPol2->GetParError(1));    
    muParD.addValue(fitPol2->GetParameter(0));
    muParD.addError(fitPol2->GetParError(0));

    //    ccSigmaGraph[sector]->Fit("fitPol3","RQ");
    ccSigmaGraph[sector]->Fit("fitPol2","RQ");
    //    sigmaParA.addValue(fitPol3->GetParameter(3));
    //    sigmaParA.addError(fitPol3->GetParError(3));
    sigmaParA.addValue(0.0);
    sigmaParA.addError(0.0);
    sigmaParB.addValue(fitPol2->GetParameter(2));
    sigmaParB.addError(fitPol2->GetParError(2));    
    sigmaParC.addValue(fitPol2->GetParameter(1));
    sigmaParC.addError(fitPol2->GetParError(1));    
    sigmaParD.addValue(fitPol2->GetParameter(0));
    sigmaParD.addError(fitPol2->GetParError(0));
  }


  nSigma.addValueAndError(2.00, 0.00);

  // Save those values
  eidParameters.addParameterSet(nSigma);
  eidParameters.addParameterSet(muParA);
  eidParameters.addParameterSet(muParB);
  eidParameters.addParameterSet(muParC);
  eidParameters.addParameterSet(muParD);
  eidParameters.addParameterSet(sigmaParA);
  eidParameters.addParameterSet(sigmaParB);
  eidParameters.addParameterSet(sigmaParC);
  eidParameters.addParameterSet(sigmaParD);
}

void findElectronIDParameters::CalculateSamplingFractionCut(){

  // Routine is supposed to cut slices out of the etot/p vs. p 
  // then fit with gaussians in each bin and further fit the mean 
  // with pol3. 
  int binsPerSlice = (int) ecSamplingFraction[0]->GetXaxis()->GetNbins()/numberSlices; 
  double binWidth  = ecSamplingFraction[0]->GetXaxis()->GetBinWidth(0)*binsPerSlice;
  double pMin      = ecSamplingFraction[0]->GetXaxis()->GetBinLowEdge(0);

  // For data use 0.25-0.45
  double yMin = 0.25;
  double yMax = 0.45;

  //  TF1 *fitGauss = new TF1("fitGauss","gaus",0.1,0.45);  
  TF1 *fitGauss = new TF1("fitGauss","gaus",yMin,yMax);  
  //  TF1 *fitPol3  = new TF1("fitPol3","pol3", 0.50,3.5);


  // Used for TGraphErrors constructor
  double dummyAxis[numberSlices];
  double dummyAxisError[numberSlices];
  for(int d=0; d<numberSlices; d++){ dummyAxis[d] = (double) pMin+d*binWidth; dummyAxisError[d] = 0.00; } 

  double mu[numberSlices];
  double sigma[numberSlices];

  double muError[numberSlices];
  double sigmaError[numberSlices];

  ParameterSet muParA, muParB, muParC, muParD, 
    sigmaParA, sigmaParB, sigmaParC, sigmaParD;

  muParA   .setName("EL_SF_MU_A");
  muParB   .setName("EL_SF_MU_B");
  muParC   .setName("EL_SF_MU_C");
  muParD   .setName("EL_SF_MU_D");
  sigmaParA.setName("EL_SF_SIGMA_A");
  sigmaParB.setName("EL_SF_SIGMA_B");
  sigmaParC.setName("EL_SF_SIGMA_C");
  sigmaParD.setName("EL_SF_SIGMA_D");


  for (int sector=0; sector<6; sector++){
    for(int slice=0; slice<numberSlices; slice++){

      std::string name = Form("ecSlices_%d_sector_%d",slice,sector);
      std::cout << "[ECSamplingCuts] Projecting Bin=" << slice*binsPerSlice << " to Bin=" << slice*binsPerSlice + binsPerSlice << std::endl;
      ecSlices[sector][slice] = ecSamplingFraction[sector]->ProjectionY(name.c_str(),slice*binsPerSlice,binsPerSlice+slice*binsPerSlice);
      ecSlices[sector][slice]->Fit("fitGauss","RQ");
 
      mu[slice]         = fitGauss->GetParameter(1);
      sigma[slice]      = fitGauss->GetParameter(2);
      muError[slice]    = fitGauss->GetParError(1);
      sigmaError[slice] = fitGauss->GetParError(2);
    }
    
    ecMuGraph[sector]    = new TGraphErrors(numberSlices, dummyAxis,    mu, dummyAxisError, muError);
    ecSigmaGraph[sector] = new TGraphErrors(numberSlices, dummyAxis, sigma, dummyAxisError, sigmaError);
    ecPol3Mu[sector]     = new TF1(Form("ecPol3Mu%d",sector),   "pol3",0.50,4.5);
    ecPol3Sigma[sector]  = new TF1(Form("ecPol3Sigma%d",sector),"pol3",0.50,4.5);

    ecMuGraph[sector]   ->SetName(Form("ecMuGraph%d",sector));
    ecSigmaGraph[sector]->SetName(Form("ecSigmaGraph%d",sector));

    ecMuGraph[sector]->Fit(Form("ecPol3Mu%d",sector),"RQ");
    muParA.addValue(ecPol3Mu[sector]->GetParameter(3));
    muParA.addError(ecPol3Mu[sector]->GetParError(3));
    muParB.addValue(ecPol3Mu[sector]->GetParameter(2));
    muParB.addError(ecPol3Mu[sector]->GetParError(2));    
    muParC.addValue(ecPol3Mu[sector]->GetParameter(1));
    muParC.addError(ecPol3Mu[sector]->GetParError(1));    
    muParD.addValue(ecPol3Mu[sector]->GetParameter(0));
    muParD.addError(ecPol3Mu[sector]->GetParError(0));

    ecSigmaGraph[sector]->Fit(Form("ecPol3Sigma%d",sector),"RQ");
    sigmaParA.addValue(ecPol3Sigma[sector]->GetParameter(3));
    sigmaParA.addError(ecPol3Sigma[sector]->GetParError(3));
    sigmaParB.addValue(ecPol3Sigma[sector]->GetParameter(2));
    sigmaParB.addError(ecPol3Sigma[sector]->GetParError(2));    
    sigmaParC.addValue(ecPol3Sigma[sector]->GetParameter(1));
    sigmaParC.addError(ecPol3Sigma[sector]->GetParError(1));    
    sigmaParD.addValue(ecPol3Sigma[sector]->GetParameter(0));
    sigmaParD.addError(ecPol3Sigma[sector]->GetParError(0));
  }


  // Save those values
  eidParameters.addParameterSet(muParA);
  eidParameters.addParameterSet(muParB);
  eidParameters.addParameterSet(muParC);
  eidParameters.addParameterSet(muParD);
  eidParameters.addParameterSet(sigmaParA);
  eidParameters.addParameterSet(sigmaParB);
  eidParameters.addParameterSet(sigmaParC);
  eidParameters.addParameterSet(sigmaParD);

}

void findElectronIDParameters::WriteHardCodedParameters(){
  ParameterSet zVertexLimitMin, zVertexLimitMax, ecUCoordMin, ecVCoordMin, ecWCoordMin,
    ecUCoordMax, ecVCoordMax, ecWCoordMax, ecEdepMin, ecNSigma, dcr1FidH, dcr1FidA, 
    dcr3FidH, dcr3FidA, ccFidA, ccFidB, ccFidC, NPheMin, NPheMax; 

  NPheMin.setName("EL_NPHE_MIN");
  NPheMin.addValueAndError(25.0, 0.0);

  NPheMax.setName("EL_NPHE_MAX");
  NPheMax.addValueAndError(999.0, 0.0);

  ecNSigma.setName("EL_EC_NSIGMA");
  ecNSigma.addValue(1.00);
  ecNSigma.addError(0.00);

  dcr1FidH.setName("EL_DCR1_FIDH");
  dcr1FidH.addValue(22.0);
  dcr1FidH.addError(0.0);

  dcr1FidA.setName("EL_DCR1_FIDA");
  dcr1FidA.addValue(60.0);
  dcr1FidA.addError(0.0);

  dcr3FidH.setName("EL_DCR3_FIDH");
  dcr3FidH.addValue(80.0);
  dcr3FidH.addError(0.0);

  dcr3FidA.setName("EL_DCR3_FIDA");
  dcr3FidA.addValue(49.0);
  dcr3FidA.addError(0.0);

  ccFidA.setName("EL_CC_FIDA");
  ccFidA.addValue(46.0);
  ccFidA.addError(0.0);
 
  ccFidB.setName("EL_CC_FIDB");
  ccFidB.addValue(35.0);
  ccFidB.addError(0.0);

  ccFidC.setName("EL_CC_FIDC");
  ccFidC.addValue(360.0);
  ccFidC.addError(0.0);

  zVertexLimitMin.setName("EL_VZ_MIN"); 
  zVertexLimitMin.addValue(-27.7302); 
  zVertexLimitMin.addError(0.00); 

  zVertexLimitMax.setName("EL_VZ_MAX"); 
  zVertexLimitMax.addValue(-22.6864); 
  zVertexLimitMax.addError(0.00); 

  ecUCoordMin.setName("EL_ECU_MIN"); 
  ecUCoordMin.addValue(70.0);
  ecUCoordMin.addError(0.0);
  ecUCoordMax.setName("EL_ECU_MAX"); 
  ecUCoordMax.addValue(400.0);
  ecUCoordMax.addError(0.0);

  ecVCoordMin.setName("EL_ECV_MIN"); 
  ecVCoordMin.addValue(0.0);
  ecVCoordMin.addError(0.0);
  ecVCoordMax.setName("EL_ECV_MAX"); 
  ecVCoordMax.addValue(362.0);
  ecVCoordMax.addError(0.0);

  ecWCoordMin.setName("EL_ECW_MIN"); 
  ecWCoordMin.addValue(0.0);
  ecWCoordMin.addError(0.0);
  ecWCoordMax.setName("EL_ECW_MAX"); 
  ecWCoordMax.addValue(395.0);
  ecWCoordMax.addError(0.0);

  ecEdepMin.setName("EL_EC_EDEP_MIN");
  ecEdepMin.addValue(0.05); 
  ecEdepMin.addError(0.00); 

  eidParameters.addParameterSet(ecNSigma);
  eidParameters.addParameterSet(NPheMin);
  eidParameters.addParameterSet(NPheMax);
  eidParameters.addParameterSet(ccFidA);
  eidParameters.addParameterSet(ccFidB);
  eidParameters.addParameterSet(ccFidC);
  eidParameters.addParameterSet(dcr1FidA);
  eidParameters.addParameterSet(dcr1FidH);
  eidParameters.addParameterSet(dcr3FidA);
  eidParameters.addParameterSet(dcr3FidH);
  eidParameters.addParameterSet(zVertexLimitMin);
  eidParameters.addParameterSet(zVertexLimitMax);
  eidParameters.addParameterSet(ecUCoordMin);
  eidParameters.addParameterSet(ecUCoordMax);
  eidParameters.addParameterSet(ecVCoordMin);
  eidParameters.addParameterSet(ecVCoordMax);
  eidParameters.addParameterSet(ecWCoordMin);
  eidParameters.addParameterSet(ecWCoordMax);
  eidParameters.addParameterSet(ecEdepMin); 
}

void findElectronIDParameters::SaveParameters(std::string outputFilename){
  eidParameters.saveParameters(outputFilename); 
}
 
void findElectronIDParameters::SaveHistograms(std::string outputFilename){

  TFile *outputFile = new TFile(outputFilename.c_str(), "recreate"); 

  for(int sect=0; sect<6; sect++){
    ecSamplingFraction[sect]->Write();
    ccThetaCorrelation[sect]->Write();

    ecMuGraph[sect]   ->Write();
    ecSigmaGraph[sect]->Write();
    ccMuGraph[sect]   ->Write();
    ccSigmaGraph[sect]->Write();
    ecPol3Mu[sect]    ->Write();
    ecPol3Sigma[sect] ->Write();
    

    for (int slice=0; slice<numberSlices; slice++){
      ecSlices[sect][slice]->Write();
    }

    for (int segment=0; segment<18; segment++){
      ccSlices[sect][segment]->Write();
    }
  }

  outputFile->Write();
  outputFile->Close(); 

}

#endif
