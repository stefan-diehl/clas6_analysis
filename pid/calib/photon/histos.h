#ifndef histo_h
#define histo_h

#include <iostream>

#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"

#include "h22Event.h"

class Histograms{
	 public:
	 Histograms(){}
	 ~Histograms(){}

	 // Define histograms/graphs here
	 TH2I *betaVsP[7];

	 void Init();
	 void Fill(h22Event event, int ipart);
	 void Save(std::string outputFilename);
};

#endif
