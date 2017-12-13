#ifndef Photon_h
#define Photon_h


 // Put your includes here 
#include "DBins.h" 
#include "DCut.h" 
#include "DEvent.h" 
#include "DSelection.h" 
#include "h22Event.h" 
#include "h22Reader.h" 

#include "histos.h"
#include "histos.cxx"

class Photon : public h22Reader {
    public:
        Photon();
        ~Photon();

    Histograms histos;

    void Loop();
    void ProcessEvent();
    void WriteParameters();
};
#endif
