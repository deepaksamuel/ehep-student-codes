////////////////////////////////////////////////////////////////////////////////
//
//  Event class for event information and store in root file as 
//  tuple and/or histograms
//
////////////////////////////////////////////////////////////////////////////////
#ifndef MULTISIM_H
#define MULTISIM_H
#include <vector>
using std::vector;
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TTree.h"
#include "TFile.h"
#include "globals.hh"
#include "TProfile.h"
//#include "G4SIunits.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <fstream>
using namespace std;

class Serc19SimAnalysis {
public:
  Serc19SimAnalysis();
  void OpenRootfiles(G4String out);
  void CloseRootfiles();
  ~Serc19SimAnalysis();
  
public:
  static Serc19SimAnalysis* AnPointer;
  G4String    text_outputFile;
  
  TFile *pVisFile;
  TFile *pRootFile;
  TFile *inputRootFile;  
  
  TTree *pEventTree;
  TTree *inputEventTree;
  
  TTree *rtree; 
  
  TH1F  *pPosRT;
  TH1F  *pPosET;  
  TH1F  *pPosPT;
  
  TH1F  *pPosRE;
  TH1F  *pPosEE;  
  TH1F  *pPosPE;
  
  TH1F  *pPosRH;
  TH1F  *pPosEH;  
  TH1F  *pPosPH;

  static const unsigned int nsilayer=13;
  static const unsigned int nhclwedge=36;
  static const unsigned int nhcalLayer=17;
  static const unsigned int nhcalEtaDiv=34;

  TH2F* h_etaphi_trkenergy[nsilayer];
  TH1F* h_trkenergy[nsilayer];
  TH1F* h_trkabsenergy[nsilayer];

  TH1F* h_parfenergy;
  TH1F* h_ecalenergy;
  TH2F* h_2decalenergy;
  TH1F* h_g10energy;

  TH1F* h_hcalenergy[nhcalLayer];
  TH1F* h_hcalabsenergy[nhcalLayer];
  
  TH2F* h2d_hcalenergy[nhcalLayer];
  TH2F* h2d_hcalabsenergy[nhcalLayer];

  int  ievent; //Event counter
  int  FirstEvt; //First to read from root file
  
  unsigned   irun;                // Run number of these events
  unsigned   ievt;                //Event number
  unsigned   ngent;
  
  G4float	ievt_wt;		//*GMa
  
  static const unsigned int ngenmx=50;
  G4int   pidin[ngenmx]; 	  //PID of incident particle
  G4float momin[ngenmx]; 	  //Energy of incident particle
  G4float thein[ngenmx];	  //Initial polar angle of incident particle
  G4float phiin[ngenmx];     //Initial azimuthal angle of incident particle 
  G4float posxin[ngenmx];	  //Initial X-position
  G4float posyin[ngenmx];     //Initial Y-position
  G4float poszin[ngenmx];     //Initial Z-position
  
  // For Simulation output of tracker
  unsigned int nsimhtTk;
  static const unsigned int nsimhtmxTk=2000;
  unsigned int detidTk[nsimhtmxTk];

  float simtimeTk[nsimhtmxTk];
  float simenrTk[nsimhtmxTk];

  int   simpdgidTk[nsimhtmxTk];
  float simvxTk[nsimhtmxTk]; 
  float simvyTk[nsimhtmxTk];
  float simvzTk[nsimhtmxTk];
  float simpxTk[nsimhtmxTk];
  float simpyTk[nsimhtmxTk];
  float simpzTk[nsimhtmxTk];

  float simloctheTk[nsimhtmxTk];
  float simlocphiTk[nsimhtmxTk];

  // For Simulation output of ECAL
  unsigned int nsimhtEC;
  static const unsigned int nsimhtmxEC=2000;
  unsigned int detidEC[nsimhtmxEC];
  float energyEC[nsimhtmxEC];
  float thetaEC[nsimhtmxEC];
  float phiEC[nsimhtmxEC]; 

  // For Simulation output of HCAL
  unsigned int nsimhtHL;
  static const unsigned int nsimhtmxHL=2000;
  unsigned long int detidHL[nsimhtmxHL];
  unsigned int timeHL[nsimhtmxHL];

  int fNtColId[10];
private:
  
};

#endif
