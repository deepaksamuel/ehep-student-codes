
// $Id: Serc19EventAction.hh,v 1.8 2003/09/15 15:38:14 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Serc19EventAction_h
#define Serc19EventAction_h 1

#include "Serc19SimAnalysis.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "TCanvas.h"

class Serc19EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Serc19EventAction : public G4UserEventAction
{
public:
  Serc19EventAction();
  ~Serc19EventAction();

public:
  void  BeginOfEventAction(const G4Event*);
  void    EndOfEventAction(const G4Event*);
    
  void Addcal0(G4double de, G4double dl) {Energycal0 += de; TrackLcal0 += dl;};
  void AddGap(G4double de, G4double dl) {EnergyGap += de; TrackLGap += dl;};
      
  //   void Addcal1(G4double de, G4double dl) {Energycal1 += de; TrackLcal1 += dl;};
  //   void Addcal2(G4double de, G4double dl) {Energycal2 += de; TrackLcal2 += dl;};
               
  void SetDrawFlag   (G4String val)  {drawFlag = val;};
  void SetPrintModulo(G4int    val)  {printModulo = val;};

  //  G4double GetEnergycal0 () { return Energycal0;}
  //  G4double GetEnergyGap () { return EnergyGap;}
  //  G4double GetTrackLcal0 () { return TrackLcal0;}
  //  G4double GetTrackLGap () { return TrackLGap;}

private:

  G4double  Energycal0, EnergyGap;
  G4double  TrackLcal0, TrackLGap;
     
  G4String  drawFlag;
  G4int     printModulo;
                             
  Serc19EventActionMessenger*  eventMessenger;

  G4int TrkCollID;
  G4int EcalCollID;
  G4int HclCollID;

  G4int nevent;
  //  const TCanvas* c0; 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
