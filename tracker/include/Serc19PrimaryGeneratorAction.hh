//
// $Id: Serc19PrimaryGeneratorAction.hh,v 1.6 2003/09/15 15:38:14 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Serc19PrimaryGeneratorAction_h
#define Serc19PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "Serc19SimAnalysis.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
//#include <iostream.h>
//#include <fstream.h>
#include <iostream>
#include <fstream>
#include "TFile.h"
class G4ParticleGun;
class G4Event;
class Serc19DetectorConstruction;
class Serc19PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Serc19PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
	Serc19PrimaryGeneratorAction(Serc19SimAnalysis *panalysis);    
   ~Serc19PrimaryGeneratorAction();

  public:
    static Serc19PrimaryGeneratorAction* AnPointer;
    void GeneratePrimaries(G4Event*);

    void SetRunNumber(G4int p) { RunNumber = p;}

    void SetRndmPIDFlag(G4String p) { rndmPIDFlag = p;}
    void SetRndmFlag(G4String p) { rndmFlag = p;}

    void SetPartId(G4int p);
    void  SetMultiplicity(G4int p);
    void SetIncEnergy(G4double p);
    void SetIncEnergySmr(G4double p);

    void SetIncDirection(G4ThreeVector p); 
    void SetIncThetaSmr(G4double p);
    void SetIncPhiSmr(G4double p);

    void SetIncPosition(G4ThreeVector p); 
    void SetIncVxSmr(G4double p);
    void SetIncVySmr(G4double p);
    void SetIncVzSmr(G4double p);

    G4int                         InputFlag;        //Flag G4(0), Nuance (1), GENIE(2)
  private:
    G4ParticleGun*                particleGun;	  //pointer a to G4  class
   
    Serc19SimAnalysis *pAnalysis;                  //Pointer to rootfile & informations
    Serc19PrimaryGeneratorMessenger* gunMessenger; //messenger of this class
  
  G4int                         g_nevt;
  G4int                         initialise;
  G4int                         RunNumber;        //Run number for these events
  G4String                      OutputFile;	  //Output root file
  
  G4String                      rndmFlag;	  //flag for a rndm impact point and 4-momentum
  G4String                      rndmPIDFlag;	  //flag for a rndm partile(e/gamma/mu/pion)
  G4int                         partId; //STDHEP id of incident particle
  G4int                         nMultiplicity; //Number of particles to be simulated in an event
  G4double                      incEnergy;    //incident energy
  G4double                      incEnergySmr;    //smearing of incident energy

  G4ThreeVector                 incDirection; //Incident particle direction
  G4double                      incThetaSmr; //smearing of incident angle (in mrad)               
  G4double                      incPhiSmr; //smearing of incident angle (in mrad)    

  G4ThreeVector                 incPosition; //Position of incident particle
  G4double                      incVxSmr; //smearing of incident x-position 
  G4double                      incVySmr; //smearing of incident y-position 
  G4double                      incVzSmr; //smearing of incident z-position 

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


