// $Id: Serc19PrimaryGeneratorMessenger.hh,v 1.6 2002/12/16 16:37:26 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Serc19PrimaryGeneratorMessenger_h
#define Serc19PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

#include "Serc19PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Serc19PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    Serc19PrimaryGeneratorMessenger(Serc19PrimaryGeneratorAction*);
   ~Serc19PrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:

    Serc19PrimaryGeneratorAction* Serc19Action;

    G4UIdirectory*               IcalDir;
    G4UIdirectory*               gunDir; 

    G4UIcmdWithAnInteger*        RunNumberCmd;

    G4UIcmdWithAnInteger*        FirstEvtCmd;
    G4UIcmdWithAString*          RndmCmd;
    G4UIcmdWithAString*          RndmPIDCmd;
    G4UIcmdWithAnInteger*        partIdCmd;
    G4UIcmdWithAnInteger*        partMultiCmd;
    G4UIcmdWithADoubleAndUnit*   incEnergyCmd;
    G4UIcmdWithADoubleAndUnit*   incEnergySmrCmd;

    G4UIcmdWith3Vector*          incDirectionCmd;
    G4UIcmdWithADoubleAndUnit*   incThetaSmrCmd;
    G4UIcmdWithADoubleAndUnit*   incPhiSmrCmd;

    G4UIcmdWith3VectorAndUnit*   incPositionCmd;
    G4UIcmdWithADoubleAndUnit*   incVxSmrCmd;
    G4UIcmdWithADoubleAndUnit*   incVySmrCmd;
    G4UIcmdWithADoubleAndUnit*   incVzSmrCmd;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

