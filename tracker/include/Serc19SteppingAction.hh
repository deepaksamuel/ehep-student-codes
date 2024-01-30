
// $Id: Serc19SteppingAction.hh,v 1.8 2003/09/15 15:38:14 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Serc19SteppingAction_h
#define Serc19SteppingAction_h 1

#include "G4UserSteppingAction.hh"

class Serc19EventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Serc19SteppingAction : public G4UserSteppingAction
{
  public:
    Serc19SteppingAction(Serc19EventAction*);
   ~Serc19SteppingAction();

    void UserSteppingAction(const G4Step*);
    
  private:
    Serc19EventAction*          eventaction;  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
