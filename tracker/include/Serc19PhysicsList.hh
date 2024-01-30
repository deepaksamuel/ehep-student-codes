#ifndef Serc19PhysicsList_h
#define Serc19PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "CompileTimeConstraints.hh"
#include "Serc19DetectorConstruction.hh"
#include "G4ProductionCutsTable.hh"


class Serc19DetectorConstruction;
class G4VPhysicsConstructor;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Serc19PhysicsList: public G4VModularPhysicsList
{
public:
    Serc19PhysicsList();
   ~Serc19PhysicsList();
   
    void Serc19AddPhysicsList(const G4String& name);
    
   
   public:
  // SetCuts() 
  virtual void SetCuts();
  void SetRegionCut(G4double){;}
  void AddPhysicsList(const G4String& name);
  void GetRange(G4double);
  
   private:
 // enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };
  //   Serc19DetectorConstruction* pDet;
   G4ProductionCuts* cuts;
};


#endif
