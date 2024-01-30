
#ifndef Serc19HclSD_h
#define Serc19HclSD_h 1

#include "G4VSensitiveDetector.hh"
#include "Serc19HclHit.hh"
#include "Serc19SimAnalysis.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class Serc19HclSD : public G4VSensitiveDetector
{

  public:
      Serc19HclSD(G4String name);
      ~Serc19HclSD();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      Serc19HclHitsCollection *HclCollection;
  Serc19SimAnalysis *pAnalysis;
      int numberInCell, InCell;
  unsigned int CellDetID[20000];

  public:
};




#endif

