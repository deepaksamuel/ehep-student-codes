#ifndef Serc19TrkSD_h
#define Serc19TrkSD_h 1
#include "G4VSensitiveDetector.hh"
#include "Serc19TrkHit.hh"
#include "Serc19SimAnalysis.hh"
#include <vector>

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class Serc19TrkSD : public G4VSensitiveDetector
{
  
public:
  Serc19TrkSD(G4String name);
  ~Serc19TrkSD();
  
  void Initialize(G4HCofThisEvent*HCE);
  G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
  void EndOfEvent(G4HCofThisEvent*HCE);
  void clear();
  void DrawAll();
  void PrintAll();
  
private:

  Serc19TrkHitsCollection *TrkCollection;
  Serc19SimAnalysis *pAnalysis;

  double ShiftInR;
  double ShiftInT;
  double ShiftInP;

  int numberInCell, InCell;
  unsigned int CellDetID[20000];

  public:
  float histxmn, histxmx, histymn, histymx, histzmn, histzmx;

};




#endif

