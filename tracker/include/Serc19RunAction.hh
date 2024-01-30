// $Id: Serc19RunAction.hh,v 1.8 2003/09/15 15:38:14 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Serc19RunAction_h
#define Serc19RunAction_h 1
#include "TFile.h"
#include "TH1.h"
#include "G4UserRunAction.hh"
#include "Serc19SimAnalysis.hh"
#include "globals.hh"
#include "iostream"
#include <fstream>
#include "G4VAnalysisManager.hh"
#include "Serc19RunActionMessenger.hh"

using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;

class Serc19RunAction : public G4UserRunAction
{
  public:
    Serc19RunAction();
   ~Serc19RunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
  void SetRunID(G4int run_id) { run_ID = run_id ;}

  void SetOutputFile(G4String p) { r_file_title = p;}
  void SetOutputFile2(G4String p) { r_file_title2 = p;}

    Serc19RunActionMessenger *theRunActMessenger;

    private:
    G4String r_file_title;
    G4String r_file_title2;

    G4int run_no;
    G4int evt_no;
    G4int run_ID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

