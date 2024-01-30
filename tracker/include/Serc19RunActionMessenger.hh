#ifndef Serc19RunActionMessenger_h
#define Serc19RunActionMessenger_h 1

class G4UIdirectory;
class Serc19RunAction;
#include "G4UImessenger.hh"
#include "globals.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

class Serc19RunActionMessenger: public G4UImessenger
{
public:
    Serc19RunActionMessenger(Serc19RunAction*);
    ~Serc19RunActionMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);

  private:
    Serc19RunAction * theRunAction;
    
  private: //commands
    G4UIdirectory *             runDirectory;
    G4UIdirectory*              runDir;
    G4UIcommand *               runIDCmd ;

    G4UIcmdWithAString*          OutputFileCmd;
    G4UIcmdWithAString*          OutputFileCmd2;


};

#endif

