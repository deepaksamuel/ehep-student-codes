#include "Serc19RunActionMessenger.hh"
#include "Serc19RunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include <sstream>

Serc19RunActionMessenger::Serc19RunActionMessenger(Serc19RunAction* aRunAction)
:theRunAction(aRunAction)
{
  runDirectory = new G4UIdirectory("/Serc19/");
  runDirectory->SetGuidance("UI commands of this example");

  runDir = new G4UIdirectory("/Serc19/run/");
  runDir->SetGuidance("RunAction control");


  //runDirectory = new G4UIdirectory("/mysetrun/");
  //runDirectory->SetGuidance("My set commands.");

  runIDCmd = new G4UIcommand("/Serc19/run/SetRunID",this);
  runIDCmd->SetGuidance("Set run ID");
  runIDCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  G4UIparameter* p1 = new G4UIparameter("runID",'i',true);
  p1->SetDefaultValue(1000);
  runIDCmd->SetParameter(p1);

  OutputFileCmd = new G4UIcmdWithAString("/Serc19/run/output_file",this);
  OutputFileCmd->SetGuidance(" OutputFile name");
  OutputFileCmd->SetParameterName("output_file",true, true);
  OutputFileCmd->SetDefaultValue("output_file");

  OutputFileCmd2 = new G4UIcmdWithAString("/Serc19/run/output_file2",this);
  OutputFileCmd2->SetGuidance(" OutputFile2 name");
  OutputFileCmd2->SetParameterName("output_file2",true, true);
  OutputFileCmd2->SetDefaultValue("output_file2");

}

Serc19RunActionMessenger::~Serc19RunActionMessenger() {
  delete runIDCmd;
  delete runDirectory;
}

void Serc19RunActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue) {
  const char* nv = (const char*)newValue;
  if( command==runIDCmd ) {
    G4int id;
    std::istringstream is(nv);
    is >> id;
    
    theRunAction->SetRunID(id);
  }
  
  if( command == OutputFileCmd )
    {theRunAction->SetOutputFile(newValue);}
  
  if( command == OutputFileCmd2 )
    { theRunAction->SetOutputFile2(newValue);}
  


}

