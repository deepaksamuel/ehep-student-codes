
// $Id: Serc19PrimaryGeneratorMessenger.cc,v 1.8 2002/12/16 16:37:27 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Serc19PrimaryGeneratorMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Serc19PrimaryGeneratorMessenger::Serc19PrimaryGeneratorMessenger(
                                          Serc19PrimaryGeneratorAction* Serc19Gun)
:Serc19Action(Serc19Gun)
{
  
  IcalDir = new G4UIdirectory("/Serc19/");
  IcalDir->SetGuidance("UI commands of this example");

  gunDir = new G4UIdirectory("/Serc19/gun/");
  gunDir->SetGuidance("PrimaryGenerator control");

  RunNumberCmd = new G4UIcmdWithAnInteger("/Serc19/gun/runNumber",this);
  RunNumberCmd->SetGuidance("Run number for these events");
  RunNumberCmd->SetParameterName("Run_number",true, true);
  RunNumberCmd->SetDefaultValue(1);

  RndmCmd = new G4UIcmdWithAString("/Serc19/gun/rndm",this);
  RndmCmd->SetGuidance("Shoot randomly the incident particle.");
  RndmCmd->SetGuidance("  Choice : on(default), off");
  RndmCmd->SetParameterName("choice",true, true);
  RndmCmd->SetDefaultValue("on");
  RndmCmd->SetCandidates("on off");

  RndmPIDCmd = new G4UIcmdWithAString("/Serc19/gun/rndmPID",this);
  RndmPIDCmd->SetGuidance("Shoot randomly different typ of particle.");
  RndmPIDCmd->SetGuidance("  Choice : on(default), off");
  RndmPIDCmd->SetParameterName("choice",true, true);
  RndmPIDCmd->SetDefaultValue("on");
  RndmPIDCmd->SetCandidates("on off charge");

  partIdCmd = new G4UIcmdWithAnInteger("/Serc19/gun/pid",this);
  partIdCmd->SetGuidance("Set incident energy of particle ID (PDG)");
  partIdCmd->SetParameterName("partID",true, true);
  partIdCmd->SetDefaultValue(13);
  //  partIdCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  partMultiCmd = new G4UIcmdWithAnInteger("/Serc19/gun/mult",this);
  partMultiCmd->SetGuidance("Set number of particles in a event");
  partMultiCmd->SetParameterName("Multiplicity",true, true);
  partMultiCmd->SetDefaultValue(1);

  incEnergyCmd = new G4UIcmdWithADoubleAndUnit("/Serc19/gun/energy",this);
  incEnergyCmd->SetGuidance("Set incident energy of particle.");
  incEnergyCmd->SetParameterName("Energy",true, true);
  incEnergyCmd->SetDefaultUnit("GeV");
  incEnergyCmd->SetDefaultValue(15.);
  incEnergyCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");
  //  incEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  incEnergySmrCmd = new G4UIcmdWithADoubleAndUnit("/Serc19/gun/ensmear",this);
  incEnergySmrCmd->SetGuidance("Set incident energy of particle.");
  incEnergySmrCmd->SetParameterName("EnergySmr",true, true);
  incEnergySmrCmd->SetDefaultUnit("GeV");
  incEnergySmrCmd->SetDefaultValue(0.1);
  incEnergySmrCmd->SetUnitCandidates("eV keV MeV GeV TeV");
  //  incEnergySmrCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  incDirectionCmd = new G4UIcmdWith3Vector("/Serc19/gun/incdir",this);
  incDirectionCmd->SetGuidance("Set the incident direction ");
  incDirectionCmd->SetParameterName("cosx","cosy","cosz",true, true);
  incDirectionCmd->SetDefaultValue(G4ThreeVector(1., 0., 0.));
  //  incDirectionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  incThetaSmrCmd = new G4UIcmdWithADoubleAndUnit("/Serc19/gun/thsmear",this);
  incThetaSmrCmd->SetGuidance("Set incident energy of particle.");
  incThetaSmrCmd->SetParameterName("ThetaSmr",true, true);
  incThetaSmrCmd->SetDefaultUnit("mrad");
  incThetaSmrCmd->SetDefaultValue(10.0);
  incThetaSmrCmd->SetUnitCandidates("degree rad mrad");
  //  incThetaSmrCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  incPhiSmrCmd = new G4UIcmdWithADoubleAndUnit("/Serc19/gun/phsmear",this);
  incPhiSmrCmd->SetGuidance("Set incident energy of particle.");
  incPhiSmrCmd->SetParameterName("PhiSmr",true, true);
  incPhiSmrCmd->SetDefaultUnit("mrad");
  incPhiSmrCmd->SetDefaultValue(10.0);
  incPhiSmrCmd->SetUnitCandidates("degree rad mrad");
  //  incPhiSmrCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  incPositionCmd = new G4UIcmdWith3VectorAndUnit("/Serc19/gun/incpos",this);
  incPositionCmd->SetGuidance("Set the incident position ");
  incPositionCmd->SetParameterName("x","y","z",true, true);
  incPositionCmd->SetDefaultValue(G4ThreeVector(0.,0.,-200.));
  incPositionCmd->SetDefaultUnit("cm");
  incPositionCmd->SetUnitCandidates("microm mm cm m km");
  //  incPositionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  incVxSmrCmd = new G4UIcmdWithADoubleAndUnit("/Serc19/gun/vxsmear",this);
  //  incVxSmrCmd->SetGuidance("Set Gaussian smearing in X-position");
  incVxSmrCmd->SetGuidance("Set uniform smearing in X-position");
  incVxSmrCmd->SetParameterName("VxSmr",true, true);
  incVxSmrCmd->SetDefaultUnit("cm");
  incVxSmrCmd->SetDefaultValue(0.01);
  //  incVxSmrCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  incVySmrCmd = new G4UIcmdWithADoubleAndUnit("/Serc19/gun/vysmear",this);
  //  incVySmrCmd->SetGuidance("Set Gaussian smearing in Y-position");
  incVySmrCmd->SetGuidance("Set uniform smearing in Y-position");
  incVySmrCmd->SetParameterName("VySmr",true, true);
  incVySmrCmd->SetDefaultUnit("cm");
  incVySmrCmd->SetDefaultValue(0.01);
  //  incVySmrCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  incVzSmrCmd = new G4UIcmdWithADoubleAndUnit("/Serc19/gun/vzsmear",this);
  incVzSmrCmd->SetGuidance("Set uniform smearing in Z-position");
  incVzSmrCmd->SetParameterName("VzSmr",true, true);
  incVzSmrCmd->SetDefaultUnit("cm");
  incVzSmrCmd->SetDefaultValue(5);
  //  incVzSmrCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  //InputOutputCmd = new G4UIcmdWithAnInteger("/Serc19/gun/inout",this);
  //InputOutputCmd->SetGuidance("Input/output option, 0:GEN->RECO, 1:GEN->DIGI, 2:GEN->SIM, 3: SIM -> RECO, 4: SIM -> DIFI, 5 : DIGI -> RECO");
  //InputOutputCmd->SetParameterName("InputOutput",true, true);
  //InputOutputCmd->SetDefaultValue(0);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Serc19PrimaryGeneratorMessenger::~Serc19PrimaryGeneratorMessenger() {
  cout<<"Closing Serc19PrimaryGeneratorMessenger"<<endl;  

  if (RunNumberCmd){     delete RunNumberCmd;}

  if (RndmCmd){          delete RndmCmd;}
  if (partIdCmd){        delete  partIdCmd;}
  if (partMultiCmd){     delete partMultiCmd;}

  if (incEnergyCmd){     delete incEnergyCmd;}
  if (incEnergySmrCmd){  delete incEnergySmrCmd;}

  if (incDirectionCmd){  delete  incDirectionCmd;}
  if (incThetaSmrCmd){   delete  incThetaSmrCmd;}
  if (incPhiSmrCmd){     delete  incPhiSmrCmd;}

  if (incPositionCmd){   delete  incPositionCmd;}
  if (incVxSmrCmd){      delete incVxSmrCmd;}
  if (incVySmrCmd){      delete incVySmrCmd;}
  if (incVzSmrCmd){      delete incVzSmrCmd;}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Serc19PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) { 
  
  if( command == RunNumberCmd )
    { Serc19Action->SetRunNumber(RunNumberCmd->GetNewIntValue(newValue));}
  
  if( command == RndmCmd )
    { Serc19Action->SetRndmFlag(newValue);}
  
  if( command == RndmPIDCmd )
    { Serc19Action->SetRndmPIDFlag(newValue);}
  

  if( command == partIdCmd )
    { Serc19Action->SetPartId(partIdCmd->GetNewIntValue(newValue));}

  if( command == partMultiCmd )
    { Serc19Action->SetMultiplicity(partMultiCmd->GetNewIntValue(newValue));}
  
  if( command == incEnergyCmd )
    { Serc19Action->SetIncEnergy(incEnergyCmd->GetNewDoubleValue(newValue));}
  
  if( command == incEnergySmrCmd )
    { Serc19Action->SetIncEnergySmr(incEnergySmrCmd->GetNewDoubleValue(newValue));}
  
  if( command == incDirectionCmd )
    { Serc19Action->SetIncDirection(incDirectionCmd->GetNew3VectorValue(newValue));}
  
  if( command == incThetaSmrCmd )
    { Serc19Action->SetIncThetaSmr(incThetaSmrCmd->GetNewDoubleValue(newValue));}
  
  if( command == incPhiSmrCmd )
    { Serc19Action->SetIncPhiSmr(incPhiSmrCmd->GetNewDoubleValue(newValue));}
  
  if( command == incPositionCmd )
    { Serc19Action->SetIncPosition(incPositionCmd->GetNew3VectorValue(newValue));}
  
  if( command == incVxSmrCmd )
    { Serc19Action->SetIncVxSmr(incVxSmrCmd->GetNewDoubleValue(newValue));}
  
  if( command == incVySmrCmd )
    { Serc19Action->SetIncVySmr(incVySmrCmd->GetNewDoubleValue(newValue));}
  
  if( command == incVzSmrCmd )
    { Serc19Action->SetIncVzSmr(incVzSmrCmd->GetNewDoubleValue(newValue));}
  
  //if( command == InputOutputCmd )
  // { Serc19Action->SetInputOutput(InputOutputCmd->GetNewIntValue(newValue));}
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

