// $Id: Serc19PrimaryGeneratorAction.cc,v 1.7 2003/09/15 15:38:18 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Serc19PrimaryGeneratorAction.hh"

#include "Serc19DetectorConstruction.hh"
#include "Serc19PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4Box.hh"

#include "math.h"
#include "CLHEP/Random/RandGauss.h"
#include "G4UnitsTable.hh"
#include "G4RunManager.hh"

//using namespace std;
//#include "Serc19DetectorParameterDef.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Serc19PrimaryGeneratorAction *Serc19PrimaryGeneratorAction::AnPointer;
Serc19PrimaryGeneratorAction::Serc19PrimaryGeneratorAction(
							 Serc19SimAnalysis *panalysis)
  :pAnalysis(panalysis)
{ 
  AnPointer =this;
  // G4cout<<" Initialized Serc19PrimaryGeneratorAction Constructor"<<endl;
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  //Default settings :
  
  SetRunNumber(0);
  //  SetOutputFile("simulation.root");
  //  SetFirstEvt(1);

  SetRndmFlag("off");
  SetRndmPIDFlag("off");

  SetPartId(211);

  SetIncEnergy(15.0*GeV);
  SetIncEnergySmr(100*MeV);

  SetIncDirection(G4ThreeVector(1.0,0.0,0.0));
  SetIncThetaSmr(2*mrad);
  SetIncPhiSmr(2*mrad);

  SetIncPosition(G4ThreeVector(0.0*cm,0.0*cm, 0.0*cm));
  SetIncVxSmr(0.03*mm);
  SetIncVySmr(0.03*mm);
  SetIncVzSmr(4.0*cm);
  nMultiplicity = 1;
  initialise = 0;
  
  //create a messenger for this class
  gunMessenger = new Serc19PrimaryGeneratorMessenger(this);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Serc19PrimaryGeneratorAction::~Serc19PrimaryGeneratorAction() {
  cout<<"Closing Serc19PrimaryGeneratorAction"<<endl;  
  if (particleGun)	{delete particleGun;}
  if (gunMessenger) {delete gunMessenger;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Serc19PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

  G4double   vx=0, vy=0, vz=0;
  
  if (initialise==0) {
    gunMessenger = new Serc19PrimaryGeneratorMessenger(this);
    g_nevt=-1;
    initialise = 1;
  }

  pAnalysis->irun = RunNumber;  //Keep an option that in a file, one may have more than two run number
  g_nevt++;

  pAnalysis->ievt=g_nevt;
  pAnalysis->ngent = (nMultiplicity <=int(pAnalysis->ngenmx)) ? nMultiplicity : pAnalysis->ngenmx;
  cout<<"multi "<< pAnalysis->ngent<<" "<<nMultiplicity<<endl;

  for (unsigned int ij=0; ij<pAnalysis->ngent; ij++)	{
    //Option to have multiple particle
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    int sign=1;
    if (rndmPIDFlag=="on") { 
      int idxx = 7*G4UniformRand();
      switch(idxx) { 
      case 0 :  partId=11; break;
      case 1 :  partId=-11; break;
      case 2 :  partId=22; break;
      case 3 :  partId=13; break;
      case 4 :  partId=-13; break;
      case 5 :  partId=211; break;
      case 6 :  partId=-211; break;
      default : partId=22; break;
      }
    } else if ( rndmPIDFlag=="charge") {
      if (G4UniformRand()>0.5) sign=-1;
    }

    //    partId=111;
    //    incEnergy = 5000;
    
    G4ParticleDefinition* particle = particleTable->FindParticle(partId*sign);
    particleGun->SetParticleDefinition(particle);
    
    G4ThreeVector ini_Dir(incDirection);
    double in_Energy = incEnergy*MeV;
    vx = incPosition.x()*mm;
    vy = incPosition.y()*mm;
    vz = incPosition.z()*mm;
    
    if (rndmFlag=="on") {
      G4double theta=ini_Dir.theta();
      
      if(incThetaSmr>=0){
	theta = G4RandGauss::shoot(0,incThetaSmr);
      } else if(incThetaSmr<0) {
	theta = incThetaSmr*(2*G4UniformRand()-1);
      }

      G4double phi=0;
      if(incPhiSmr>=0    && incDirection.theta() !=0) {
	phi = G4RandGauss::shoot(0,incPhiSmr);
      } else if(incPhiSmr<0 || incDirection.theta()==0) {
	phi = incPhiSmr*(2*G4UniformRand()-1);
      }
      
      ini_Dir.setTheta(ini_Dir.theta()+theta);
      ini_Dir.setPhi(ini_Dir.phi()+phi);
      
      if(incEnergySmr>=0) {
	in_Energy += G4RandGauss::shoot(0,incEnergySmr);
	if (in_Energy <1*MeV) in_Energy=1*MeV;
      } else if(incEnergySmr<0) {
	in_Energy += incEnergySmr*(2*G4UniformRand()-1);
	if (in_Energy <1*MeV) in_Energy=1*MeV;
      }
      
      if(incVxSmr>=0&&incVySmr>=0&&incVzSmr>=0) {
	vx += G4RandGauss::shoot(0,incVxSmr*mm);
	vy += G4RandGauss::shoot(0,incVySmr*mm);
	vz += G4RandGauss::shoot(0,incVzSmr*mm);
      } else {
	vx += incVxSmr*(2*G4UniformRand()-1.)*mm;
	vy += incVySmr*(2*G4UniformRand()-1.)*mm;
	vz += incVzSmr*(2*G4UniformRand()-1.)*mm;
      }
    }
    
    //    if (g_nevt >G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed()) {
    //      in_Energy = 0;
    //    } 

    particleGun->SetParticleMomentumDirection(ini_Dir);
    //	particleGun->SetParticleEnergy(in_Energy);
    particleGun->SetParticleMomentum(in_Energy);
    
    particleGun->SetParticlePosition(G4ThreeVector(vx, vy, vz));
    particleGun->GeneratePrimaryVertex(anEvent);
    
    cout<<"Gen Pid "<<particle->GetPDGEncoding()<<" "<<G4BestUnit(in_Energy, "Energy")<<" Dir "<<ini_Dir<<" Pos "<<G4BestUnit(vx, "Length")<<" "<<G4BestUnit(vy, "Length")<<" "<<G4BestUnit(vz, "Length")<<endl;
    
    if (ij<(int)pAnalysis->ngenmx) { 
      pAnalysis->pidin[ij] = particle->GetPDGEncoding();
      pAnalysis->posxin[ij] = vx;
      pAnalysis->posyin[ij] = vy;
      pAnalysis->poszin[ij] = vz;
      //can not use charge for neutral particle *(particle->GetPDGCharge());
      if (particle->GetPDGCharge()<0) { 
	pAnalysis->momin[ij] = -in_Energy/1000.0; //convert in GeV 
      } else {
	pAnalysis->momin[ij] = in_Energy/1000.0; ////convert in GeV
      }
      pAnalysis->thein[ij] = particleGun->GetParticleMomentumDirection().theta();
      pAnalysis->phiin[ij] = particleGun->GetParticleMomentumDirection().phi();
    }
    //  G4RunManager::GetRunManager()->RunTermination()
//     if (G4UniformRand()<0.0) {
//       cout<<" G4RunManager::GetRunManager()->TerminateOneEvent() crashed "<<  G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed()<<" "<<g_nevt<<endl;
//       G4RunManager::GetRunManager()->SetVerboseLevel(1);

//       G4RunManager::GetRunManager()->SetNumberOfEventsToBeProcessed(G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed()-1);

//       G4RunManager::GetRunManager()->TerminateEventLoop();
//       anEvent->SetEventID(anEvent->GetEventID()+1);
//       //      G4RunManager::GetRunManager()->RunTermination(); //Geomtry is not closed AbortEvent(); //RunTermination();

//       G4RunManager::GetRunManager()->BeamOn(G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed()-1-g_nevt);
//       G4RunManager::GetRunManager()->AbortEvent();
//       g_nevt++;
//     }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void Serc19PrimaryGeneratorAction::SetMultiplicity(G4int p) {
  nMultiplicity = p;
}


void Serc19PrimaryGeneratorAction::SetIncEnergy(G4double p) {
  incEnergy = p;
}

void Serc19PrimaryGeneratorAction::SetPartId(G4int p) {
  partId = p;
}


void Serc19PrimaryGeneratorAction::SetIncEnergySmr(G4double p) {
  incEnergySmr = p;
}

void Serc19PrimaryGeneratorAction::SetIncDirection(G4ThreeVector p){
  incDirection = p;
} 

void Serc19PrimaryGeneratorAction::SetIncThetaSmr(G4double p) {
  incThetaSmr =p;
}

void Serc19PrimaryGeneratorAction::SetIncPhiSmr(G4double p) {
  incPhiSmr =p;
}

void Serc19PrimaryGeneratorAction::SetIncPosition(G4ThreeVector p){
  incPosition = p;
} 

void Serc19PrimaryGeneratorAction::SetIncVxSmr(G4double p) {
  incVxSmr =p;
}

void Serc19PrimaryGeneratorAction::SetIncVySmr(G4double p) {
  incVySmr =p;
}

void Serc19PrimaryGeneratorAction::SetIncVzSmr(G4double p) {
  incVzSmr =p;
}

//void Serc19PrimaryGeneratorAction::SetOutputFile(G4String p) {
//  OutputFile = p; pAnalysis->text_outputFile=p;
//}
