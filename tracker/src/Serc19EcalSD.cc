#include "Serc19EcalSD.hh"
#include "Serc19EcalHit.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

Serc19EcalSD::Serc19EcalSD(G4String name)
  :G4VSensitiveDetector(name), numberInCell(20000), InCell(0) 
{
  G4String HCname;
  collectionName.insert(HCname="EcalCollect");
  pAnalysis = Serc19SimAnalysis::AnPointer;
}

Serc19EcalSD::~Serc19EcalSD(){;}

void Serc19EcalSD::Initialize(G4HCofThisEvent* HCE)
{
  static int HCID = -1;
  EcalCollection = new Serc19EcalHitsCollection
    (SensitiveDetectorName,collectionName[0]); 
  if(HCID<0)
    { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection(HCID,EcalCollection);
}

G4bool Serc19EcalSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
  
  G4double edep = aStep->GetTotalEnergyDeposit(); // /MeV;
 
  G4TouchableHistory* theTouchable = 
    (G4TouchableHistory*)( aStep->GetPreStepPoint()->GetTouchable() );
  
  //  int level = theTouchable->GetHistoryDepth();
  G4ThreeVector glbpos = 0.5*(aStep->GetPreStepPoint()->GetPosition() + aStep->GetPostStepPoint()->GetPosition());
  G4ThreeVector localpos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(glbpos);
  
  pAnalysis->pPosRE->Fill(glbpos.rho(), edep);
  pAnalysis->pPosEE->Fill(glbpos.eta(), edep);
  pAnalysis->pPosPE->Fill(glbpos.phi(), edep);

  //  pAnalysis->h_ecalenergy->Fill(edep/GeV);
  pAnalysis->h_2decalenergy->Fill(glbpos.eta(), glbpos.phi(), edep/GeV);
  
  //  edep /=GeV; //Store in MeV binning

  int InThe = int(localpos.theta()/degree-45);  
  int InPhi = int(localpos.phi()/degree + 45);

  //  G4cout <<"hcal0EcalSD "<<aStep->GetTrack()->GetVolume()->GetName()<<" x "<<aStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetName()<<" y "<<aStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName()<<" "<<edep<<localpos<<" "<<localpos.theta()/degree<<" "<<InThe<<" "<<localpos.phi()/degree + 45<<" "<<InPhi<<endl;


  if (InThe <0 || InThe>127) return false;
  if (InPhi <0 || InPhi>127) return false;
  
  unsigned int detid = 128*InThe + InPhi;
  int oldCellId = -1;
  for (int ij=0; ij<InCell; ij++) {
    if (detid ==CellDetID[ij]) {oldCellId = ij; break;}
  }

  G4double atime = aStep->GetPreStepPoint()->GetGlobalTime()/(ns);
  if (oldCellId ==-1 && InCell <numberInCell -1 ) {
    Serc19EcalHit* newHit = new Serc19EcalHit();
    
    newHit->SetHitId(detid);
    newHit->SetEdep( edep );
    newHit->SetPos(localpos);
    newHit->SetTime(atime);
    
    InCell = EcalCollection->insert( newHit );
    CellDetID[InCell-1] = detid;
    
  } else {
    (*EcalCollection)[oldCellId]->AddEdep(edep);
    if (atime <(*EcalCollection)[oldCellId]->GetTime()) {
      (*EcalCollection)[oldCellId]->SetTime(atime);
    }
  }
  return true;
}

void Serc19EcalSD::EndOfEvent(G4HCofThisEvent*) {
  InCell = 0;
  pAnalysis->nsimhtEC = 0; //EcalCollection->entries();

  double toten = 0;
  for (int ij=0; ij< EcalCollection->entries(); ij++) {
    if ((*EcalCollection)[ij]->GetEdep() >1.e-3 && pAnalysis->nsimhtEC < pAnalysis->nsimhtmxEC) { //More than 1 MeV
      pAnalysis->detidEC[pAnalysis->nsimhtEC] =  ((*EcalCollection)[ij]->GetHitId()<<18)+min(262143,int((*EcalCollection)[ij]->GetEdep())); //2^18-1
      pAnalysis->energyEC[pAnalysis->nsimhtEC] = (*EcalCollection)[ij]->GetEdep();
      pAnalysis->thetaEC[pAnalysis->nsimhtEC] = (*EcalCollection)[ij]->GetPos().theta();
      pAnalysis->phiEC[pAnalysis->nsimhtEC] = (*EcalCollection)[ij]->GetPos().phi();
      pAnalysis->nsimhtEC++;
      toten +=(*EcalCollection)[ij]->GetEdep();
    }
  }
  cout<<"Total Energy in ECAL "<<G4BestUnit(toten, "Energy")<<" in "<<EcalCollection->entries()<<" towers"<<endl;
  
}


void Serc19EcalSD::clear()
{
} 

void Serc19EcalSD::DrawAll()
{
} 

void Serc19EcalSD::PrintAll()
{
} 
