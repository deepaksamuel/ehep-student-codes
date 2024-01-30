
#include "Serc19HclSD.hh"
#include "Serc19HclHit.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"


Serc19HclSD::Serc19HclSD(G4String name)
  :G4VSensitiveDetector(name),numberInCell(20000), InCell(0) 
{
  G4String HCname;
  collectionName.insert(HCname="HclCollect");
  pAnalysis = Serc19SimAnalysis::AnPointer;

}

Serc19HclSD::~Serc19HclSD(){;}

void Serc19HclSD::Initialize(G4HCofThisEvent* HCE)
{
  static int HCID = -1;
  HclCollection = new Serc19HclHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  if(HCID<0)
  { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection(HCID,HclCollection);
}

G4bool Serc19HclSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
  
  G4double edep = aStep->GetTotalEnergyDeposit();
  
  //  if(edep==0.) return false;
  
  G4TouchableHistory* theTouchable = 
    (G4TouchableHistory*)( aStep->GetPreStepPoint()->GetTouchable() );
  
  //  int level = theTouchable->GetHistoryDepth();
  G4ThreeVector glbpos = 0.5*(aStep->GetPreStepPoint()->GetPosition() + aStep->GetPostStepPoint()->GetPosition());

  //  cout <<"hcal0HclSD "<<aStep->GetTrack()->GetVolume()->GetName()<<" x "<<aStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetName()<<" y "<<aStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName()<<" "<<level<<" "<<edep<<" "<<glbpos<<G4endl;

  pAnalysis->pPosRH->Fill(glbpos.rho(), edep);
  pAnalysis->pPosEH->Fill(glbpos.eta(), edep);
  pAnalysis->pPosPH->Fill(glbpos.phi(), edep);  
  
  int iphi = theTouchable->GetCopyNumber( 3 );
  int idepth = theTouchable->GetCopyNumber( 2 );
  int ieta = theTouchable->GetCopyNumber( 0 );

  pAnalysis->h2d_hcalenergy[idepth]->Fill(glbpos.eta(), glbpos.phi(), edep/MeV);

  //  cout<<"cal2 : aStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetName() "<<aStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetName()<<" "<<edep<<" "<<level<<" "<<glbpos.phi()/degree<<" "<<glbpos.theta()/degree<<" "<<iphi<<" "<<ieta<<" "<<idepth<<" "<< theTouchable->GetCopyNumber(5)<<" "<< theTouchable->GetCopyNumber(4)<<" "<< theTouchable->GetCopyNumber(1)<<endl;

  if (edep <10*eV) return false;
  
  edep /=keV;

  G4int nInT = G4int(10*aStep->GetPreStepPoint()->GetGlobalTime()/ns); //in 100 ps
  unsigned int detid = ieta;

  detid<<=6;
  detid +=iphi;
  detid<<=5;
  detid +=idepth;
  //  depth<<=13;

  int oldCellId = -1;
  for (int ij=0; ij<InCell; ij++) {
    if (detid ==CellDetID[ij]) {oldCellId = ij; break;}
  }
  
  if (oldCellId ==-1 && InCell <numberInCell -1 ) {
    Serc19HclHit* newHit = new Serc19HclHit();
    newHit->SetHitId(detid);
    newHit->SetEdep( edep );
    newHit->SetPos(glbpos);
    newHit->SetTime(nInT);

    InCell = HclCollection->insert( newHit );
    CellDetID[InCell-1] = detid;

  } else {
    (*HclCollection)[oldCellId]->AddEdep(edep);
    if (nInT <(*HclCollection)[oldCellId]->GetTime()) {
      (*HclCollection)[oldCellId]->SetTime(nInT);
    }
  }

  return true;
}

void Serc19HclSD::EndOfEvent(G4HCofThisEvent*) {

  InCell = 0;
  pAnalysis->nsimhtHL = 0; // HclCollection->entries();
  
  for (int ij=0; ij< HclCollection->entries(); ij++) {
    if ((*HclCollection)[ij]->GetEdep() >1.0 && pAnalysis->nsimhtHL< pAnalysis->nsimhtmxHL) { //More than 1 KeV
      pAnalysis->detidHL[pAnalysis->nsimhtHL] =  ((*HclCollection)[ij]->GetHitId()<<24)+min(16777215,int((*HclCollection)[ij]->GetEdep())); //2^24-1
      pAnalysis->timeHL[pAnalysis->nsimhtHL] = (*HclCollection)[ij]->GetEdep();
      pAnalysis->nsimhtHL++;
    }
  }
}

void Serc19HclSD::clear()
{
} 

void Serc19HclSD::DrawAll()
{;
} 

void Serc19HclSD::PrintAll()
{;
} 
