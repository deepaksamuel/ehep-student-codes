#include "Serc19TrkSD.hh"
#include "Serc19TrkHit.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "TRandom.h"
#include "TMath.h"

Serc19TrkSD::Serc19TrkSD(G4String name)
  :G4VSensitiveDetector(name),
   numberInCell(20000), InCell(0) 
{
  G4String HCname;
  collectionName.insert(HCname="TrkCollect");
  pAnalysis = Serc19SimAnalysis::AnPointer;

}

Serc19TrkSD::~Serc19TrkSD() {

}

void Serc19TrkSD::Initialize(G4HCofThisEvent* HCE) {
  static int HCID = -1;

  TrkCollection = new Serc19TrkHitsCollection
    (SensitiveDetectorName,collectionName[0]); 
  if(HCID<0) { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection(HCID,TrkCollection);
}

G4bool Serc19TrkSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
  G4double edep = aStep->GetTotalEnergyDeposit();

  G4TouchableHistory* theTouchable = (G4TouchableHistory*)( aStep->GetPreStepPoint()->GetTouchable() );
  //  int level = theTouchable->GetHistoryDepth();
  G4ThreeVector parmom = aStep->GetTrack()->GetMomentum();

  //  cout<<"Trk : aStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetName() "<<aStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetName()<<" "<<edep<<" "<<parmom<<" "<<parmom.mag()<<level<<" "<<aStep->GetTrack()->GetDefinition()->GetPDGEncoding()<<endl;

  if(edep==0.) return false;
  
  G4ThreeVector glbpos = 0.5*(aStep->GetPreStepPoint()->GetPosition() + aStep->GetPostStepPoint()->GetPosition()); //aStep->GetPreStepPoint()->GetPosition();

  //  G4int nLayer = theTouchable->GetCopyNumber(0) ;
  
  G4double atime = aStep->GetPreStepPoint()->GetGlobalTime()/(ns);
  
  G4ThreeVector localpos = theTouchable->GetHistory()->GetTopTransform().TransformPoint(glbpos); 
  
  if (parmom.mag()<10) { aStep->GetTrack()->SetTrackStatus(fStopAndKill);} //Manually kill a track

  G4double zz = 750.0 + localpos.z();
  int nIThe = int(zz/0.2); //200micon strip length //13bit
  G4double phi = acos(-1.0)/4. + localpos.phi();
  int nIPhi = int(phi*20000); // 50 micon strip width : 15bit
  
  pAnalysis->pPosRT->Fill(glbpos.rho(), edep);
  pAnalysis->pPosET->Fill(glbpos.eta(), edep);
  pAnalysis->pPosPT->Fill(glbpos.phi(), edep);
  
  edep /=keV;  //Convert in 
  
  unsigned int detid = theTouchable->GetCopyNumber(0);
  //  cout<<"detid "<<detid<< glbpos<<" "<<edep<<" "<<level<<" "<<theTouchable->GetCopyNumber(0)<< endl;
  pAnalysis->h_etaphi_trkenergy[detid]->Fill(glbpos.eta(), glbpos.phi(), edep);

  detid<<=13;  //For layer number
  detid +=nIThe;
  detid<<=15; 
  detid +=nIPhi;
  
  int oldCellId = -1;
  for (int ij=0; ij<InCell; ij++) {
    if (detid ==CellDetID[ij]) {oldCellId = ij; break;}
  }
      
  if (oldCellId ==-1 && InCell <numberInCell -1 ) {
    Serc19TrkHit* newHit = new Serc19TrkHit();
    newHit->SetHitId(detid);
    int pdgid = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
    newHit->SetpdgId(pdgid);
    newHit->SetEdep(edep);
    newHit->SetTime(atime);
    newHit->SetPos(glbpos);
    newHit->SetLocalTheta(localpos.theta());
    newHit->SetLocalPhi(localpos.phi());		     
    
    newHit->SetMom( aStep->GetTrack()->GetMomentum());
    
    InCell = TrkCollection->insert( newHit );
    CellDetID[InCell-1] = detid;
  }
  if (oldCellId >=0) {
    (*TrkCollection)[oldCellId]->AddEdep(edep);
    if (atime <(*TrkCollection)[oldCellId]->GetTime()) {
      (*TrkCollection)[oldCellId]->SetTime(atime);
    }
  }
  return true;
}

void Serc19TrkSD::EndOfEvent(G4HCofThisEvent*) { 
  InCell = 0;

  pAnalysis->nsimhtTk = 0;
  for (int ij=0; ij< TrkCollection->entries(); ij++) {
    if ((*TrkCollection)[ij]->GetEdep()>1.0 && pAnalysis->nsimhtTk < (int)pAnalysis->nsimhtmxTk) { //More than 1 KeV
      pAnalysis->detidTk[pAnalysis->nsimhtTk] =  (*TrkCollection)[ij]->GetHitId();
      pAnalysis->simpdgidTk[pAnalysis->nsimhtTk] =  (*TrkCollection)[ij]->GetpdgId();
      pAnalysis->simtimeTk[pAnalysis->nsimhtTk] = (*TrkCollection)[ij]->GetTime();
      pAnalysis->simenrTk[pAnalysis->nsimhtTk] = (*TrkCollection)[ij]->GetEdep();
      
      G4ThreeVector posvec1 = (*TrkCollection)[ij]->GetPos();
      pAnalysis->simvxTk[pAnalysis->nsimhtTk] = posvec1.x();
      pAnalysis->simvyTk[pAnalysis->nsimhtTk] = posvec1.y();
      pAnalysis->simvzTk[pAnalysis->nsimhtTk] = posvec1.z();
      
      G4ThreeVector momvec = (*TrkCollection)[ij]->GetMom();
      pAnalysis->simpxTk[pAnalysis->nsimhtTk] = momvec.x();
      pAnalysis->simpyTk[pAnalysis->nsimhtTk] = momvec.y();
      pAnalysis->simpzTk[pAnalysis->nsimhtTk] = momvec.z();
      
      pAnalysis->simloctheTk[pAnalysis->nsimhtTk] = (*TrkCollection)[ij]->GetLocalTheta();
      pAnalysis->simlocphiTk[pAnalysis->nsimhtTk] = (*TrkCollection)[ij]->GetLocalPhi();
      pAnalysis->nsimhtTk++;
      
    }
  }
}

void Serc19TrkSD::clear()
{
} 

void Serc19TrkSD::DrawAll()
{
} 

void Serc19TrkSD::PrintAll()
{
}   

