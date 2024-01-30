//
// $Id: Serc19EventAction.cc,v 1.24 2005/05/30 14:24:31 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 
//GMAA Store begin and end hit position and muo energy at those points for the resolution of reconstruction algorithms, which is much better than true muon energy resolution

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// #include "g4root.hh"
//#include "Serc19Analysis.hh"
#include "Serc19EventAction.hh"

#include "Serc19TrkHit.hh"
#include "Serc19EcalHit.hh"
#include "Serc19HclHit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"

#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

#include "Randomize.hh"
#include <iomanip>
#include <utility>
#include "G4SDManager.hh"
#include "G4Trajectory.hh"
//#include <iomanip.h>

#include "TStyle.h"
#include "TMatrixD.h"
#include "TMath.h"
//#include "TVectorD.h"
//#include "TMatrixTBase.h"
#include "TMatrixDEigen.h"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Serc19EventAction::Serc19EventAction()
  :drawFlag("all"),printModulo(1000)
{
  TrkCollID = -1;
  EcalCollID = -1;
  HclCollID = -1;
  nevent = 0;

  cout<<"Print "<<drawFlag<<" "<< printModulo<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Serc19EventAction::~Serc19EventAction() {
  cout<<"Closing Serc19EventAction"<<endl;  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Serc19EventAction::BeginOfEventAction(const G4Event* evt) {  
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) {
    
    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
    //HepRandom::showEngineStatus();
  }
  
  //initialisation per event
  Energycal0 = EnergyGap = 0.;
  TrackLcal0 = TrackLGap = 0.;
  //  Energycal1 = Energycal2 = 0.;
  //  TrackLcal1 = TrackLcal2 = 0.;
  
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  if(TrkCollID<0||EcalCollID<0||HclCollID<0) {
  
    //  if(TrkCollID<0) {
    G4String colNam;
    TrkCollID = SDman->GetCollectionID(colNam="TrkCollect");
    EcalCollID = SDman->GetCollectionID(colNam="EcalCollect");
    HclCollID = SDman->GetCollectionID(colNam="HclCollect");

    cout<<"TrkCollID "<<TrkCollID<<" "<<EcalCollID<<" "<<HclCollID<<endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Serc19EventAction::EndOfEventAction(const G4Event* evt) {

  cout<<"1drawFlag0 "<<endl;
  nevent++;
  cout<<"1drawFlag0 "<<drawFlag<<endl;
  Serc19SimAnalysis *pAnalysis = Serc19SimAnalysis::AnPointer;
  if (abs(pAnalysis->momin[0])>1.e-3) {
  // get analysis manager
//   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
//   cout<<"2drawFlag0 "<<drawFlag<<endl;
//   analysisManager->FillH1(1, Energycal0);
//   analysisManager->FillH1(2, EnergyGap);
//   analysisManager->FillH1(3, TrackLcal0);
//   analysisManager->FillH1(4, TrackLGap);
  
//   // fill ntuple
//   //  analysisManager->FillNtupleDColumn(0, Energycal0);
//   //  analysisManager->FillNtupleDColumn(1, EnergyGap);
//   //  analysisManager->FillNtupleDColumn(2, TrackLcal0);
//   //  analysisManager->FillNtupleDColumn(3, TrackLGap);
//   //  analysisManager->AddNtupleRow();  
  
//   //  double Energycal0a  = Energycal0;
//   //  double Energycal0b = Energycal0;

//   //  double EnergyGapa = EnergyGap;
//   //  double TrackLcal0a = TrackLcal0;
//   //  double TrackLGapa = TrackLGap;

//   analysisManager->FillNtupleDColumn(1, pAnalysis->fNtColId[0], Energycal0);
//   analysisManager->FillNtupleDColumn(1, pAnalysis->fNtColId[1], EnergyGap);
//   analysisManager->FillNtupleDColumn(1, pAnalysis->fNtColId[2], TrackLcal0);
//   analysisManager->FillNtupleDColumn(1, pAnalysis->fNtColId[3], TrackLGap);
//   analysisManager->AddNtupleRow(1);  

//   analysisManager->FillNtupleDColumn(2, pAnalysis->fNtColId[4], Energycal0);
//   analysisManager->AddNtupleRow(2); 

//   analysisManager->FillNtupleDColumn(3, pAnalysis->fNtColId[5], TrackLcal0);
//   analysisManager->AddNtupleRow(3); 
//   analysisManager->AddNtupleRow(3); 

//   analysisManager->FillNtupleDColumn(4, pAnalysis->fNtColId[6], EnergyGap);
//   analysisManager->FillNtupleDColumn(4, pAnalysis->fNtColId[7], TrackLcal0);
//   analysisManager->FillNtupleDColumn(4, pAnalysis->fNtColId[8], TrackLGap);
//   analysisManager->AddNtupleRow(4); 
//   analysisManager->AddNtupleRow(4); 
//   analysisManager->AddNtupleRow(4);  

  G4int evtNb = evt->GetEventID();
  //  cout<<"# "<<evtNb<<endl;
  if (evtNb%printModulo == 0) {
    G4cout<< G4endl
	  << "    Absrober: total energy: " << std::setw(7)
	  << G4BestUnit(EnergyGap,"Energy")
	  << "        total track length: " << std::setw(7)
	  << G4BestUnit(TrackLGap,"Length")
	  << G4endl
	  << "RPC_gas: cal0 total energy: " << std::setw(7)
	  << G4BestUnit(Energycal0,"Energy")
	  << "        total track length: " << std::setw(7)
	  << G4BestUnit(TrackLcal0,"Length")
	  << G4endl;
  }
  //    pAnalysis->simtotabenr = EnergyGap;
  //    pAnalysis->simtotrpcenr= Energycal0;
  //    pAnalysis->simtotablen = TrackLGap;
  //    pAnalysis->simtotrpclen= TrackLcal0;
  //GMA Possible point of error
//   cout<<"drawFlag0 "<<drawFlag<<endl;
  G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();
  if (pVisManager) {
    G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
    
		cout<<"drawFlag1 "<<drawFlag<<" "<<n_trajectories<<endl;
     for (G4int ij=0; ij<n_trajectories; ij++) {
       G4VTrajectory* trj = ((*(evt->GetTrajectoryContainer()))[ij]);
       if (drawFlag == "all") { pVisManager->Draw(*trj); //
       } else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.)) {
 	pVisManager->Draw(*trj); // GMA14 ,100);
       } else if ((drawFlag == "neutral")&&(trj->GetCharge() == 0.)) {
 	pVisManager->Draw(*trj); //GMA14 ,100);
       }
     }
  } //  // if (pVisManager)

      
  if(TrkCollID>0 || EcalCollID>0 || HclCollID>0) {
      
    G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
    Serc19TrkHitsCollection* EHC0 = 0;
    Serc19EcalHitsCollection* EHC1 = 0;
    Serc19HclHitsCollection* EHC2 = 0;
    if(HCE) {
      EHC0 = (Serc19TrkHitsCollection*)(HCE->GetHC(TrkCollID));
      EHC1 = (Serc19EcalHitsCollection*)(HCE->GetHC(EcalCollID));
      EHC2 = (Serc19HclHitsCollection*)(HCE->GetHC(HclCollID));
    }
    
    if(EHC0 && TrkCollID>=0) {
      int n_hit = EHC0->entries();
      //    pAnalysis->ascii_output <<n_hit<<" ";
      
      G4double totET[pAnalysis->nsilayer] = {0};
      for (int ij=0; ij<n_hit; ij++) {
	unsigned il = ((*EHC0)[ij]->GetHitId()>>28)&0xF;
	if (il<pAnalysis->nsilayer) {
	  totET[il] +=(*EHC0)[ij]->GetEdep();
	} else {
	  cout<<"Wrong layer in Tracker "<< il<<endl;
	}
      }
      for (unsigned ij=0; ij<pAnalysis->nsilayer; ij++) {
	pAnalysis->h_trkenergy[ij]->Fill(log10(totET[ij]));
      }
    }

    if(EHC1 && EcalCollID>=0) {
      G4double totEE = 0; 
      for (int ij=0; ij<EHC1->entries(); ij++) { 
	totEE +=(*EHC1)[ij]->GetEdep();
      }
      pAnalysis->h_ecalenergy->Fill(log10(totEE));
    }
    //    cout<<"3Eventaction "<<TrkCollID<<" "<<EcalCollID<<" "<<HclCollID<<endl;
    if(EHC2 && HclCollID>=0) {
      G4double totEH[pAnalysis->nhcalLayer] = {0};
      for (int ij=0; ij<EHC2->entries(); ij++) {
	unsigned il = ((*EHC2)[ij]->GetHitId())&0x1F;
	if (il<pAnalysis->nhcalLayer) {
	  totEH[il] +=(*EHC2)[ij]->GetEdep();
	} else {
	  cout<<"Wrong layer in HCAL "<< il<<endl;
	}
      }

      for (unsigned ij=0; ij<pAnalysis->nhcalLayer; ij++) {
	pAnalysis->h_hcalenergy[ij]->Fill(log10(totEH[ij]));
      }
    }
  } // if(TrkCollID>0 || EcalCollID>0 || HclCollID>0)
  //  cout<<" (evt->GetEventID())%(printModulo/100) "<<evt->GetEventID()<<" "<<drawFlag<<" "<<printModulo<<" "<<(printModulo/100.)<<" "<<(evt->GetEventID())%(printModulo/100)<<endl;
  //  if ((evt->GetEventID())%(printModulo/100)==0) {
    cout<<"evt# "<<evt->GetEventID()<< " TrkHit "<<pAnalysis->nsimhtTk<<" EcalHit "<<pAnalysis->nsimhtEC<<" HcalHit "<<pAnalysis->nsimhtHL<<endl;
    //}

  pAnalysis->pEventTree->Fill();
  }
}
