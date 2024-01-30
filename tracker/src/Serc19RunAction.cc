//
// $Id: Serc19RunAction.cc,v 1.15 2003/11/25 16:50:13 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//#include "Serc19Analysis.hh"

#include "Serc19RunAction.hh"
// #include "g4root.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Serc19RunAction::Serc19RunAction() {
  Serc19SimAnalysis *pAnalysis = Serc19SimAnalysis::AnPointer;
  theRunActMessenger = new Serc19RunActionMessenger(this);

  for (int ij=0; ij<10; ij++) {pAnalysis->fNtColId[ij]=-1; }

  SetOutputFile("test");
  SetOutputFile2("g4_test.root");

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Serc19RunAction::~Serc19RunAction() {
  
delete theRunActMessenger; theRunActMessenger=0;
 cout<<"Closing Serc19RunAction"<<endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Serc19RunAction::BeginOfRunAction(const G4Run* aRun) {

  Serc19SimAnalysis *pAnalysis = Serc19SimAnalysis::AnPointer;

  G4cout << "### RunID "<< aRun->GetRunID() <<" Nevt " <<aRun->GetNumberOfEvent()<< G4endl;
  //inform the runManager to save random number seed
  
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  //  gRandom->SetSeed(1327512);
  char output_title[100];
  char name[100];
  sprintf(name, r_file_title);
  int ij=aRun->GetRunID();

  sprintf(output_title, "%s_run%d.root", name, ij*10); //body of the file name
  
  pAnalysis->OpenRootfiles(output_title);
  cout <<"output_title = "<<output_title<<endl;
  //cout<<"#############################"<<endl;

//   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
//   cout<<"r_file_title2 r_file_title2 "<<r_file_title2<<" yyy "<<r_file_title<<endl;
//   analysisManager->OpenFile(r_file_title2); //r_file_title2);
//   // Create directories
//   analysisManager->SetHistoDirectoryName("serc19_histograms");
//   analysisManager->SetNtupleDirectoryName("serc19_ntuple");
//   analysisManager->SetVerboseLevel(1);
//   analysisManager->SetFirstHistoId(1);
//   analysisManager->SetFirstNtupleId(1);

//   // Book histograms, rootuple
//   //
 
//   // Creating histograms
//   analysisManager->CreateH1("Eabs","Edep in absorber", 100, 0., 800*MeV);
//   analysisManager->CreateH1("Egap","Edep in gap", 100, 0., 100*MeV);
//   analysisManager->CreateH1("Labs","trackL in absorber", 100, 0., 1*m);
//   analysisManager->CreateH1("Lgap","trackL in gap", 100, 0., 50*cm);

//   // Creating ntuple
//   //
//   //  analysisManager->CreateNtuple("B4", "Edep and TrackL");
//   //  analysisManager->CreateNtupleDColumn("Eabs");
//   //  analysisManager->CreateNtupleDColumn("Egap");
//   //  analysisManager->CreateNtupleDColumn("Labs");
//   //  analysisManager->CreateNtupleDColumn("Lgap");
//   //  analysisManager->FinishNtuple();

//   analysisManager->CreateNtuple("B4", "Edep and TrackL");
//   pAnalysis->fNtColId[0] = analysisManager->CreateNtupleDColumn("Eabs");
//   pAnalysis->fNtColId[1] = analysisManager->CreateNtupleDColumn("Egap");
//   pAnalysis->fNtColId[2] = analysisManager->CreateNtupleDColumn("Labs");
//   pAnalysis->fNtColId[3] = analysisManager->CreateNtupleDColumn("Lgap");
//   analysisManager->FinishNtuple();

//  //Create one ntuple
//   analysisManager -> CreateNtuple("b4a", "Primary");
//   pAnalysis->fNtColId[4] = analysisManager -> CreateNtupleDColumn("Eabs1");
//   analysisManager -> FinishNtuple();

//   //Create Second ntuple
//   analysisManager-> CreateNtuple("b4b", "Secondary");
//   pAnalysis->fNtColId[5] = analysisManager -> CreateNtupleDColumn("Labs1");
//   analysisManager -> FinishNtuple();

//   //creating third ntuple
//   analysisManager -> CreateNtuple("b4c", "Tertiary");
//   pAnalysis->fNtColId[6] = analysisManager -> CreateNtupleDColumn("Egap2");
//   pAnalysis->fNtColId[7] = analysisManager -> CreateNtupleDColumn("Labs2");
//   pAnalysis->fNtColId[8] = analysisManager -> CreateNtupleDColumn("Lgap2");
//   analysisManager -> FinishNtuple();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Serc19RunAction::EndOfRunAction(const G4Run* ) {

  Serc19SimAnalysis *pAnalysis = Serc19SimAnalysis::AnPointer;
  cout<<"Serc19SimAnalysis *pAnalysis = Serc19SimAnalysis::AnPointer;"<<endl;
  pAnalysis->CloseRootfiles();

//   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
//   analysisManager->Write();
//   analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
