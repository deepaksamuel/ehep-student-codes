
#include "Serc19DetectorConstruction.hh"
// #include "B2ActionInitialization.hh"
#include "Serc19PhysicsList.hh"
#include "Serc19PrimaryGeneratorAction.hh"
#include "Serc19RunAction.hh"
#include "Serc19EventAction.hh"
#include "Serc19SteppingAction.hh"
#include "Serc19SimAnalysis.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

// #ifdef G4VIS_USE
// #include "G4VisExecutive.hh"
// #endif

// #ifdef G4UI_USE
// #include "G4UIExecutive.hh"
// #endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {
  G4cout <<"argc "<<argc<<" "<<argv[0]<<" "<<argv[1]<<" "<<argv[2]<<G4endl;
  // Choose the Random engine
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) { ui = new G4UIExecutive(argc, argv); }

  long seeds[2]={12334457,1239075};
  time_t systime = time(NULL);
  seeds[0] = (long) systime;
  seeds[1] = (long) (systime*G4UniformRand());
  
  CLHEP::HepRandom::setTheSeeds(seeds);
  CLHEP::HepRandom::showEngineStatus();

  // G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
  // Construct the default run manager
  
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
#else
  G4RunManager * runManager = new G4RunManager;
  G4RunManager::GetRunManager()->SetRunIDCounter(1);
#endif

  // Set mandatory initialization classes
  Serc19SimAnalysis *panalysis = new Serc19SimAnalysis();
  Serc19DetectorConstruction* detector = new Serc19DetectorConstruction(panalysis);
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Serc19PhysicsList());
  // G4VModularPhysicsList* physicsList = new FTFP_BERT;
  // physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  // runManager->SetUserInitialization(physicsList);
    
  // Set user action classes
  runManager->SetUserAction(new Serc19PrimaryGeneratorAction(panalysis));
  runManager->SetUserAction(new Serc19RunAction);
  Serc19EventAction* eventaction = new Serc19EventAction;
  runManager->SetUserAction(eventaction);
  runManager->SetUserAction(new Serc19SteppingAction(eventaction));
  // runManager->SetUserInitialization(new B2ActionInitialization());
  
  // Initialize G4 kernel
  runManager->Initialize();
  
// #ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
// #endif

  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( ! ui ) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  // delete eventaction;
  // delete detector;
  delete panalysis;
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
