#ifndef Serc19DetectorConstruction_h
#define Serc19DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4SubtractionSolid.hh"
#include "G4Ellipsoid.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "Serc19TrkSD.hh"
#include "Serc19EcalSD.hh"
#include "Serc19HclSD.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4SDManager.hh"
#include "Serc19SimAnalysis.hh"

class G4Box;
class G4Trd;
class G4Tubs;
class G4Sphere;
class G4Ellipsoid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class G4MaterialCutsCouple;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Serc19DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
  Serc19DetectorConstruction(Serc19SimAnalysis *panalysis);
   ~Serc19DetectorConstruction();

  public:
  void DefineMaterials();
     void SetUniformMagField(G4double);
     
     G4VPhysicalVolume* Construct();

  //     void UpdateGeometry();
     
  public:
  
  
  const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};           

private:
  Serc19SimAnalysis *pAnalysis;   

  G4Material*        Air;
  G4Material*        Brass;
  G4Material*        SiliconStr;
  G4Material*        Scintillator;
  G4Material*        G10;
  G4Material*        CarbonFRP;
  G4Material*        Paraffin;
  G4Material*        PbWO4;
  G4Material*        Iron;

  G4Tubs*            solidWorld;
  G4Tubs*            solidTrack;
  G4Tubs*            solidTrackSpt;
  G4Sphere*          solidEcal;
  G4Sphere*          solidElectronics;
  G4Box*             solidHcal;
  G4Trd*             solidHcalAbs;
  G4Tubs*            solidHcalBox;  
  G4Box*            solidHcalSci; 

  G4LogicalVolume*   logicWorld;
  G4VPhysicalVolume* physiWorld;
  
  G4LogicalVolume*   logicTrack;
  G4VPhysicalVolume* physiTrack;

  G4LogicalVolume*   logicTrackSpt;
  G4VPhysicalVolume* physiTrackSpt;

  G4LogicalVolume*   logicEcal;
  G4VPhysicalVolume* physiEcal;

  G4LogicalVolume*   logicElectronics;
  G4VPhysicalVolume* physiElectronics;

  G4LogicalVolume*   logicHcal;
  G4VPhysicalVolume* physiHcal;

  G4LogicalVolume*   logicHcalAbs;
  G4VPhysicalVolume* physiHcalAbs;

  G4LogicalVolume*   logicHcalBox;
  G4VPhysicalVolume* physiHcalBox;

  G4LogicalVolume*   logicHcalSci;
  G4VPhysicalVolume* physiHcalSci;

  G4VPhysicalVolume* physiHcalSci_div;
  G4VPhysicalVolume* physiParf;


  G4UniformMagField* magField;      //pointer to the magnetic field
  G4SDManager* SDman;
  G4String  ecalName;
  Serc19TrkSD * TrkSD;
  G4Region* aRegion0;

  Serc19EcalSD * EcalSD;
  G4Region* aRegion1;

  Serc19HclSD * HclSD;
  G4Region* aRegion2;
  
  G4VisAttributes* visWorld;
  G4VisAttributes* visTrack;
  G4VisAttributes* visTrackSpt;
  G4VisAttributes* visEcal;
  G4VisAttributes* visElectronics;
  G4VisAttributes* visHcal;
  G4VisAttributes* visHcalAbs;
  G4VisAttributes* visHcalBox;
  G4VisAttributes* visHcalSci;
  G4VisAttributes* visParf;
  G4VisAttributes* visNull;
protected:
  G4FieldManager*         fieldMgr;
  G4FieldManager*         localfldMgr;
private:

  G4VPhysicalVolume* ConstructCalorimeter();     

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

