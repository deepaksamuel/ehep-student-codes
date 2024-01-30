//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR"); // the world volume is Air
  G4Material* box_mat = nist->FindOrBuildMaterial("G4_WATER");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;


  // construct a world volume box of dimensions 30 cm x 30 cm x 30 cm
  G4Box* solidWorld =  new G4Box("World", 15*cm, 15*cm, 15*cm);     //Note convention: Half length is used
  G4LogicalVolume* logicWorld =   new G4LogicalVolume(solidWorld, world_mat,"World");            //its name
  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking


  // Big box: construct a water box of dimensions 20 cm x 20 cm x 20 cm
  G4Box* solidBox =  new G4Box("Box", 10*cm, 10*cm, 10*cm);     //Note convention: Half length is used
  G4LogicalVolume* logicBox =   new G4LogicalVolume(solidBox, box_mat,"Box");            //its name
  G4RotationMatrix* rotation = new G4RotationMatrix();
  rotation->rotateX(0*deg);
  rotation->rotateY(0*deg);
  rotation->rotateZ(0*deg);
  
  G4VPhysicalVolume* physBox =
                      new G4PVPlacement(rotation,     // rotation matrix
                      G4ThreeVector(0,0,0),      //at // translation
                      logicBox,             //its logical volume
                      "Box",                 //its name
                      logicWorld,            //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking



  // Second smaller box
  // G4Box* solidBox2 =  new G4Box("Box", 2*cm, 2*cm, 2*cm);     //Note convention: Half length is used
  // G4LogicalVolume* logicBox2 =   new G4LogicalVolume(solidBox2, box_mat,"Box");            //its name
  
  // G4VPhysicalVolume* physBox2 =
  //                     new G4PVPlacement(0,     // no rotation 
  //                     G4ThreeVector(0,0,0),      //at // translation
  //                     logicBox2,             //its logical volume
  //                     "Box",                 //its name
  //                     logicWorld,            //its mother  volume
  //                     false,                 //no boolean operation
  //                     0,                     //copy number
  //                     checkOverlaps);        //overlaps checking



  // Set Shape2 as scoring volume
  //
  fScoringVolume = logicBox;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
