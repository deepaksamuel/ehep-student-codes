//////////////////////////////////
//
//  Define "uniformField" to have uniform field and use field value
//  in SetUniformMagField(1*tesla);<-
//  and choose field fdirection in 
//  magField = new G4UniformMagField(G4ThreeVector(fieldValue,0., 0.)); //Field along x-axis
//
///////////////////////
//
// Define "arbField" for arbtraty magnetic field map
// give input field value as a function of x in file_magField.dat, 
// which is read in constructor
// Serc19MagneticField::Serc19MagneticField(const G4String &file_magField)
// Define proper filed value in memberfunction
//void Serc19MagneticField::MagneticField(const double x[3], double B[3]) const
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//#define uniformField
#define arbField 

#include "Serc19DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Trd.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVDivision.hh"

#ifdef uniformField
#include "G4UniformMagField.hh"
#endif

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4NistManager.hh"
#include "G4PVParameterised.hh"
#include "G4RotationMatrix.hh"

#include "G4UniformMagField.hh"

#include "G4FieldManager.hh"
#include "Serc19Field.hh"
#define debug

#ifdef debug
#include "G4Timer.hh"
#endif

#include "G4ProductionCuts.hh"
#include "G4MaterialCutsCouple.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Serc19DetectorConstruction::Serc19DetectorConstruction(Serc19SimAnalysis *p)
  : pAnalysis(p),
    Air(0),
    Brass(0),
    SiliconStr(0),
    Scintillator(0),
    G10(0),
    CarbonFRP(0),
    Paraffin(0),
    PbWO4(0),    

    solidWorld(0),
    solidTrack(0),
    solidTrackSpt(0),
    solidEcal(0),
    solidElectronics(0),
    solidHcal(0),
    solidHcalAbs(0),
    solidHcalBox(0),
    solidHcalSci(0),

    logicWorld(0), physiWorld(0),
    logicTrack(0), physiTrack(0),
    logicTrackSpt(0), physiTrackSpt(0),
    logicEcal(0), physiEcal(0),
    logicElectronics(0), physiElectronics(0),
    logicHcal(0), physiHcal(0),
    logicHcalAbs(0), physiHcalAbs(0),
    logicHcalBox(0), physiHcalBox(0),
    logicHcalSci(0), physiHcalSci(0),
    physiHcalSci_div(0), physiParf(0),
    magField(0), 
    visWorld(0),
    visTrack(0),
    visTrackSpt(0),
    visEcal(0),
    visElectronics(0),
    visHcal(0),
    visHcalAbs(0),
    visHcalBox(0),
    visHcalSci(0),
    visParf(0),
    visNull(0),
    fieldMgr(0),
    localfldMgr(0)

{

  // materials
  DefineMaterials();

  SDman = G4SDManager::GetSDMpointer();
  ecalName = "serc19Trk";
  TrkSD = new Serc19TrkSD(ecalName);
  SDman->AddNewDetector(TrkSD);

  aRegion0 = new G4Region("Tracker_block");

  ecalName = "serc19Ecal";
  EcalSD = new Serc19EcalSD(ecalName);
  SDman->AddNewDetector(EcalSD);

  aRegion1 = new G4Region("ECAL_Block");

  ecalName = "serc19Hcl";
  HclSD = new Serc19HclSD(ecalName);
  SDman->AddNewDetector(HclSD);

  aRegion2 = new G4Region("HCAL_Block");


  // create commands for interactive definition of the calorimeter

  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Serc19DetectorConstruction::~Serc19DetectorConstruction() {

//    if (  Air) {                cout<<"1xxx"<<endl; delete  Air;}            
//    if (  Brass) { 	       cout<<"xxx"<<endl; delete  Brass;}          
//    if (  SiliconStr) { 	       cout<<"xxx"<<endl; delete  SiliconStr;}     
//    if (  Scintillator) {       cout<<"xxx"<<endl; delete  Scintillator;}   
//    if (  G10) { 	       cout<<"xxx"<<endl; delete  G10;}            
//    if (  CarbonFRP) { 	       cout<<"xxx"<<endl; delete  CarbonFRP;}      
//    if (  Paraffin) { 	       cout<<"xxx"<<endl; delete  Paraffin;}       
//    if (  PbWO4) { 	       cout<<"xxx"<<endl; delete  PbWO4;}          
			                               
//  if ( solidWorld) { 	     cout<<"2xxx"<<endl; delete solidWorld;}        
//  if ( solidTrack) { 	     cout<<"xxx"<<endl; delete solidTrack;}        
//  if ( solidTrackSpt) { 	     cout<<"xxx"<<endl; delete solidTrackSpt;}     
//  if ( solidEcal) { 	     cout<<"xxx"<<endl; delete solidEcal;}         
//  if ( solidElectronics) {    cout<<"xxx"<<endl; delete solidElectronics;}  
//  if ( solidHcal) { 	     cout<<"xxx"<<endl; delete solidHcal;}         
//  if ( solidHcalAbs) { 	     cout<<"xxx"<<endl; delete solidHcalAbs;}      
//  if ( solidHcalBox) { 	     cout<<"xxx"<<endl; delete solidHcalBox;}      
//  if ( solidHcalSci) { 	     cout<<"xxx"<<endl; delete solidHcalSci;}      
			                               
//  if ( logicWorld) { 	     cout<<"3xxx"<<endl; delete logicWorld;}  
//  cout<<"8xxxxx"<<endl;
//  if ( physiWorld) { 	     cout<<"xxx"<<endl; delete physiWorld;}        
//  cout<<"7xxxxx"<<endl;
//  if ( logicTrack) { 	     cout<<"xxx"<<endl; delete logicTrack;}        
//  cout<<"6xxxxx"<<endl;
//  if ( physiTrack) { 	     cout<<"yyxxx"<<endl; delete physiTrack;}        
//  cout<<"4xxxxx"<<endl;			                               
//  if ( logicTrackSpt) { 	     cout<<"xxx"<<endl; delete logicTrackSpt;}     
// if ( physiTrackSpt) { 	     cout<<"xxx"<<endl; delete physiTrackSpt;}     
			                               
//  if ( logicElectronics) {    cout<<"4xxx"<<endl; delete logicElectronics;}  
//  if ( physiElectronics) {    cout<<"xxx"<<endl; delete physiElectronics;}  
//  cout<<"4xxxxx"<<endl;			                               
//  if ( logicEcal) { 	     cout<<"xxx"<<endl; delete logicEcal;}         
//  if ( physiEcal) { 	     cout<<"xxx"<<endl; delete physiEcal;}         
			                               
//  if ( logicHcal) { 	     cout<<"xxx"<<endl; delete logicHcal;}         
//  if ( physiHcal) { 	     cout<<"xxx"<<endl; delete physiHcal;}         
//  cout<<"3xxxxx"<<endl;			                               
//  if ( logicHcalAbs) { 	     cout<<"xxx"<<endl; delete logicHcalAbs;}      
//  if ( physiHcalAbs) { 	     cout<<"xxx"<<endl; delete physiHcalAbs;}      
			                               
//  if ( logicHcalBox) { 	     cout<<"5xxx"<<endl; delete logicHcalBox;}      
//  if ( physiHcalBox) { 	     cout<<"xxx"<<endl; delete physiHcalBox;}      
			                               
//  if ( logicHcalSci) { 	     cout<<"xxx"<<endl; delete logicHcalSci;}      
//  if ( physiHcalSci) { 	     cout<<"xxx"<<endl; delete physiHcalSci;}      
			                               
//  if ( physiHcalSci_div) {    cout<<"xxx"<<endl; delete physiHcalSci_div;}  
//  if ( physiParf) { 	     cout<<"xxx"<<endl; delete physiParf;}         
//  cout<<"2xxxxx"<<endl;			                               
//   if (  visWorld) { 	      cout<<"6xxx"<<endl; delete  visWorld;}        
//   if (  visTrack) { 	      cout<<"xxx"<<endl; delete  visTrack;}        
//   if (  visTrackSpt) { 	      cout<<"xxx"<<endl; delete  visTrackSpt;}     
//   if (  visEcal) { 	      cout<<"xxx"<<endl; delete  visEcal;}         
//   if (  visElectronics) {     cout<<"xxx"<<endl; delete  visElectronics;}  
//   if (  visHcal) { 	      cout<<"xxx"<<endl; delete  visHcal;}         
//   if (  visHcalAbs) { 	      cout<<"xxx"<<endl; delete  visHcalAbs;}      
//   if (  visHcalBox) { 	      cout<<"xxx"<<endl; delete  visHcalBox;}      
//   if (  visHcalSci) { 	      cout<<"xxx"<<endl; delete  visHcalSci;}      
//   if (  visParf) { 	      cout<<"xxx"<<endl; delete  visParf;}         
//   if (  visNull) {            cout<<"xxx"<<endl; delete  visNull;}   
     
  cout<<"xxxxx"<<endl;
  //delete   cal0SD ; cal0SD=0;
}
  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Serc19DetectorConstruction::Construct() {
  
  // SetUniformMagField(3.8*tesla);
  
  // physiWorld = ConstructCalorimeter();

  // return physiWorld;


  //return ConstructCalorimeter();
// Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  auto solidWorld = new G4Box("World",                           // its name
    0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
    world_mat,                                       // its material
    "World");                                        // its name

  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                           // at (0,0,0)
    logicWorld,                                // its logical volume
    "World",                                   // its name
    nullptr,                                   // its mother  volume
    false,                                     // no boolean operation
    0,                                         // copy number
    checkOverlaps);                            // overlaps checking

  //
  // Envelope
  //
  auto solidEnv = new G4Box("Envelope",                    // its name
    0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ);  // its size

  auto logicEnv = new G4LogicalVolume(solidEnv,  // its solid
    env_mat,                                     // its material
    "Envelope");                                 // its name

  new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),          // at (0,0,0)
    logicEnv,                 // its logical volume
    "Envelope",               // its name
    logicWorld,               // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking

  //
  // Shape 1
  //
  G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  G4ThreeVector pos1 = G4ThreeVector(0, 2*cm, -7*cm);

  // Conical section shape
  G4double shape1_rmina =  0.*cm, shape1_rmaxa = 2.*cm;
  G4double shape1_rminb =  0.*cm, shape1_rmaxb = 4.*cm;
  G4double shape1_hz = 3.*cm;
  G4double shape1_phimin = 0.*deg, shape1_phimax = 360.*deg;
  auto solidShape1 = new G4Cons("Shape1", shape1_rmina, shape1_rmaxa, shape1_rminb, shape1_rmaxb,
    shape1_hz, shape1_phimin, shape1_phimax);

  auto logicShape1 = new G4LogicalVolume(solidShape1,  // its solid
    shape1_mat,                                        // its material
    "Shape1");                                         // its name

  new G4PVPlacement(nullptr,  // no rotation
    pos1,                     // at position
    logicShape1,              // its logical volume
    "Shape1",                 // its name
    logicEnv,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking

  //
  // Shape 2
  //
  G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  G4ThreeVector pos2 = G4ThreeVector(0, -1*cm, 7*cm);

  // Trapezoid shape
  G4double shape2_dxa = 12*cm, shape2_dxb = 12*cm;
  G4double shape2_dya = 10*cm, shape2_dyb = 16*cm;
  G4double shape2_dz  = 6*cm;
  auto solidShape2 = new G4Trd("Shape2",  // its name
    0.5 * shape2_dxa, 0.5 * shape2_dxb, 0.5 * shape2_dya, 0.5 * shape2_dyb,
    0.5 * shape2_dz);  // its size

  auto logicShape2 = new G4LogicalVolume(solidShape2,  // its solid
    shape2_mat,                                        // its material
    "Shape2");                                         // its name

  new G4PVPlacement(nullptr,  // no rotation
    pos2,                     // at position
    logicShape2,              // its logical volume
    "Shape2",                 // its name
    logicEnv,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking

  // Set Shape2 as scoring volume
  //
  // fScoringVolume = logicShape2;

  //
  //always return the physical World
  //
  return physWorld;



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Serc19DetectorConstruction::DefineMaterials()
{ 

	//This function illustrates the possible ways to define materials
  //  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4String name, symbol;             //a=mass of a mole;
  G4double density;      //z=mean number of protons;  
  // n=number of nucleons in an isotope;
  
  G4int ncomponents,natoms;
  G4double fractionmass;
  
  G4NistManager* mat = G4NistManager::Instance();  
  mat->SetVerbose(1);

  // define Elements
  G4Element* C = mat->FindOrBuildElement("C");
  G4Element* H = mat->FindOrBuildElement("H");
  G4Element* O = mat->FindOrBuildElement("O");
  G4Element* Si = mat->FindOrBuildElement("Si");
  G4Element* Cu = mat->FindOrBuildElement("Cu");
  G4Element* Zn = mat->FindOrBuildElement("Zn");
  
  Iron =  mat->FindOrBuildMaterial("G4_Fe");
  Air =  mat->FindOrBuildMaterial("G4_AIR");

  Brass =  new G4Material("Brass", density=8.4*g/cm3, ncomponents=2);
  Brass->AddElement(Cu, fractionmass=70.0*perCent);
  Brass->AddElement(Zn, fractionmass=30.0*perCent);  
  
  SiliconStr = new G4Material("Silicon", density=2.32*g/cm3, ncomponents=2);
  SiliconStr->AddElement(Si, fractionmass=99.99*perCent);
  SiliconStr->AddElement(O, fractionmass=0.01*perCent);

  Scintillator = new G4Material(name="Scintillator", density=1.032*g/cm3, ncomponents=2);
  Scintillator->AddElement(C, 9);
  Scintillator->AddElement(H, 10);

  // //  Scintillator = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  G10 = new G4Material("G10", density= 1.09*g/cm3, ncomponents=4);
  
  G10->AddElement(Si, natoms=1);
  G10->AddElement(O , natoms=2);
  G10->AddElement(C , natoms=3);
  G10->AddElement(H , natoms=3);
  
  CarbonFRP = new G4Material("FRPCarbon",density = 2.26*g/cm3,ncomponents=2);
  CarbonFRP->AddElement(C, fractionmass=99.*perCent);
  CarbonFRP->AddElement(H, fractionmass= 1.*perCent);

  Paraffin =  new G4Material("Paraffin",density = 0.9*g/cm3, ncomponents=2);
  Paraffin->AddElement(C, natoms=31);
  Paraffin->AddElement(H, natoms=64);
 
  //  G4Element* elO = new G4Element(name="Oxygen",symbol="O",z=8.,16.00*g/mole);
  G4Element* elW = new G4Element(name="Tungsten",symbol="W",z=74.,183.84*g/mole);
  G4Element* elPb = new G4Element(name="Lead",symbol="Pb",z=82.,207.20*g/mole);
  
  PbWO4 = new G4Material(name="PbWO4",density=8.280*g/cm3, ncomponents=3);
  PbWO4->AddElement(elPb, natoms=1); 
  PbWO4->AddElement(elW, natoms=1);
  PbWO4->AddElement(O, natoms=4);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Serc19DetectorConstruction::ConstructCalorimeter()
{
  // cout <<"G4VPhysicalVolume* Serc19DetectorConstruction::ConstructCalorimeter()"<<endl;
  // Clean old geometry, if any//
  
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();


  //13 layer silicon detector with silicon support structure
  const int nsilayer=pAnalysis->nsilayer;
  //  double sirad[pAnalysis->nsilayer]={2.9, 6.8, 10.9, 16.0, 23, 32, 41, 50, 60, 70, 80, 90, 100};
  double sirad[pAnalysis->nsilayer]={2.9, 6.8, 10.9, 16.0, 22, 29, 37, 45, 52, 59, 66, 73, 80}; 

  double sithickness = 0.3*mm;
  double supportwidth=0.5*mm;
  double epsilon=1.e-6*mm;
  for (int ij=0; ij<nsilayer; ij++) { sirad[ij] *=cm;}

  double pival = acos(-1)*rad;
  //     
  // World
  //
  //  solidWorld = 0; logicWorld=0; physiWorld=0;
  solidWorld = new G4Tubs("WorldSolid", 0.0, 3.2*m, 4.5*m, -pival/2.-5*epsilon, pival+10*epsilon);
  
  logicWorld = new G4LogicalVolume(solidWorld,
																	 Air,
																	 "WorldLogic");
  physiWorld = new G4PVPlacement(0,
																 G4ThreeVector(0,0,0),
																 logicWorld,
																 "WorldPhysi",
																 0,
																 false,
																 0);
  
	
  //  G4Tubs* solidTrackx[pAnalysis->nsilayer]={0};
  //  G4LogicalVolume* logicTrackx[pAnalysis->nsilayer]={0};
  //  G4VPhysicalVolume* physiTrackx[pAnalysis->nsilayer]={0};
	
  for (int ij=0; ij<nsilayer; ij++) { 
    solidTrack = new G4Tubs("Silicontrk", sirad[ij], sirad[ij]+sithickness, 75*cm, -pival/4., pival/2.);
    logicTrack  = new G4LogicalVolume(solidTrack, SiliconStr, "SitrkLog");
    physiTrack  = new G4PVPlacement(0, 
																		G4ThreeVector(0,0,0), 
																		logicTrack,
																		"solidTrack",
																		logicWorld,
																		false,
																		ij);
    logicTrack->SetSensitiveDetector(TrkSD);
		logicTrack->SetRegion(aRegion0);
		aRegion0->AddRootLogicalVolume(logicTrack);
		
	}
	
	
  //  G4Tubs* solidTrackSptx[pAnalysis->nsilayer]={0};
  //  G4LogicalVolume* logicTrackSptx[pAnalysis->nsilayer]={0};
  //  G4VPhysicalVolume* physiTrackSptx[pAnalysis->nsilayer]={0};
  
  for (int ij=0; ij<nsilayer; ij++) { 
    //    solidTrackSpt = new G4Tubs("TrackSpt", sirad[ij], sirad[ij]+sithickness, 80*cm, -pival/4., pival/2);
    //    logicTrackSpt  = new G4LogicalVolume(solidTrackSpt, SiliconStr, "TrackSptLog");
    solidTrackSpt = new G4Tubs("TrackSpt", sirad[ij]-supportwidth, sirad[ij]-epsilon, 80*cm, -pival/4., pival/2);
    logicTrackSpt  = new G4LogicalVolume(solidTrackSpt, CarbonFRP, "TrackSptLog");
		
		
    physiTrackSpt  = new G4PVPlacement(0, 
																			 G4ThreeVector(0,0,0), 
																			 logicTrackSpt,
																			 "solidTrackSpt",
																			 logicWorld,
																			 false,
																			 nsilayer+ij);
  }
	
	
	//Segmentation fault  
  //Paraffin in between tracker and ECAL
	G4Tubs* solidParf = new G4Tubs("SolParaffin", (sirad[nsilayer-1]+2.0*cm), (sirad[nsilayer-1]+6.0*cm), 1.2*m, -pival/4., pival/2.);
  
  G4LogicalVolume* logicParf = new G4LogicalVolume(solidParf,
																									 Paraffin,
																									 "ParfLogic");
  physiParf = new G4PVPlacement(0,
																G4ThreeVector(0,0,0),
																logicParf,
																"ParfPhysi",
																logicWorld,
																false,
																0);  
  
	//  //ECAL : for simplicity used simple sphreical shell, 
  //       but each crystal in different theta/eta should have different trapizium 
	
  solidEcal = new G4Sphere("SolECAL", 123.0*cm, 152.6*cm, -pival/4., pival/2., pival/4., pival/2.);
  logicEcal = new G4LogicalVolume(solidEcal, PbWO4, "EcalLogic");
  physiEcal = new G4PVPlacement(0,
																G4ThreeVector(0,0,0),
																logicEcal,
																"EcalPhysi",
																logicWorld,
																false,
																0); 
	
	logicEcal->SetSensitiveDetector(EcalSD);
	logicEcal->SetRegion(aRegion1);
	aRegion1->AddRootLogicalVolume(logicEcal);  
  
  //Electronics behind ECAL using G10
	//Segmentation fault
  solidElectronics = new G4Sphere("SolElectronics", 155.0*cm, 160*cm, -pival/4., pival/2., pival/4., pival/2.);
  logicElectronics = new G4LogicalVolume( solidElectronics, G10, "ElectronicsLogic");
  physiElectronics = new G4PVPlacement(0,
																			 G4ThreeVector(0,0,0),
																			 logicElectronics,
																			 "ElectronicsPhysi",
																			 logicWorld,
																			 false,
																			 0);
	
  //HCAL : Make volume of Brass and then placed 4mm scintillator inside
  double r1 = 181.0*cm;  //Inner radius
  double r2 = 295.0*cm;
  double dz = 260.0*cm;
	
  int nhclwedge=pAnalysis->nhclwedge; //36
  int nhcalLayer=pAnalysis->nhcalLayer; //17;
  int nhcalEtaDiv=pAnalysis->nhcalEtaDiv; //34;
	
  double active_dr=0.4*cm; //Width of air box is twice the schintillator width
  double passive_dr=5.0*cm; //Thickness of Brass absorber  
  double wedgeangle=5.0*degree; //Width of scintillator supporting wedge
  double scintangle=4.9*degree; //Width of scintillator in phi
	
  localfldMgr = new G4FieldManager(new Serc19Field());
	
  solidHcalBox = new G4Tubs("SolHcalBox", 0.95*r1, 1.05*r2, 4.3*m, -pival/2.-2*epsilon, pival+4*epsilon);
	logicHcalBox = new G4LogicalVolume( solidHcalBox, Brass, "HcalBoxLogic", localfldMgr);
	//  logicHcalBox = new G4LogicalVolume( solidHcalBox, Brass, "HcalBoxLogic");
	//  logicHcalBox->setFieldManager(localFldMgr, true);
	
	
  physiHcalBox = new G4PVPlacement(0,
																	 G4ThreeVector(0,0,0),
																	 logicHcalBox,
																	 "HcalBoxPhysi",
																	 logicWorld,
																	 false,
																	 0);    
  double angle = atan((dz-5*cm)/r1); //Side is also made with 5cm thick iron 
  double dx1 = r1*tan(angle);
  double dx2 = r2*tan(angle);
  double dy1 = r1*tan(scintangle)/2.;
  double dy2 = r2*tan(scintangle)/2.;
  double dr = (r1+r2)/2.; //Middle of the wedge
	
  solidHcalAbs = new G4Trd("SolHcalAbs", dx1, dx2, dy1, dy2, (r2-r1)/2.);
  logicHcalAbs = new G4LogicalVolume(solidHcalAbs, Brass, "HcalAbsLogic");
	
  for (int ij=0; ij<nhclwedge; ij++) {
    double rotang = -pival/2. + (ij+0.5)*wedgeangle;
    //    for (int ij=0; ij<1; ij++) { 
    //      double rotang = 0; //-pival/2. + (ij+0.5)*wedgeangle;
    
    G4ThreeVector xAxis(0.0, 0.0, -1.0);    
    G4ThreeVector yAxis(-sin(rotang), cos(rotang), 0.0);
    G4ThreeVector zAxis(cos(rotang), sin(rotang), 0.0);    
    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateAxes(xAxis, yAxis, zAxis);
    rot->invert();
    
    physiHcalAbs = new G4PVPlacement(rot,
																		 G4ThreeVector(dr*cos(rotang), dr*sin(rotang),0),
																		 logicHcalAbs,
																		 "HcalAbsPhysi",
																		 logicHcalBox,
																		 false,
																		 ij);   
  }
	
	//   //Put 17 layer 0f scintllator box, which are made of air and scintillator inside the absorber
	
	
	// //   G4Box* solidHcalx[pAnalysis->nhcalLayer];
	// //   G4LogicalVolume* logicHcalx[pAnalysis->nhcalLayer];
	
	// //   G4Box* solidHcalScix[pAnalysis->nhcalLayer];
	// //   G4LogicalVolume* logicHcalScix[pAnalysis->nhcalLayer];
	
	// //   G4Box* solidHcalSci_divx[pAnalysis->nhcalLayer];
	// //   G4LogicalVolume* logicHcalSci_divx[pAnalysis->nhcalLayer];
  
	// //   G4VPhysicalVolume* physiHcalx[pAnalysis->nhcalLayer];
	// //   G4VPhysicalVolume* physiHcalScix[pAnalysis->nhcalLayer];
	// //   G4VPhysicalVolume* physiHcalSci_divx[pAnalysis->nhcalLayer];
	

  for (int ij=0; ij<nhcalLayer; ij++) { 
    //  for (int ij=0; ij<1/*nhcalLayer*/; ij++) { 
		
    double drsc = r1+8.0*cm+passive_dr/2.0+ij*(2*active_dr+passive_dr);
    double dxsc = 0.995*drsc*tan(angle);
    double dysc = 0.95*drsc*tan(scintangle)/2.;
    double tilesize = dxsc/nhcalEtaDiv;
		
    solidHcal = new G4Box("SolHcal", dxsc+epsilon, dysc, active_dr);
    logicHcal = new G4LogicalVolume(solidHcal, Air, "HcalLogic");
		
    solidHcalSci = new G4Box("SolHcalSci", 0.995*dxsc, 0.95*dysc-epsilon, active_dr/2.0);
    logicHcalSci = new G4LogicalVolume(solidHcalSci, Air, "HcalLogicSci");
    
    G4Box* solidHcalSci_div = new G4Box("SolHcalSci_div", tilesize, dysc-epsilon, active_dr/2.0);
    G4LogicalVolume* logicHcalSci_div =  new G4LogicalVolume(solidHcalSci_div, Scintillator, "HcalLogicSci_div");
    //    physiHcalSci_div = new G4PVReplica("PhyHcalSci_div", logicHcalSci_div[ij], logicHcalSci, kXAxis, nhcalEtaDiv, tilesize); 
    physiHcalSci_div = new G4PVDivision("PhyHcalSci_div", logicHcalSci_div, logicHcalSci, kXAxis, nhcalEtaDiv, tilesize);     
    
    physiHcal = new G4PVPlacement(0,
																	G4ThreeVector(0,0,drsc-dr),
																	logicHcal,
																	"HcalPhysi",
																	logicHcalAbs,
																	false,
																	ij);   
    physiHcalSci = new G4PVPlacement(0,
																		 G4ThreeVector(0,0,0),
																		 logicHcalSci,
																		 "HcalPhysi",
																		 logicHcal,
																		 false,
																		 0);     
    logicHcalSci_div->SetVisAttributes(visNull);
		
    logicHcalSci_div->SetSensitiveDetector(HclSD);
		
    //    logicHcalSci_div->SetRegion(aRegion2);
    //    aRegion2->AddRootLogicalVolume(logicHcalSci_div);  
  }
  
//   // SDman->AddNewDetector(cal0SD);
  


  

	//   if(!visWorld){visWorld= new G4VisAttributes(true, G4Colour(0.0,0.0,1.0));}//blue
	if(!visTrack){ visTrack= new G4VisAttributes(true, G4Colour(1.0,0.0,0.0));}//red
	visTrack->SetForceAuxEdgeVisible(true);
  
	//   if(!visTrackSpt){ visTrackSpt= new G4VisAttributes(true, G4Colour(1.0,1.0,0.0));}//yellow
	//   //visTrackSpt->SetForceSolid(true);
	if(!visEcal){
		visEcal= new G4VisAttributes(true, G4Colour(0.5,0.5,0.5));
		visEcal->SetForceSolid(true);
	}//gray
	//   if(!visElectronics){visElectronics= new G4VisAttributes(true, G4Colour(0.0,1.0,1.0));}//cyan
	//   if(!visHcal){visHcal= new G4VisAttributes(true, G4Colour(1.0,0.0,1.0));}//cyan
	if(!visHcalAbs){
		visHcalAbs= new G4VisAttributes(true, G4Colour(0.0,1.0,0.0));
		visHcalAbs->SetForceWireframe(1);
		//      visHcalAbs->SetForceSolid(true);
	}//green
	if(!visHcalBox){visHcalBox= new G4VisAttributes(true, G4Colour(0.0,0.0,1.0));}//blue
	//    visHcalBox->SetForceSolid(true);
  
	//   if(!visParf){visParf= new G4VisAttributes(true, G4Colour(1.0,1.0,1.0));}//white
	
	//   //  logicWorld->SetVisAttributes(visWorld);
	//   logicWorld->SetVisAttributes(visNull);
  
	logicTrack->SetVisAttributes(visTrack);
	//   //   logicTrackSpt->SetVisAttributes(visTrackSpt);
	//  logicTrackSpt->SetVisAttributes(visNull);
	logicEcal->SetVisAttributes(visEcal);
	//   //  logicEcal->SetVisAttributes(visNull);
	//   logicElectronics->SetVisAttributes(visNull);
	logicElectronics->SetVisAttributes(visElectronics);
  
	//   logicHcal->SetVisAttributes(visNull);
	logicHcal->SetVisAttributes(visHcal);
	logicHcalAbs->SetVisAttributes(visHcalAbs);
	//    logicHcalAbs->SetVisAttributes(visNull);
	logicHcalBox->SetVisAttributes(visHcalBox);
	//       logicHcalBox->SetVisAttributes(visNull);
	logicHcalSci->SetVisAttributes(visHcalSci);
	//    logicHcalSci->SetVisAttributes(visNull);
	
	
	logicParf->SetVisAttributes(visParf);
	//   logicParf->SetVisAttributes(visNull);
  
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

void Serc19DetectorConstruction::SetUniformMagField(G4double fieldValue) {
  //apply a global uniform magnetic field along Z axis

  fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  cout<<"fieldvalue "<< fieldValue/tesla<<" Tesla"<<endl;
  if(magField) { delete magField;}		//delete the existing magn field

  magField = new G4UniformMagField(G4ThreeVector(0.0, 0.0, fieldValue)); //Field along z-axis
  fieldMgr->SetDetectorField(magField);
  fieldMgr->CreateChordFinder(magField);

}

//void Serc19DetectorConstruction::UpdateGeometry()
//{
//  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

/*
*** glibc detected *** bin/inoical0_field: double free or corruption (out): 0x00000000014ad6b0 ***
======= Backtrace: =========
/lib64/libc.so.6[0x315bc75f4e]
/lib64/libc.so.6[0x315bc78cf0]
bin/inoical0_field(_ZN25Serc19DetectorConstructionD1Ev+0x3e6)[0x433d16]
bin/inoical0_field(_ZN25Serc19DetectorConstructionD0Ev+0x18)[0x4341ca]
/usr/local/physics/geant4.10.00.p01-install/lib64/libG4run.so(_ZN12G4RunManager25DeleteUserInitializationsEv+0x13)[0x7f0a0c064c93]
/usr/local/physics/geant4.10.00.p01-install/lib64/libG4run.so(_ZN12G4RunManagerD1Ev+0xd1)[0x7f0a0c0662b1]
/usr/local/physics/geant4.10.00.p01-install/lib64/libG4run.so(_ZN12G4RunManagerD0Ev+0x9)[0x7f0a0c066539]
bin/inoical0_field(main+0x8cf)[0x42b55f]
/lib64/libc.so.6(__libc_start_main+0xfd)[0x315bc1ed5d]
bin/inoical0_field[0x41e189]

//    logicHcalSci_div->SetRegion(aRegion2);
//    aRegion2->AddRootLogicalVolume(logicHcalSci_div);  
#2  0x00007f9dd3dba4e8 in TUnixSystem::StackTrace() () from /usr/local/lib/libCore.so
#3  0x00007f9dd3db9363 in TUnixSystem::DispatchSignals(ESignals) () from /usr/local/lib/libCore.so
#4  <signal handler called>
#5  0x00007f9dce4e5932 in G4VEmProcess::PostStepGetPhysicalInteractionLength(G4Track const&, double, G4ForceCondition*) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4processes.so
#6  0x00007f9dcf59d638 in G4SteppingManager::DefinePhysicalStepLength() () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4tracking.so

#0  0x000000315bcac65e in waitpid () from /lib64/libc.so.6
#1  0x000000315bc3e609 in do_system () from /lib64/libc.so.6
#2  0x00007fb46ab664e8 in TUnixSystem::StackTrace() () from /usr/local/lib/libCore.so
#3  0x00007fb46ab65363 in TUnixSystem::DispatchSignals(ESignals) () from /usr/local/lib/libCore.so
#4  <signal handler called>
#5  0x00007fb4640e1e59 in G4PhysicalVolumeStore::DeRegister(G4VPhysicalVolume*) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4geometry.so
#6  0x00007fb4640ee60b in G4VPhysicalVolume::~G4VPhysicalVolume() () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4geometry.so
#7  0x00007fb464266dd9 in G4PVPlacement::~G4PVPlacement() () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4geometry.so
#8  0x000000000042bede in Serc19DetectorConstruction::~Serc19DetectorConstruction (this=0x2818660, __in_chrg=<value optimized out>) at src/Serc19DetectorConstruction.cc:160 //delete physiTrack;
#9  0x000000000042c392 in Serc19DetectorConstruction::~Serc19DetectorConstruction (this=0x2818660, __in_chrg=<value optimized out>) at src/Serc19DetectorConstruction.cc:199 //end of detructor
#10 0x00000000004247b6 in main (argc=2, argv=0x7ffd712bda88) at src/inoical0_field.cc:138


===========================================================
#0  0x000000315bcac65e in waitpid () from /lib64/libc.so.6
#1  0x000000315bc3e609 in do_system () from /lib64/libc.so.6
#2  0x00007fbeb830f4e8 in TUnixSystem::StackTrace() () from /usr/local/lib/libCore.so
#3  0x00007fbeb830e363 in TUnixSystem::DispatchSignals(ESignals) () from /usr/local/lib/libCore.so
#4  <signal handler called>
#5  0x00007fbeb18d7379 in G4VoxelNavigation::ComputeStep(CLHEP::Hep3Vector const&, CLHEP::Hep3Vector const&, double, double&, G4NavigationHistory&, bool&, CLHEP::Hep3Vector&, bool&, bool&, G4VPhysicalVolume**, int&) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4geometry.so
#6  0x00007fbeb18b8554 in G4Navigator::ComputeStep(CLHEP::Hep3Vector const&, CLHEP::Hep3Vector const&, double, double&) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4geometry.so
#7  0x00007fbeb301128d in G4Transportation::AlongStepGetPhysicalInteractionLength(G4Track const&, double, double, double&, G4GPILSelection*) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4processes.so
#8  0x00007fbeb3af27ef in G4SteppingManager::DefinePhysicalStepLength() () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4tracking.so
#9  0x00007fbeb3af0f88 in G4SteppingManager::Stepping() () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4tracking.so
#10 0x00007fbeb3af99cc in G4TrackingManager::ProcessOneTrack(G4Track*) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4tracking.so
#11 0x00007fbeb3d279ca in G4EventManager::DoProcessing(G4Event*) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4event.so
#12 0x00007fbeb3f96497 in G4RunManager::ProcessOneEvent(int) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4run.so
#13 0x00007fbeb3f95c15 in G4RunManager::DoEventLoop(int, char const*, int) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4run.so
#14 0x00007fbeb3f95bb2 in G4RunManager::BeamOn(int, char const*, int) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4run.so
#15 0x00007fbeb3fa9925 in G4RunMessenger::SetNewValue(G4UIcommand*, G4String) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4run.so
#16 0x00007fbeb10b67ab in G4UIcommand::DoIt(G4String) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4intercoms.so
#17 0x00007fbeb10c2e97 in G4UImanager::ApplyCommand(char const*) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4intercoms.so
#18 0x00007fbeb10a80b7 in G4UIbatch::ExecCommand(G4String const&) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4intercoms.so
#19 0x00007fbeb10a904b in G4UIbatch::SessionStart() () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4intercoms.so
#20 0x00007fbeb10c0a23 in G4UImanager::ExecuteMacroFile(char const*) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4intercoms.so
#21 0x00007fbeb10bcf1f in G4UIcontrolMessenger::SetNewValue(G4UIcommand*, G4String) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4intercoms.so
#22 0x00007fbeb10b67ab in G4UIcommand::DoIt(G4String) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4intercoms.so
#23 0x00007fbeb10c2e97 in G4UImanager::ApplyCommand(char const*) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4intercoms.so
#24 0x0000000000422fc7 in main (argc=2, argv=0x7ffe2692e1e8) at src/serc19_cmsmodel.cc:123
===========================================================


The lines below might hint at the cause of the crash.
If they do not help you then please submit a bug report at
http://root.cern.ch/bugs. Please post the ENTIRE stack trace
from above as an attachment in addition to anything else
that might help us fixing this issue.
===========================================================
#5  0x00007fbeb18d7379 in G4VoxelNavigation::ComputeStep(CLHEP::Hep3Vector const&, CLHEP::Hep3Vector const&, double, double&, G4NavigationHistory&, bool&, CLHEP::Hep3Vector&, bool&, bool&, G4VPhysicalVolume**, int&) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4geometry.so
#6  0x00007fbeb18b8554 in G4Navigator::ComputeStep(CLHEP::Hep3Vector const&, CLHEP::Hep3Vector const&, double, double&) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4geometry.so
#7  0x00007fbeb301128d in G4Transportation::AlongStepGetPhysicalInteractionLength(G4Track const&, double, double, double&, G4GPILSelection*) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4processes.so
#8  0x00007fbeb3af27ef in G4SteppingManager::DefinePhysicalStepLength() () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4tracking.so
#9  0x00007fbeb3af0f88 in G4SteppingManager::Stepping() () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4tracking.so
#10 0x00007fbeb3af99cc in G4TrackingManager::ProcessOneTrack(G4Track*) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4tracking.so
#11 0x00007fbeb3d279ca in G4EventManager::DoProcessing(G4Event*) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4event.so
#12 0x00007fbeb3f96497 in G4RunManager::ProcessOneEvent(int) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4run.so
#13 0x00007fbeb3f95c15 in G4RunManager::DoEventLoop(int, char const*, int) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4run.so
#14 0x00007fbeb3f95bb2 in G4RunManager::BeamOn(int, char const*, int) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4run.so
#15 0x00007fbeb3fa9925 in G4RunMessenger::SetNewValue(G4UIcommand*, G4String) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4run.so
#16 0x00007fbeb10b67ab in G4UIcommand::DoIt(G4String) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4intercoms.so
#17 0x00007fbeb10c2e97 in G4UImanager::ApplyCommand(char const*) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4intercoms.so
#18 0x00007fbeb10a80b7 in G4UIbatch::ExecCommand(G4String const&) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4intercoms.so
#19 0x00007fbeb10a904b in G4UIbatch::SessionStart() () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4intercoms.so
#20 0x00007fbeb10c0a23 in G4UImanager::ExecuteMacroFile(char const*) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4intercoms.so
#21 0x00007fbeb10bcf1f in G4UIcontrolMessenger::SetNewValue(G4UIcommand*, G4String) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4intercoms.so
#22 0x00007fbeb10b67ab in G4UIcommand::DoIt(G4String) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4intercoms.so
#23 0x00007fbeb10c2e97 in G4UImanager::ApplyCommand(char const*) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4intercoms.so
#24 0x0000000000422fc7 in main (argc=2, argv=0x7ffe2692e1e8) at src/serc19_cmsmodel.cc:123
===========================================================


#0  0x000000315bcac65e in waitpid () from /lib64/libc.so.6
#1  0x000000315bc3e609 in do_system () from /lib64/libc.so.6
#2  0x00007f8fd03be4e8 in TUnixSystem::StackTrace() () from /usr/local/lib/libCore.so
#3  0x00007f8fd03bd363 in TUnixSystem::DispatchSignals(ESignals) () from /usr/local/lib/libCore.so
#4  <signal handler called>
#5  0x00007f8fc9939e59 in G4PhysicalVolumeStore::DeRegister(G4VPhysicalVolume*) () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4geometry.so
#6  0x00007f8fc994660b in G4VPhysicalVolume::~G4VPhysicalVolume() () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4geometry.so
#7  0x00007f8fc9abedd9 in G4PVPlacement::~G4PVPlacement() () from /usr/local/physics/geant4.10.00.p01-install/lib64/libG4geometry.so
#8  0x000000000041ab72 in Serc19DetectorConstruction::~Serc19DetectorConstruction (this=0x1c084f0, __in_chrg=<value optimized out>) at src/Serc19DetectorConstruction.cc:172  if ( physiTrackSpt) { 	     cout<<"xxx"<<endl; delete physiTrackSpt;} 
#9  0x000000000041b438 in Serc19DetectorConstruction::~Serc19DetectorConstruction (this=0x1c084f0, __in_chrg=<value optimized out>) at src/Serc19DetectorConstruction.cc:208
#10 0x00000000004201ce in main (argc=3, argv=0x7fffdc68e298) at src/serc19_cmsmodel.cc:139

WARNING - Attempt to delete the solid store while geometry closed !
WARNING - Attempt to delete the logical volume store while geometry closed !
WARNING - Attempt to delete the physical volume store while geometry closed !
WARNING - Attempt to delete the region store while geometry closed !





*/
