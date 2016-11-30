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
// $Id: DetectorConstruction.cc 87359 2014-12-01 16:04:27Z gcosmo $
// 
/// \file B4DetectorConstruction.cc
/// \brief Implementation of the B4DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
   //fAbsorberPV(0),
   //fGapPV(0),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  G4NistManager* nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_WATER");
  nistManager->FindOrBuildMaterial("G4_AIR");
  
  G4String name ,symbol;
  G4int n, z, ncomponents; 
  G4double a, density, fractionmass;
  
  G4Isotope * Gd156_iso = new G4Isotope(name="Gd156_iso",z=64,n=156,a=155.922*g/mole);
  G4Element * Gd156_el = new G4Element(name="Gd156_el",symbol="Gd156",ncomponents=1);
  Gd156_el->AddIsotope(Gd156_iso,fractionmass=100.*perCent);
  G4Material * Gd156_mat = new G4Material(name="Gd156_mat",density=7.90*g/cm3,ncomponents=1);
  Gd156_mat->AddElement(Gd156_el,fractionmass=100.*perCent);

  G4Isotope * C12_iso = new G4Isotope(name="C12_iso",z=6,n=6,a=12.*g/mole);
  G4Element * C12_el = new G4Element(name="C12_el",symbol="C12",ncomponents=1);
  C12_el->AddIsotope(C12_iso,fractionmass=100.*perCent);
  G4Material * C12_mat = new G4Material(name="C12_mat",density=7.90*g/cm3,ncomponents=1);
  C12_mat->AddElement(C12_el,fractionmass=100.*perCent);

  G4Isotope * Tc98_iso = new G4Isotope(name="Tc98_iso",z=43,n=98,a=98.*g/mole);
  G4Element * Tc98_el = new G4Element(name="Tc98_el",symbol="Tc98",ncomponents=1);
  Tc98_el->AddIsotope(Tc98_iso,fractionmass=100.*perCent);
  G4Material * Tc98_mat = new G4Material(name="Tc98_mat",density=11.5*g/cm3,ncomponents=1);
  Tc98_mat->AddElement(Tc98_el,fractionmass=100.*perCent);
  
  G4Isotope * Eu152_iso = new G4Isotope(name="Eu152_iso",z=63,n=152,a=152.*g/mole);
  G4Element * Eu152_el = new G4Element(name="Eu152_el",symbol="Eu152",ncomponents=1);
  Eu152_el->AddIsotope(Eu152_iso,fractionmass=100.*perCent);
  G4Material * Eu152_mat = new G4Material(name="Eu152_mat",density=7.90*g/cm3,ncomponents=1);
  Eu152_mat->AddElement(Eu152_el,fractionmass=100.*perCent);
  
  // Liquid argon material
  //new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Vacuum
  //new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
  //                kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  G4double boxSizeXY = 10.*cm;
  G4double boxSizeZ = 30.*cm;
  
  G4double targetSizeXY = 5.*mm;
  G4double targetSizeZ = 10.*mm;
  
  G4double worldSizeXYZ = 1.2 * boxSizeZ;

  // Get materials
  G4Material* airMaterial = G4Material::GetMaterial("G4_AIR");
  G4Material* boxMaterial = G4Material::GetMaterial("G4_WATER");
  G4Material* tracerMaterial = G4Material::GetMaterial("Tc98_mat");

  G4Material* targetMaterial = new G4Material("target_mat", boxMaterial->GetDensity(), 2);
  double tracerConcentration = 90.*perCent;
  targetMaterial->AddMaterial(tracerMaterial, tracerConcentration);
  targetMaterial->AddMaterial(boxMaterial, 100.*perCent - tracerConcentration);

  if( !airMaterial || !boxMaterial || !tracerMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("DetectorConstruction::DefineVolumes()", "MyCode0001", FatalException, msg);
  }
   
  //     
  // World
  //
  G4VSolid* worldS = new G4Box("WorldS",           // its name
                 worldSizeXYZ/2, worldSizeXYZ/2, worldSizeXYZ/2); // its size
                         
  G4LogicalVolume* worldLV = new G4LogicalVolume(
                 worldS,           // its solid
                 airMaterial,      // its material
                 "WorldLV");       // its name
                                   
  fWorldPV = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "WorldPV",        // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                               
  // Box
  //  
  G4VSolid* boxS = new G4Box("BoxS",     // its name
                 boxSizeXY/2, boxSizeXY/2, boxSizeZ/2); // its size
                         
  G4LogicalVolume* boxLV = new G4LogicalVolume(
                 boxS,     // its solid
                 boxMaterial,  // its material
                 "BoxLV");   // its name
                                   
  //                               
  // target
  //  
  G4VSolid* targetS = new G4Box("TargetS",     // its name
                 targetSizeXY/2, targetSizeXY/2, targetSizeZ/2); // its size
                         
  G4LogicalVolume* targetLV = new G4LogicalVolume(
                 targetS,     // its solid
                 targetMaterial,  // its material
                 "TargetLV");   // its name
  
  // placing the box in the world and then the target in the box
  fBoxPV = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 boxLV,            // its logical volume
                 "BoxPV",          // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  fTargetPV = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0*cm,0*cm,11.1*cm),  // at (0,0,0)
                 targetLV,         // its logical volume
                 "TargetPV",       // its name
                 boxLV,            // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
 
  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::Invisible);

  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(0.0,0.5,1.0)); // blue
  simpleBoxVisAtt->SetVisibility(true);
  boxLV->SetVisAttributes(simpleBoxVisAtt);
  
  G4VisAttributes* simpleTargetVisAtt= new G4VisAttributes(G4Colour(0.0,1.0,0.0)); // green
  simpleTargetVisAtt->SetVisibility(true);
  targetLV->SetVisAttributes(simpleTargetVisAtt);

  //
  // Always return the physical World
  //
  return fWorldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
