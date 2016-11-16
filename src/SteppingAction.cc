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
// $Id: B4aSteppingAction.cc 68058 2013-03-13 14:47:43Z gcosmo $
// 
/// \file B4aSteppingAction.cc
/// \brief Implementation of the B4aSteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

//#include "G4AnalysisManager.hh"
#include "Analysis.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(
                      const DetectorConstruction* detectorConstruction,
                      EventAction* eventAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
// Collect energy and track length step by step

  // get volume of the current step
  G4VPhysicalVolume* firstVolume 
    = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  G4VPhysicalVolume* secondVolume 
    = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
  
  // energy deposit
  //G4double edep = step->GetTotalEnergyDeposit();
  G4double ekin = step->GetPreStepPoint()->GetKineticEnergy();
  G4double time = step->GetPreStepPoint()->GetGlobalTime();
  G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();
  G4int particle;
  if(particleName == "gamma")        particle = 1;
  else if(particleName == "proton")  particle = 2;
  else if(particleName == "e-")      particle = 3; 
  else if(particleName == "e+")      particle = 4; 
  else if(particleName == "neutron") particle = 5; 
  else if(particleName == "alpha")   particle = 6; 
  else if(particleName == "C12")     particle = 7;
  else if(particleName == "O16")     particle = 8; 
  else particle = 0; 

  if( firstVolume == fDetConstruction->GetWaterPV() && secondVolume == fDetConstruction->GetWorldPV() ) {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillNtupleDColumn(0, ekin/MeV);
    analysisManager->FillNtupleDColumn(1, time/ms);
    analysisManager->FillNtupleIColumn(2, particle);
    analysisManager->AddNtupleRow();
  }

  // step length
  //G4double stepLength = 0.;
  //if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
  //  stepLength = step->GetStepLength();
  //}
      
  //if ( volume == fDetConstruction->GetAbsorberPV() ) {
  //  fEventAction->AddAbs(edep,stepLength);
  //}
  
  //if ( volume == fDetConstruction->GetGapPV() ) {
  //  fEventAction->AddGap(edep,stepLength);
  //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
