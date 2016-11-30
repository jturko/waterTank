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
// $Id: B4RunAction.hh 74265 2013-10-02 14:41:20Z gcosmo $
// 
/// \file B4RunAction.hh
/// \brief Definition of the B4RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

/// Run action class
///
/// It accumulates statistic and computes dispersion of the energy deposit 
/// and track lengths of charged particles with use of analysis tools:
/// H1D histograms are created in BeginOfRunAction() for the following 
/// physics quantities:
/// - Edep in absorber
/// - Edep in gap
/// - Track length in absorber
/// - Track length in gap
/// The same values are also saved in the ntuple.
/// The histograms and ntuple are saved in the output file in a format
/// accoring to a selected technology in B4Analysis.hh.
///
/// In EndOfRunAction(), the accumulated statistic and computed 
/// dispersion is printed.
///

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

	 void Add(G4double energy, G4double time, G4int particle, G4double vertex_x, G4double vertex_y, G4double vertex_z, G4double pos_x, G4double pos_y, G4double pos_z, G4int trackId) {
		 fEnergy.push_back(energy);
		 fTime.push_back(time);
		 fParticle.push_back(particle);
		 fVertexX.push_back(vertex_x);
		 fVertexY.push_back(vertex_y);
		 fVertexZ.push_back(vertex_z);
		 fPosX.push_back(pos_x);
		 fPosY.push_back(pos_y);
		 fPosZ.push_back(pos_z);
		 fTrackId.push_back(trackId);
	 }

	 void Clear() {
		 fEnergy.clear();
		 fTime.clear();
		 fParticle.clear();
		 fVertexX.clear();
		 fVertexY.clear();
		 fVertexZ.clear();
		 fPosX.clear();
		 fPosY.clear();
		 fPosZ.clear();
		 fTrackId.clear();
	 }

  private:
	 std::vector<G4double> fEnergy;
	 std::vector<G4double> fTime;
	 std::vector<G4int>    fParticle;
	 std::vector<G4double> fVertexX;
	 std::vector<G4double> fVertexY;
	 std::vector<G4double> fVertexZ;
	 std::vector<G4double> fPosX;
	 std::vector<G4double> fPosY;
	 std::vector<G4double> fPosZ;
	 std::vector<G4int>    fTrackId;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

