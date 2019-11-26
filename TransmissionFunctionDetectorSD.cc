#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4ios.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "TrackInformation.hh"
#include "TransmissionFunctionDetectorSD.hh"

TransmissionFunctionDetectorSD::TransmissionFunctionDetectorSD(G4String name)
    : G4VSensitiveDetector(name) {}

TransmissionFunctionDetectorSD::~TransmissionFunctionDetectorSD() {}

void TransmissionFunctionDetectorSD::Initialize(G4HCofThisEvent *) {}

G4bool TransmissionFunctionDetectorSD::ProcessHits(G4Step *aStep,
                                                   G4TouchableHistory *ROhist) {
  G4Track *currentTrack = aStep->GetTrack();
  PrimaryGeneratorAction *genAction =
      (PrimaryGeneratorAction *)G4RunManager::GetRunManager()
          ->GetUserPrimaryGeneratorAction();
  DetectorConstruction *detConst =
      (DetectorConstruction *)G4RunManager::GetRunManager()
          ->GetUserDetectorConstruction();

  if (currentTrack->GetDefinition() == G4Gamma::Definition()) {
    G4double energy = aStep->GetPreStepPoint()->GetTotalEnergy();

    if (energy / keV == 0.) {
      currentTrack->SetTrackStatus(fStopAndKill);
      return false;
    } else if (energy > maxE) {
      currentTrack->SetTrackStatus(fStopAndKill);
      return false;
    } else if (energy < minE) {
      currentTrack->SetTrackStatus(fStopAndKill);
      return false;
    }

    int hitX;
    int hitY;

    hitX = (int)((aStep->GetPreStepPoint()->GetPosition().x() / pixelSizeX) +
                 (resX / 2));
    hitY = (int)((aStep->GetPreStepPoint()->GetPosition().y() / pixelSizeY) +
                 (resY / 2));
    PrimaryGeneratorAction *genAction =
        (PrimaryGeneratorAction *)G4RunManager::GetRunManager()
            ->GetUserPrimaryGeneratorAction();

    intBeta[(int)((energy - minE) / (maxE - minE) * (resE - 1))][hitX][hitY] +=
        genAction->GetIntBeta();
    intDelta[(int)((energy - minE) / (maxE - minE) * (resE - 1))][hitX][hitY] +=
        genAction->GetIntDelta();
  }
  currentTrack->SetTrackStatus(fStopAndKill);
  return true;
}

void TransmissionFunctionDetectorSD::InitImages() {
  intBeta.resize(resE);
  for (int i = 0; i < resE; i++) {
    intBeta[i].resize(resX);
    for (int j = 0; j < resX; j++) {
      intBeta[i][j].assign(resY, 0);
    }
  }

  intDelta.resize(resE);
  for (int i = 0; i < resE; i++) {
    intDelta[i].resize(resX);
    for (int j = 0; j < resX; j++) {
      intDelta[i][j].assign(resY, 0);
    }
  }
}

void TransmissionFunctionDetectorSD::ResetImages() {
  for (int ie = 0; ie < resE; ie++) {
    for (int ix = 0; ix < resX; ix++) {
      for (int iy = 0; iy < resY; iy++) {
        intBeta[ie][ix][iy] = 0;
      }
    }
  }

  for (int ie = 0; ie < resE; ie++) {
    for (int ix = 0; ix < resX; ix++) {
      for (int iy = 0; iy < resY; iy++) {
        intDelta[ie][ix][iy] = 0;
      }
    }
  }
}

void TransmissionFunctionDetectorSD::EndOfEvent(G4HCofThisEvent *HCE) {}

void TransmissionFunctionDetectorSD::clear() {}

void TransmissionFunctionDetectorSD::DrawAll() {}

void TransmissionFunctionDetectorSD::PrintAll() {}
