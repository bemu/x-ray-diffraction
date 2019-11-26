#include "G4EmProcessOptions.hh"
#include "G4LossTableManager.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4RunManager.hh"
#include "G4StepLimiter.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "IntegrateRefractiveIndexProcess.hh"

#include "G4EmPenelopePhysics.hh"
#include "G4PenelopeAnnihilationModel.hh"
#include "G4PenelopeBremsstrahlungModel.hh"
#include "G4PenelopeComptonModel.hh"
#include "G4PenelopeGammaConversionModel.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4PenelopePhotoElectricModel.hh"
#include "G4PenelopeRayleighModel.hh"

#include "G4EmLivermorePhysics.hh"
#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermoreComptonModifiedModel.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreRayleighModel.hh"

#include "G4KleinNishinaCompton.hh"
#include "G4KleinNishinaModel.hh"
#include "G4LowEPComptonModel.hh"

#include "G4SystemOfUnits.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "RunAction.hh"

using namespace std;

ModularPhysicsList::ModularPhysicsList() : G4VModularPhysicsList() {
  cutDefault = 10 * micrometer;
  cutGamma = 10 * micrometer;
  cutElectron = cutDefault;
  cutPositron = cutDefault;

  physicsModel = Livermore;
  comptonModel = ComptonLivermore;
}

ModularPhysicsList::~ModularPhysicsList() {}

void ModularPhysicsList::ConstructParticle() {
  ConstructBosons();
  ConstructLeptons();
}

void ModularPhysicsList::Configure(string physicsModelName,
                                   string comptonModelName) {
  if (physicsModelName == "Livermore") {
    physicsModel = Livermore;
  } else if (physicsModelName == "Penelope") {
    physicsModel = Penelope;
  }

  if (comptonModelName == "Livermore") {
    comptonModel = ComptonLivermore;
  } else if (comptonModelName == "Penelope") {
    comptonModel = ComptonPenelope;
  } else if (comptonModelName == "KleinNishina") {
    comptonModel = KleinNishina;
  } else if (comptonModelName == "LowEP") {
    comptonModel = ComptonLowEP;
  } else if (comptonModelName == "LivermoreModified") {
    comptonModel = LivermoreModified;
  }
}

void ModularPhysicsList::ConstructBosons() { G4Gamma::GammaDefinition(); }

void ModularPhysicsList::ConstructLeptons() {
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}

void ModularPhysicsList::ConstructProcess() {
  AddTransportation();
  ConstructEM();

  G4VAtomDeexcitation *de = new G4UAtomicDeexcitation();
  de->SetDeexcitationActiveRegion("World", true, false, false);
  de->SetFluo(true);
  de->SetAuger(false);
  de->SetPIXE(false);
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);
}

void ModularPhysicsList::ConstructEM() {
  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition *particle = theParticleIterator->value();
    G4ProcessManager *pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
      IntegrateRefractiveIndexProcess *theContinuousProcess =
          new IntegrateRefractiveIndexProcess("IntegrateRefractiveIndex");
      if (physicsModel == Penelope) {
        theContinuousProcess->SetModels(new G4PenelopePhotoElectricModel(),
                                        new G4PenelopeRayleighModel(),
                                        new G4PenelopeComptonModel());
      } else {
        theContinuousProcess->SetModels(new G4LivermorePhotoElectricModel(),
                                        new G4LivermoreRayleighModel(),
                                        new G4LivermoreComptonModel());
      }
      pmanager->AddContinuousProcess(theContinuousProcess);
    } else if (particleName == "e-") {
      G4eIonisation *theIonisation = new G4eIonisation();
      if (physicsModel == Penelope) {
        theIonisation->SetEmModel(new G4PenelopeIonisationModel());
        pmanager->AddProcess(theIonisation, -1, 2, 2);
      } else {
        G4LivermoreIonisationModel *theLivermoreIonisationModel =
            new G4LivermoreIonisationModel();
        theIonisation->SetEmModel(theLivermoreIonisationModel);
        pmanager->AddProcess(theIonisation, -1, 2, 2);
      }
      G4eMultipleScattering *msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      pmanager->AddProcess(msc, -1, 1, 1);
      G4eBremsstrahlung *theBremsstrahlung = new G4eBremsstrahlung();
      if (physicsModel == Penelope) {
        theBremsstrahlung->SetEmModel(new G4PenelopeBremsstrahlungModel());
        pmanager->AddProcess(theBremsstrahlung, -1, -3, 3);
      } else {
        theBremsstrahlung->SetEmModel(new G4LivermoreBremsstrahlungModel());
        pmanager->AddProcess(theBremsstrahlung, -1, -3, 3);
      }
    }
  }
}

void ModularPhysicsList::ConstructGeneral() {}

void ModularPhysicsList::SetCuts() {
  if (verboseLevel > 0) {
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(cutDefault, "Length") << G4endl;
  }

  SetCutValue(cutGamma, "gamma");
  SetCutValue(cutElectron, "e-");
  SetCutValue(cutPositron, "e+");

  G4double limitLow = 250. * eV;
  G4double limitHigh = 100. * GeV;
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(limitLow,
                                                                  limitHigh);

  if (verboseLevel > 0)
    DumpCutValuesTable();
}
