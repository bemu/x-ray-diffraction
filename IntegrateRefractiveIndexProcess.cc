#include "IntegrateRefractiveIndexProcess.hh"
#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreRayleighModel.hh"
#include "G4PEEffectFluoModel.hh"
#include "G4PenelopeComptonModel.hh"
#include "G4PenelopePhotoElectricModel.hh"
#include "G4PenelopeRayleighModel.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4PhysicalConstants.hh"
#include "G4RayleighScattering.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "PrimaryGeneratorAction.hh"

IntegrateRefractiveIndexProcess::IntegrateRefractiveIndexProcess(
    const G4String &processName, G4ProcessType type)
    : G4VContinuousProcess(processName, type) {}

IntegrateRefractiveIndexProcess::~IntegrateRefractiveIndexProcess() {}

void IntegrateRefractiveIndexProcess::SetModels(G4VEmModel *mPhotoElectric,
                                                G4VEmModel *mRayleigh,
                                                G4VEmModel *mCompton) {
  G4DataVector theCuts;

  photoElectricModel = mPhotoElectric;
  photoElectricModel->SetLowEnergyLimit(0 * keV);
  photoElectricModel->SetHighEnergyLimit(1000 * keV);
  photoElectricModel->Initialise(
      (const G4ParticleDefinition *)(G4Gamma::Gamma()), theCuts);

  rayleighModel = mRayleigh;
  rayleighModel->SetLowEnergyLimit(0 * keV);
  rayleighModel->SetHighEnergyLimit(1000 * keV);
  rayleighModel->Initialise((const G4ParticleDefinition *)(G4Gamma::Gamma()),
                            theCuts);

  comptonModel = mCompton;
  comptonModel->SetLowEnergyLimit(0 * keV);
  comptonModel->SetHighEnergyLimit(1000 * keV);
  comptonModel->Initialise((const G4ParticleDefinition *)(G4Gamma::Gamma()),
                           theCuts);
}

G4VParticleChange *
IntegrateRefractiveIndexProcess::AlongStepDoIt(const G4Track &aTrack,
                                               const G4Step &aStep) {

  PrimaryGeneratorAction *genAction =
      (PrimaryGeneratorAction *)G4RunManager::GetRunManager()
          ->GetUserPrimaryGeneratorAction();

  aParticleChange.Initialize(aTrack);

  G4double energy = genAction->particleGun->GetParticleEnergy();

  G4double stepLength = aStep.GetStepLength();
  G4double beta = 1.02135e-10; // TODO: remove harccoded values fof PMMA
  G4double delta = 2.96115e-7; // TODO: remove harccoded values fof PMMA

  const G4Material *m = aTrack.GetMaterial();
  /* TODO: calculate beta and delta from atomic form factors
  double mu = 0.0;
  int i = 0;

  const G4ElementVector *matVector = m->GetElementVector();

  for (G4ElementVector::const_iterator e = matVector->begin();
       e != matVector->end(); e++) {
    mu += Avogadro * m->GetDensity() * m->GetFractionVector()[i] /
          (*e)->GetA() *
          photoElectricModel->ComputeCrossSectionPerAtom(
              (G4Gamma::Gamma()), energy, (*e)->GetZ(), 0.0, 0.0, 0.0);
    i++;
  }

  i = 0;
  for (G4ElementVector::const_iterator e = matVector->begin();
       e != matVector->end(); e++) {
    mu +=
        Avogadro * m->GetDensity() * m->GetFractionVector()[i] / (*e)->GetA() *
        rayleighModel->ComputeCrossSectionPerAtom((G4Gamma::Gamma()), energy,
                                                  (*e)->GetZ(), 0.0, 0.0, 0.0);
    i++;
  }

  i = 0;
  for (G4ElementVector::const_iterator e = matVector->begin();
       e != matVector->end(); e++) {
    mu += Avogadro * m->GetDensity() * m->GetFractionVector()[i] /
          (*e)->GetA() *
          comptonModel->ComputeCrossSectionPerAtom((G4Gamma::Gamma()), energy,
                                                   (*e)->GetZ(), 0.0, 0.0, 0.0);
    i++;
  }
  */

  genAction->IncIntBeta(stepLength * beta);
  genAction->IncIntDelta(stepLength * delta);

  aParticleChange.ClearDebugFlag();
  aParticleChange.ProposeLocalEnergyDeposit(energy * (weightStart - weightEnd));
  aParticleChange.SetNumberOfSecondaries(0);

  return &aParticleChange;
}

G4double IntegrateRefractiveIndexProcess::GetContinuousStepLimit(
    const G4Track &aTrack, G4double, G4double currentMinimumStep, G4double &) {
  return currentMinimumStep;
}
