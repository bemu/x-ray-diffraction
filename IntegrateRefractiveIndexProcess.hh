#ifndef IntegrateRefractiveIndexProcess_hh
#define IntegrateRefractiveIndexProcess_hh

#include "G4ComptonScattering.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4VContinuousProcess.hh"
#include "G4VEmModel.hh"

class IntegrateRefractiveIndexProcess : public G4VContinuousProcess {
public:
  IntegrateRefractiveIndexProcess(const G4String &processName,
                                  G4ProcessType type = fElectromagnetic);
  ~IntegrateRefractiveIndexProcess();

  G4VParticleChange *AlongStepDoIt(const G4Track &track, const G4Step &aStep);

  G4double GetContinuousStepLimit(const G4Track &aTrack, G4double,
                                  G4double currentMinimumStep, G4double &);

  void SetModels(G4VEmModel *, G4VEmModel *, G4VEmModel *);

private:
  G4VEmModel *photoElectricModel;
  G4VEmModel *rayleighModel;
  G4VEmModel *comptonModel;
};
#endif
