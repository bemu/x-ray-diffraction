#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "G4VUserPhysicsList.hh"

class G4LowEnergyIonisation;
class G4LowEnergyPhotoElectric;
class G4LowEnergyBremsstrahlung;

enum PhysicsModel { Livermore, Penelope };
enum ComptonModel {
  ComptonLivermore,
  ComptonPenelope,
  KleinNishina,
  ComptonLowEP,
  LivermoreModified
};

class ModularPhysicsList : public G4VModularPhysicsList {
public:
  ModularPhysicsList();
  ~ModularPhysicsList();

  void Configure(string physicsModelName, string comptonModelName);

  void SetElectronCut(G4double);
  void SetPositronCut(G4double);
  void SetGammaCut(G4double);

protected:
  void ConstructParticle();
  void ConstructProcess();

  void SetCuts();

  void ConstructBosons();
  void ConstructLeptons();

  void ConstructGeneral();
  void ConstructEM();

private:
  G4double cutGamma;
  G4double cutElectron;
  G4double cutPositron;

  PhysicsModel physicsModel;
  ComptonModel comptonModel;
};
#endif
