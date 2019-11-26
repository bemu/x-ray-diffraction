#ifndef TransmissionFunctionDetectorSD_h
#define TransmissionFunctionDetectorSD_h 1

#include "G4VPrimitiveScorer.hh"
#include "G4VSensitiveDetector.hh"

class DetectorConstruction;
class G4HCofThisEvent;
class G4Step;

using namespace std;

class TransmissionFunctionDetectorSD : public G4VSensitiveDetector {
public:
  TransmissionFunctionDetectorSD(G4String);
  ~TransmissionFunctionDetectorSD();

  void Initialize(G4HCofThisEvent *);
  G4bool ProcessHits(G4Step *astep, G4TouchableHistory *ROHist);
  void EndOfEvent(G4HCofThisEvent *);
  void clear();
  void DrawAll();
  void PrintAll();

  G4double GetSizeX() { return sizeX; };
  G4double GetSizeY() { return sizeY; };
  G4double GetSizeZ() { return sizeZ; };

  G4double GetPosZ() { return posZ; };

  G4double GetPixelSizeX() { return pixelSizeX; };
  G4double GetPixelSizeY() { return pixelSizeY; };

  G4int GetResX() { return resX; };
  G4int GetResY() { return resY; };
  G4int GetResE() { return resE; };

  vector GetIntBeta() { return intBeta; };
  vector GetIntDelta() { return intDelta; };

  void SetSizeX(G4double s) { sizeX = s; };
  void SetSizeY(G4double s) { sizeY = s; };
  void SetSizeZ(G4double s) { sizeZ = s; };

  void SetPosZ(G4double p) { posZ = p; };

  void SetResX(G4int r) { resX = r; };
  void SetResY(G4int r) { resY = r; };

  void SetResE(G4int e) { resE = e; };
  void SetMinE(G4double e) { minE = e; };
  void SetMaxE(G4double e) { maxE = e; };

  G4double GetMinE(void) { return minE; };
  G4double GetMaxE(void) { return maxE; };

  void SetPixelSizeX(G4double s) { pixelSizeX = s; };
  void SetPixelSizeY(G4double s) { pixelSizeY = s; };

  void InitImages();
  void ResetImages();

private:
  DetectorConstruction *detectorConstruction;

  vector<vector<vector<G4double>>> intBeta;
  vector<vector<vector<G4double>>> intDelta;

  G4int resX;
  G4int resY;
  G4int resE;

  G4double minE;
  G4double maxE;

  G4double sizeX;
  G4double sizeY;
  G4double sizeZ;

  G4double posZ;

  G4double pixelSizeX;
  G4double pixelSizeY;
};
#endif
