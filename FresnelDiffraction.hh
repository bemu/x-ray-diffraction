#ifndef FresnelDiffraction_h
#define FresnelDiffraction_h 1

#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <vector>

using namespace std;

class FresnelDiffraction {
public:
  FresnelDiffraction(FresnelDiffraction *FD, double waveNumber);
  ~FresnelDiffraction();

  void InitSphericalWave(double sourceIntensity, double initDistance,
                         double initPixelSizeX, double initPixelSizeY,
                         double waveNumber,
                         vector<vector<vector<double>>> *intBeta,
                         vector<vector<vector<double>>> *intDelta,
                         bool reference);
  void PropagateSphericalWave(double distance);
  void PropagatePlaneWave(double distance, double deltaX, double deltaY);
  void GratingAbsorption(double gratingPitch);

  double GetDistance() { return distanceTotal; };
  double GetPixelSizeX() { return pixelSizeX; };
  double GetPixelSizeY() { return pixelSizeY; };
  double GetResE() { return resE; };
  double GetResX() { return resX; };
  double GetResY() { return resY; };

  fftw_complex *psi;
  fftw_complex *psiK;

private:
  double distanceTotal;
  double pixelSizeX;
  double pixelSizeY;
  double k;
  int resE;
  int resX;
  int resY;
  bool realSpace;
};
#endif
