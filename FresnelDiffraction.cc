#include "FresnelDiffraction.hh"

using namespace std;

FresnelDiffraction::FresnelDiffraction(FresnelDiffraction *FD,
                                       double waveNumber) {
  k = waveNumber;
  distanceTotal = FD->distanceTotal;
  pixelSizeX = FD->pixelSizeX;
  pixelSizeY = FD->pixelSizeY;
  resE = FD->resE;
  resX = FD->resX;
  resY = FD->resY;

  if (FD->realSpace) {
    psi = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * resX * resY);
    memcpy(psi, FD->psi, sizeof(fftw_complex) * resX * resY);
    realSpace = true;
  } else {
    psiK = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * resX * resY);
    memcpy(psiK, FD->psiK, sizeof(fftw_complex) * resX * resY);
    realSpace = false;
  }
}

FresnelDiffraction::~FresnelDiffraction() {}

void FresnelDiffraction::InitSphericalWave(
    double sourceIntensity, double initDistance, double initPixelSizeX,
    double initPixelSizeY, double waveNumber,
    vector<vector<vector<double>>> *intBeta,
    vector<vector<vector<double>>> *intDelta, bool reference) {

  distanceTotal = initDistance;
  pixelSizeX = initPixelSizeX;
  pixelSizeY = initPixelSizeY;
  k = waveNumber;
  resE = (*intBeta).size();
  resX = (*intBeta)[0].size();
  resY = (*intBeta)[0][0].size();

  psi = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * resX * resY);
  psiK = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * resX * resY);
  fftw_plan plan =
      fftw_plan_dft_2d(resX, resY, psi, psiK, FFTW_FORWARD, FFTW_ESTIMATE);

  for (int ie = 0; ie < resE; ie++) {
    for (int ix = 0; ix < resX; ix++) {
      for (int iy = 0; iy < resY; iy++) {
        if (reference) {
          psi[iy + resY * ix][0] = sourceIntensity / initDistance;
          psi[iy + resY * ix][1] = 0.0;
        } else {
          // TODO: use complex<double> instead of __complex__
          __complex__ double c = sourceIntensity / initDistance *
                                 cexp(-I * k * (*intDelta)[ie][ix][iy]) *
                                 cexp(-k * (*intBeta)[ie][ix][iy]);
          psi[iy + resY * ix][0] = creal(c);
          psi[iy + resY * ix][1] = cimag(c);
        }
      }
    }
  }

  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_free(psi);
  realSpace = false;
}

void FresnelDiffraction::PropagatePlaneWave(double distance, double deltaX,
                                            double deltaY) {
  double lambda = 2 * M_PI / k;
  // TODO: make configurable instead of hardcoding propagation coefficients
  double propagCoeffDelta = 0.5;
  double propagCoeffSigma = 0.3;
  double lX = propagCoeffDelta * distance / pixelSizeX / 2;
  double lY = propagCoeffDelta * distance / pixelSizeY / 2;
  double sigmaX =
      propagCoeffSigma * propagCoeffDelta * distance / pixelSizeX / 2;
  double sigmaY =
      propagCoeffSigma * propagCoeffDelta * distance / pixelSizeY / 2;

  fftw_complex *propagator =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * resX * resY);
  fftw_complex *propagatorK =
      (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * resX * resY);
  fftw_plan plan = fftw_plan_dft_2d(resX, resY, propagator, propagatorK,
                                    FFTW_FORWARD, FFTW_ESTIMATE);

  for (int ie = 0; ie < resE; ie++) {
    for (int ix = 0; ix < resX; ix++) {
      for (int iy = 0; iy < resY; iy++) {
        // take into account magnification
        double x, y;
        if (ix < resX / 2)
          x = deltaX * (double)ix;
        else
          x = deltaX * (double)(ix - resX);
        if (iy < resY / 2)
          y = deltaY * (double)iy;
        else
          y = deltaY * (double)(iy - resY);

        double ratioX = 0;
        double ratioY = 0;
        if (fabs(x) > lambda * lX)
          ratioX = (fabs(x) - lambda * lX) / sigmaX / lambda;
        if (fabs(y) > lambda * lY)
          ratioY = (fabs(y) - lambda * lY) / sigmaY / lambda;

        __complex__ double c =
            cexp(I * k / 2 / distance * (pow(x, 2) + pow(y, 2))) *
            cexp(-ratioX * ratioX / 2) * cexp(-ratioY * ratioY / 2);
        propagator[iy + resY * ix][0] = creal(c);
        propagator[iy + resY * ix][1] = cimag(c);
      }
    }
  }

  fftw_execute(plan);

  for (int ie = 0; ie < resE; ie++) {
    for (int ikx = 0; ikx < resX; ikx++) {
      for (int iky = 0; iky < resY; iky++) {
        __complex__ double c1 = propagatorK[iky + resY * ikx][0] +
                                I * propagatorK[iky + resY * ikx][1];
        __complex__ double c2 =
            psiK[iky + resY * ikx][0] + I * psiK[iky + resY * ikx][1];
        __complex__ double c = c1 * c2 * deltaX * deltaY / resX / resY;
        psiK[iky + resY * ikx][0] = creal(c);
        psiK[iky + resY * ikx][1] = cimag(c);
      }
    }
  }

  fftw_destroy_plan(plan);
  fftw_free(propagator);
  fftw_free(propagatorK);

  psi = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * resX * resY);
  plan = fftw_plan_dft_2d(resX, resY, psiK, psi, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_free(psiK);
  realSpace = true;
}

void FresnelDiffraction::PropagateSphericalWave(double distance) {
  // take into account magnification
  double M = (distanceTotal + distance) / distanceTotal;
  PropagatePlaneWave(distance / M, pixelSizeX, pixelSizeY);
  pixelSizeX *= M;
  pixelSizeY *= M;
  distanceTotal += distance;
  // TODO: correct the total intensity taking into account energy conservation
}

void FresnelDiffraction::GratingAbsorption(double gratingPitch) {
  psiK = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * resX * resY);
  fftw_plan plan =
      fftw_plan_dft_2d(resX, resY, psi, psiK, FFTW_FORWARD, FFTW_ESTIMATE);

  for (int ie = 0; ie < resE; ie++) {
    for (int ix = 0; ix < resX; ix++) {
      for (int iy = 0; iy < resY; iy++) {
        __complex__ double c =
            (psi[iy + resY * ix][0] + I * psi[iy + resY * ix][1]) *
            ((int)((double)ix * pixelSizeX / gratingPitch) % 2);
        psi[iy + resY * ix][0] = creal(c);
        psi[iy + resY * ix][1] = creal(c);
      }
    }
  }
  // TODO: implement exp(-x^2) decay at the edges to avoid oscillations

  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_free(psi);
  realSpace = false;
}
