#include "TMath.h"
#include "TF1.h"
#include <iostream>
#include <gsl/gsl_integration.h>

const double esp=1e-7;




class cvlPoisson{
private:
  double s_rate; // signal rate.
  double bkg_mean; //  mean_bkg rate.
  double bkg_sigma; // Gaussian uncertainty of bkg.
  static double fcvlPoisson(double bkg, void *para);
  // gsl_integration_workspace * w;
  //double muRoughScan();
  double result, error;

public:
  cvlPoisson(double fsrate, double fbkg_mean, double fbkg_sigma);
  ~cvlPoisson();


  void setSRate(double fsrate);
  void setBkgMean(double fbkg_mean);
  void setBkgSigma(double fbkg_sigma);
  double getSRate();
  double getBkgMean();
  double getBkgSigma();

  double pmf(int n0);
  double cdf(int n0);
  double findBestMu(int n0);
};
