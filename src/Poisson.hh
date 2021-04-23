#pragma once
#include "TMath.h"
#include <iostream>

class Poisson{
private:
  double s_rate;
  double bkg_mean;

public:
  Poisson(double fsrate, double fbkg_mean);
  ~Poisson();
  double pmf(int n0);
  double cdf(int n0);
  void setSRate(double &fsrate);
  void setBkgMean(double &fbkg_mean);
  double getBkgMean();
  double getSRate();
  double findBestMu(int &n0);
  // void printInfo();
};