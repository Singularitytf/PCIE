#pragma once
#include "cvlPoisson.hh"
#include "Poisson.hh"
#include <algorithm>
#include <vector>

struct constructData{
  int n0;
  double prob;
  double prob_best;
  double rate;
};

struct n0Limit{
  int lower;
  int upper;
};

struct roughMu{
  double lower = 0.;
  double upper;
};

class Belt{
private:
  cvlPoisson *PoisObj = NULL;// = NULL;
  double cfdent_level;
  double mu_scan_min;
  double mu_scan_max;
  double mu_scan_esp;


  // Define Ratio Sort Rule.
  static bool rSortRule(constructData &a, constructData &b);
  // Define n0 Sort Rule.
  static bool n0SortRule(constructData &a, constructData &b);

  constructData fillStructure(int &fn0);


public:
  Belt(cvlPoisson *fPoisObj, double fcfdent_level,
       double fmu_scan_min, double fmu_scan_max,
       double fmu_scan_esp=0.05);
  ~Belt();
  n0Limit findHInterval(double &mu);
  void printConstruction(double &mu);
  std::vector<n0Limit> constructBelt();
  void outPutBelt(const std::string filename = "beltPlot.txt");
};
