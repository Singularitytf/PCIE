#include "Belt.hh"


class Interval{
private:
  double bkg_mean;
  double cfdent_level;
  // Poisson *PoisObj = NULL;
  //Belt *belt = NULL;
public:
  Interval(double fbkg, double fcl,
           double fmuULimit = 60, double fmu_esp = 0.005);
  ~Interval();
  double dipModify(int n0);
  roughMu getInterval();
  Belt *belt = NULL;//Private
  Poisson *PoisObj = NULL;//Private
};
