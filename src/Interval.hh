#include "Belt.hh"


class Interval{
private:
  double cfdent_level;
  double esp;
  // Poisson *PoisObj = NULL;
  // Belt *belt = NULL;
public:
  Interval(cvlPoisson *fPoisObj, double fcl,
           double fmuULimit = 60, double fmu_esp = 0.005);
  ~Interval();
  double dipModify(int n0);
  roughMu getInterval();
  Belt *belt = NULL;//Private
  cvlPoisson *PoisObj = NULL;//Private
  roughMu roughMuScan(int n0);
  roughMu findMuInterval(int fn0);
  // void setBkg(double fbkg);
};
