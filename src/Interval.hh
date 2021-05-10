#include "Belt.hh"


class Interval{
private:
  double cfdent_level;
  double esp;
  Belt *belt = NULL;//Private
  cvlPoisson *PoisObj = NULL;//Private
  // Poisson *PoisObj = NULL;
  // Belt *belt = NULL;
public:
  Interval(cvlPoisson *fPoisObj, double fcl,
           double fmuULimit = 60, double fmu_esp = 0.005);
  ~Interval();
  double dipModify(int n0);
  roughMu getInterval();
  roughMu roughMuScan(int n0);
  roughMu findMuInterval(int fn0);
  double findUpperLimit(int fn0);
  double findLowerLimit(int fn0);
  void setCL(double fcl);
  double getCL();
  // void setBkg(double fbkg);
};
