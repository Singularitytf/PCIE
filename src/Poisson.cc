#include "Poisson.hh"

Poisson::Poisson(double fsrate, double fbkg_mean)
  : s_rate(fsrate),
    bkg_mean(fbkg_mean)
{
}

Poisson::~Poisson()
{
}

double Poisson::pmf(int n0){
  return TMath::Poisson((double)n0, s_rate+bkg_mean);
}

double Poisson::cdf(int n0){
  double temp = 0;
  for (int i=0; i<=n0; i++){
    temp += this->pmf(n0);
  }
  return temp;
}

void Poisson::setSRate(double &fsrate) {s_rate = fsrate;}
void Poisson::setBkgMean(double &fbkg_mean) {bkg_mean = fbkg_mean;}
double Poisson::getBkgMean() {return this->bkg_mean;}
double Poisson::getSRate() {return this->s_rate;}
double Poisson::findBestMu(int &n0) {return std::max(0., n0-bkg_mean);}
/*
void Poisson::printInfo(){
  using namespace std;

  printf("Signal rate is %f, bkg is %f\n", this->getSRate(),this->getBkgMean());
}
*/
