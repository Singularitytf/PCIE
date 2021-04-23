#include "Interval.hh"

//....................................................................

Interval::Interval(double fbkg, double fcl,
                   double fmuULimit, double fmu_scan)
  : bkg_mean(fbkg),
    cfdent_level(fcl)
{
  PoisObj = new Poisson(0.5, bkg_mean);
  belt = new Belt(PoisObj, fcl, 0., fmuULimit, fmu_scan);
}

Interval::~Interval()
{
  delete PoisObj;
  delete belt;
}

//....................................................................

double Interval::dipModify(int n0){
  double mu2_now = belt->findMuInterval(n0).upper;
  double bkg_now = PoisObj->getBkgMean();
  double mu2_scan;
  bool inDip = false;
  int riter_max = 1.5/0.05;
  for (int i=0; i<riter_max; i++){
    double bkg=bkg_now+i*0.05;//0.05是粗略扫描的精度；
    PoisObj->setBkgMean(bkg);
    mu2_scan = belt->findMuInterval(n0).upper;
    // printf("%f\t%.4f\t%.4f\n", bkg, mu2_scan, mu2_now);
    if (mu2_scan > mu2_now){
      inDip = true;// it shows that mu2 in a dip;
      break;
    }
  }

  if (inDip) {
      double roughBkg = PoisObj->getBkgMean();
      double mu2_tmp = mu2_scan;
      for (double bkg=roughBkg-0.05; bkg<=roughBkg+0.05; bkg+=0.001){
        PoisObj->setBkgMean(bkg);
        mu2_scan = belt->findMuInterval(n0).upper;
        if (mu2_scan > mu2_tmp) mu2_tmp = mu2_scan;
      }
      // printf("i\n");
      return mu2_tmp;
  }
  else
    return mu2_now;

}

//....................................................................

int main(){
  auto itv = new Interval(12.0, 0.6827, 60);
  // auto a = itv->belt->findMuInterval(0);
  // FILE *fp = fopen("up_bkg.txt", "w");
  // for(double bkg=0;bkg<25;bkg+=0.005){
  //   itv->PoisObj->setBkgMean(bkg);
  //   fprintf(fp, "%.4f\t%.4f\t%.4f\n", bkg, itv->belt->findMuInterval(10).upper, itv->belt->findMuInterval(11).upper);
  // }
  // fclose(fp);
  std::cout << itv->dipModify(7) << std::endl;
  return 0;
}
