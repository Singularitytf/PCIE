#include "Interval.hh"

//....................................................................

Interval::Interval(cvlPoisson *fPoisObj, double fcl,
           double fmuULimit, double fmu_esp)
  : cfdent_level(fcl)
{
  PoisObj = fPoisObj;
  belt = new Belt(PoisObj, fcl, 0., fmuULimit, 30);
}

Interval::~Interval()
{
  delete PoisObj;
  delete belt;
}

roughMu Interval::roughMuScan(int n0){
  // printf("start Rough Scan!\n");
  double bestMu = (n0-PoisObj->getBkgMean())>0 ? 
  n0-PoisObj->getBkgMean() : 0.0;
  double upper_start_p, lower_start_p;
  double fesp = 0.1;
  int ntmp;
  roughMu scanRst;

  //Init start point.
  upper_start_p = bestMu;
  lower_start_p = (bestMu - 2*TMath::Sqrt(bestMu))<=0?
    0:(bestMu - TMath::Sqrt(bestMu));

  for (double i=lower_start_p; i>=0; i-=fesp){
    // printf("lower scan mu = %f", i);
    ntmp = belt->findHInterval(i).upper;
    // std::cout << i << "\t low" << ntmp << std::endl;
    if (ntmp == n0-1){
      scanRst.lower = i;
      break;
    }
  }

  for (double i=upper_start_p; i<=upper_start_p+10; i+=fesp){
    // printf("upper scan mu = %f", i);
    ntmp = belt->findHInterval(i).lower;
    //printf("%d\n", ntmp);
    // std::cout << i << "\t up" << ntmp << std::endl;
    if (ntmp == n0+1){
      scanRst.upper = i;
      break;
    }
  }

  //std::cout << lower_start_p << "\t" << upper_start_p << std::endl;
  return scanRst;

}

//....................................................................

roughMu Interval::findMuInterval(int fn0){

  int ntmp;
  // Init scan interval;
  roughMu roughmu = roughMuScan(fn0);
  // roughMu exacMu;
  double low_region_1, low_region_2, up_region_1, up_region_2;
  double mu_tmp;
  low_region_1 = roughmu.lower+0.5;
  low_region_2 = ((roughmu.lower-0.1)>0)?(roughmu.lower-0.1):0;
  up_region_1 = ((roughmu.upper-0.5)>0)?(roughmu.upper-0.5):0;
  up_region_2 = roughmu.upper+1.2;
  uint l_iter = (uint) (low_region_2-low_region_1)/esp;
  uint u_iter = (uint) (up_region_2-up_region_1)/esp;
  // scan upper limit

  for(int i=0; i<=u_iter; i++){
    mu_tmp = up_region_1+i*esp;
    ntmp = belt->findHInterval(mu_tmp).lower;
    if (ntmp == fn0) {
      roughmu.upper = mu_tmp;
      // printf("%.4f\t%d\n", exacMu.upper, ntmp);
    }
    
  }
  // printf("lo_lim: %.4f\n",low_scan_itv[0]);
  //scan lower limit
  for(int i=0; i<=l_iter;i++){
    mu_tmp = low_region_2-i*esp;
    ntmp = belt->findHInterval(mu_tmp).upper;
    // printf("%.4f\t%d\n", i, ntmp);
    if (ntmp == fn0){
       roughmu.lower = mu_tmp;
      //  printf("%.4f\t%d\n", exacMu.lower, ntmp);
    }
  }
  return roughmu;
}


//....................................................................

double Interval::dipModify(int n0){
  double mu2_now = findMuInterval(n0).upper;
  double bkg_now = PoisObj->getBkgMean();
  double mu2_scan;
  bool inDip = false;
  int riter_max = 1.5/0.05;
  for (int i=0; i<riter_max; i++){
    double bkg=bkg_now+i*0.05;//0.05是粗略扫描的精度；
    PoisObj->setBkgMean(bkg);
    mu2_scan = findMuInterval(n0).upper;
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
        mu2_scan = findMuInterval(n0).upper;
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
  // double bkg;
  // std::cin >> bkg;
  auto PoisObj = new cvlPoisson(0.5, 0.5, .2);
  auto bt = new Belt(PoisObj, 0.9, 0, 0.5);
  auto itv = new Interval(PoisObj, 3.0, 0.9, 60);
  double mu = 3;
  // PoisObj->setBkgMean(2);
  auto a = itv->roughMuScan(0);
  std::cout << a.lower<< "\t" << a.upper << std::endl;
  // std::cout << "---------" << std::endl;
  // PoisObj->setBkgMean(2);
  // a = itv->roughMuScan(0);
  // std::cout << a.lower<< "\t" << a.upper << std::endl;
  
  FILE *fp = fopen("up_bkg.txt", "w");
  for(double bkg=0.01;bkg<25;bkg+=0.005){
    PoisObj->setBkgMean(bkg);
    std::cout << PoisObj->getBkgMean() << std::endl;
    // fprintf(fp, "%.4f\t%.4f\t%.4f\n", bkg, 
    // itv->roughMuScan(0).upper,
    // itv->roughMuScan(1).upper);
    itv->roughMuScan(0).upper;
    itv->roughMuScan(1).upper;
  }
  fclose(fp);
  
  return 0;
}
