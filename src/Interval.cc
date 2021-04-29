#include "Interval.hh"

//....................................................................

Interval::Interval(cvlPoisson *fPoisObj, double fcl,
           double fmuULimit, double fmu_esp)
  : cfdent_level(fcl),
    esp(fmu_esp)
{
  PoisObj = fPoisObj;
  belt = new Belt(PoisObj, fcl, 0., fmuULimit, 30);
}

Interval::~Interval()
{
  // delete PoisObj;
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
    // std::cout << i << "\tlow\t" << ntmp << std::endl;
    if (ntmp == n0-1){
      scanRst.lower = i;
      break;
    }
  }

  for (double i=upper_start_p; i<=upper_start_p+50; i+=fesp){
    // printf("upper scan mu = %f", i);
    ntmp = belt->findHInterval(i).lower;
    double bkg = PoisObj->getBkgMean();
    //printf("%d\n", ntmp);
    // std::cout << i << "\t" << bkg <<"\tup\t" << ntmp << std::endl;
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
  roughMu exacMu;
  double low_region_1, low_region_2, up_region_1, up_region_2;
  double mu_tmp;
  low_region_2 = roughmu.lower+0.5;
  low_region_1 = ((roughmu.lower-0.1)>0)?(roughmu.lower-0.1):0;
  up_region_1 = ((roughmu.upper-0.5)>0)?(roughmu.upper-0.5):0;
  up_region_2 = roughmu.upper+1.2;
  uint l_iter = (uint) ((low_region_2-low_region_1)/esp);
  uint u_iter = (uint) ((up_region_2-up_region_1)/esp);
  // scan upper limit

  // printf("low scan region: %f\t%f\t%d\n", low_region_1, low_region_2, (uint)((low_region_2-low_region_1)/esp));
  // printf("up  scan region: %f\t%f\t\n", up_region_1, up_region_2);
  // printf("uiter: %f\n", esp);
  for(int i=0; i<=u_iter; i++){
    mu_tmp = up_region_1+i*esp;
    ntmp = belt->findHInterval(mu_tmp).lower;
    // printf("fine upper scan, mu:%f, n0:%d\n", mu_tmp, ntmp);
    if (ntmp == fn0) {
      exacMu.upper = mu_tmp;
      // printf("update upper limit: %.4f\tn0=%d\n", exacMu.upper, ntmp);
    }
    
  }
  // printf("lo_lim: %.4f\n",low_scan_itv[0]);
  //scan lower limit
  for(int i=0; i<=l_iter;i++){
    mu_tmp = low_region_2-i*esp;
    ntmp = belt->findHInterval(mu_tmp).upper;
    // printf("%.4f\t%d\n", i, ntmp);
    // printf("fine lower scan, mu:%f, n0:%d\n", mu_tmp, ntmp);
    if (ntmp == fn0){
       exacMu.lower = mu_tmp;
      //  printf("%.4f\t%d\n", exacMu.lower, ntmp);
    }
  }
  return exacMu;
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

