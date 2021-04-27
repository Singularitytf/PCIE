#include "debug.hh"
// using namespace std;
// #include "TGraph.h"
void PrtPoissonInfo();
void prtBeltPlottingData();

int main(){

     // double bkg;
  // std::cin >> bkg;
  auto PoisObj = new cvlPoisson(0.5, 3., .2);
  auto itv = new Interval(PoisObj, 0.9, 60, 0.005);
  // double mu = 0.2;
  // PoisObj->setBkgMean(2.6);
  // auto a = itv->findMuInterval(0);
  // std::cout << a.lower<< "\t" << a.upper << std::endl;
//   std::cout << "---------" << std::endl;
//   PoisObj->setBkgMean(2);
//   auto a = itv->roughMuScan(0);
//   std::cout << a.lower<< "\t" << a.upper << std::endl;
//   PoisObj->setBkgMean(3);
//   a = itv->roughMuScan(0);
//   std::cout << a.lower<< "\t" << a.upper << std::endl;
  
  
  FILE *fp = fopen("up_bkg.txt", "w");
  for(double bkg=2;bkg<10;bkg+=0.01){
    PoisObj->setBkgMean(bkg);
    std::cout << "Building With bkg: "<< PoisObj->getBkgMean() << std::endl;
    fprintf(fp, "%.4f\t%.4f\t%.4f\n", bkg, 
    itv->findMuInterval(9).upper,
    itv->findMuInterval(10).upper);
//     itv->roughMuScan(0).upper;
//     itv->roughMuScan(1).upper;
  }
  fclose(fp);
  
//     //PrtPoissonInfo();
//     auto ps = new cvlPoisson(.5, 3.0, 0.5);
//     ps->setBkgMean(10.);
//     auto bt = new Belt(ps, 0.9, 0, 60, 0.001);
//     roughMu a = bt->findMuInterval(5);
// //     cout << a.lower << "\t" << a.upper << endl;
//     for (int i=0;i<=20;++i){
//      a = bt->findMuInterval(i);
//      cout << a.lower << "\t" << a.upper << endl;
//     }
//     delete ps;
//     delete bt;
     // prtBeltPlottingData();
     //PrtPoissonInfo();
     //auto cvlPois = new cvlPoisson(0.5, 3., .4);
     //cout << cvlPois->findBestMu(0) << endl;
    return 0;
}

void prtBeltPlottingData(){
  //PrtPoissonInfo();
    auto ps = new cvlPoisson(.5, 2.0, 0.2);
    auto bt = new Belt(ps, 0.9, 0, 16, 0.005);
    FILE *fp = fopen("beltPlot.txt", "w");
    std::vector<n0Limit> belt = bt->constructBelt();
    const size_t sz = belt.size();
    for(int i=0; i<sz; i++){
     fprintf(fp, "%.3f\t%d\t%d\n", 0+i*0.005, belt[i].lower, belt[i].upper);
    }
     fclose(fp);
}

void PrtPoissonInfo(){
  using namespace std;
  double sr, bkg;
  sr = 0.5; bkg = 3.;
  printf("Constructing a Poisson Class with SR = %f, BKG = %f......\n",
         sr, bkg);
  auto *ps = new cvlPoisson(sr, bkg, 0.4);
  printf("The signal rate is: %f, the mean background is: %f\n",
         ps->getSRate(), ps->getBkgMean());
  sr = 1.5, bkg = 0.5;
  ps->setBkgMean(bkg);
  ps->setSRate(sr);
  cout << "Now setting sr = " << sr
       << ". bkg = " << bkg << "." << endl;
  printf("The signal rate is: %f, the mean background is: %f\n",
         ps->getSRate(), ps->getBkgMean());
  int n0 = 2;
  cout << "the prob for n0 = " << n0 << " is: "
       << ps->pmf(n0) << endl;
  cout << "And the best MU for this n0 is: "
       << ps->findBestMu(n0) << endl;

  double tmp = Poisson(4,0.5).pmf(2);
  cout << tmp << endl;
}
