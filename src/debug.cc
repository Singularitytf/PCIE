#include "debug.hh"
using namespace std;
// #include "TGraph.h"
void PrtPoissonInfo();
void prtBeltPlottingData();

int main(){
//     //PrtPoissonInfo();
//     auto ps = new Poisson(.5, 3.0);
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
     PrtPoissonInfo();
    return 0;
}

void prtBeltPlottingData(){
  //PrtPoissonInfo();
    auto ps = new Poisson(.5, 21.0);
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
  Poisson *ps = new Poisson(sr, bkg);
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
