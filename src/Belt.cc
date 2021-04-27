#include "Belt.hh"

//....................................................................

Belt::Belt(cvlPoisson *fPoisObj, double fcfdent_level,
           double fmu_scan_min, double fmu_scan_max,
           double fmu_scan_esp)
  : mu_scan_min(fmu_scan_min),
    mu_scan_max(fmu_scan_max),
    mu_scan_esp(fmu_scan_esp)
{
  //strData.clear();
  PoisObj = fPoisObj;
  if (fcfdent_level >=0 && fcfdent_level <=1)
    cfdent_level = fcfdent_level;
}

Belt::~Belt()
{
}

//....................................................................

bool Belt::rSortRule(constructData &a, constructData &b){
  // Sort ratio.
  return (a.rate >= b.rate);
}

//....................................................................

bool Belt::n0SortRule(constructData &a, constructData &b){
  return (a.n0 > b.n0);
}

//....................................................................

constructData Belt::fillStructure(int &fn0){
  constructData fstrData;
  fstrData.n0 = fn0;
  fstrData.prob = PoisObj->pmf(fn0);
  fstrData.prob_best = cvlPoisson(PoisObj->findBestMu(fn0),
                               PoisObj->getBkgMean(),
                               PoisObj->getBkgSigma()).pmf(fn0);
  // std::cout << PoisObj->findBestMu(fn0) << std::endl;
  fstrData.rate = fstrData.prob/fstrData.prob_best;
  return fstrData;
}

//....................................................................

n0Limit Belt::findHInterval(double &mu){
  std::vector<constructData> strData;
  int n0_min, n0_max;
  int max_iter_idx;
  double total_porb = 0;
  n0Limit nl;

  // Scan n0 located at the center of the belt;
  PoisObj->setSRate(mu);
  int temp = ((int)(PoisObj->getBkgMean()+PoisObj->getSRate()) - 10);
  temp > 0 ? n0_min = temp : n0_min = 0;
  n0_max = (int)(PoisObj->getBkgMean()+PoisObj->getSRate()) + 10;

  // Scan for n0.
  for (int i=n0_min; i<=n0_max; i++){
    strData.push_back(fillStructure(i));
  }

  std::sort(strData.begin(), strData.end(), rSortRule);

  for (size_t i=0;i<strData.size(); i++){
    total_porb += strData[i].prob;
    if (total_porb >= this->cfdent_level){
      max_iter_idx = i;
      break;
    }
  }
  auto result = std::minmax_element(strData.begin(),
                                    strData.begin()+max_iter_idx+1,
                                    n0SortRule);
  nl.lower = (*result.second).n0;
  nl.upper = (*result.first).n0;
  return nl;

}

//....................................................................

void Belt::printConstruction(double &mu){
  std::vector<constructData> strData;
  int n0_min, n0_max;
  double total_porb = 0;
  int hInterval[2];
  // Scan n0 located at the center of the belt;
  PoisObj->setSRate(mu);
  int temp = ((int)(PoisObj->getBkgMean()+PoisObj->getSRate()) - 10);
  temp > 0 ? n0_min = temp : n0_min = 0;
  n0_max = (int)(PoisObj->getBkgMean()+PoisObj->getSRate()) + 10;
  // Scan for n0.
  for (int i=n0_min; i<=n0_max; i++){
    strData.push_back(fillStructure(i));
  }
  std::sort(strData.begin(), strData.end(), rSortRule);
  printf("n0\tPorb\tPorb Best\tRatio\n");
  for(size_t i=0; i<strData.size(); i++){
    printf("%d\t%.3f\t%.3f\t%.3f\n",
           strData[i].n0,
           strData[i].prob,
           strData[i].prob_best,
           strData[i].rate);
  }
}

//....................................................................

std::vector<n0Limit> Belt::constructBelt(){
  std::vector<n0Limit> belt;
  for(double i=this->mu_scan_min;
      i<=this->mu_scan_max;
      i+=this->mu_scan_esp){
    belt.push_back(findHInterval(i));
    // std::cout << i <<std::endl;
  }
  return belt;
}

//....................................................................

void Belt::outPutBelt(const std::string filename){
  FILE *fp = fopen(filename.data(), "w");
  std::vector<n0Limit> belt = this->constructBelt();

  for(int i=0; i<belt.size(); i++){
    fprintf(fp, "%.3f\t%d\t%d\n",
            i*mu_scan_esp,
            belt[i].lower,
            belt[i].upper);
    }
  fclose(fp);
}


//....................................................................


//....................................................................
/*
int main(){
  auto ps = Poisson(0.5, 3.);
  auto bt= Belt(&ps, 0.9, 0, 30);
  int *a = bt.findHInterval(6.);
  std::cout << a[0] << "\t" << a[1] << std::endl;
  return 0;
}
*/
//....................................................................
