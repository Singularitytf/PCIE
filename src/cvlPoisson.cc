#include "cvlPoisson.hh"

//....................................................................
cvlPoisson::cvlPoisson(double fsrate,
                       double fbkg_mean,
                       double fbkg_sigma)
  : s_rate(fsrate),
    bkg_mean(fbkg_mean),
    bkg_sigma(fbkg_sigma)
{
    
}



cvlPoisson::~cvlPoisson()
{
  
}

//....................................................................
// Set Functions;
void cvlPoisson::setSRate(double fsrate){s_rate = fsrate;}
void cvlPoisson::setBkgMean(double fbkg_mean){bkg_mean = fbkg_mean;}
void cvlPoisson::setBkgSigma(double fbkg_sigma){bkg_sigma = fbkg_sigma;}

//....................................................................
// Get Functions;
double cvlPoisson::getSRate() {return this->s_rate;}
double cvlPoisson::getBkgMean() {return this->bkg_mean;}
double cvlPoisson::getBkgSigma() {return this->bkg_sigma;}

//....................................................................

double cvlPoisson::pmf(int n0){
  double paras[4] = {bkg_mean, bkg_sigma, s_rate, n0};
  gsl_integration_workspace * w
    = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = this->fcvlPoisson;
  F.params = paras;
  gsl_integration_qag (&F, 0, bkg_mean+5*bkg_sigma, esp, esp, 1000, 6,
                        w, &result, &error); 
  gsl_integration_workspace_free(w);
  // if ((bkg_mean - 5*bkg_sigma) <= 0)
  //   return func.Integral(0, bkg_mean+5*bkg_sigma);
  // else
  //   return func.Integral(bkg_mean-5*bkg_sigma, bkg_mean+5*bkg_sigma);
  return result;
}

//....................................................................

double cvlPoisson::fcvlPoisson(double bkg, void *fpara){
    // printf("%f\t%f\t%f\t%d\n", bkg_mean, bkg_sigma, mu, n0);
    double *para = (double*)fpara;
    return TMath::Poisson((int)para[3], para[2]+(bkg))
      *TMath::Gaus(bkg, para[0], para[1], kTRUE);
}

double cvlPoisson::cdf(int n0){
  double temp = 0.;
  for(int i=0; i<=n0; i++) {temp += this->pmf(i);}
  return temp;
}

//....................................................................

double cvlPoisson::findBestMu(int n0){
  // Get n0-bkg for init best mu;
  // double BestMu =
  //   ((n0 - this->bkg_mean)>0) ? (n0-this->bkg_mean) : 0;
  // double mu_iter0 = 0;
  // double mu_iter1 = 0;
  // bool findBestMu = false;
  // cvlPoisson tmp = cvlPoisson(BestMu, this->bkg_mean, this->bkg_sigma);
  // for (double i=0.;i<=BestMu+1;i+=0.001){
  //   tmp.setSRate(i);
  //   mu_iter0 = tmp.pmf(n0);
  //   //cout << mu_iter << endl;
  //   if (mu_iter0 > mu_iter1) mu_iter1 = mu_iter0;
  //   else {
  //     findBestMu = true;
  //     BestMu = i;
  //     break;
  //   }
  // }
  // // 检查是否找到了Bestfit
  // if (findBestMu) return BestMu;
  // else{
  //   std::cout << "Faild to find Best mu!!!" << std::endl;
  //   return -1.0;
  // }
  return (n0-bkg_mean) > 0. ? (n0-bkg_mean) : 0.;
}

//....................................................................
//TestPart

/*
void pdf(){

  double n_l[15];
  double pois_p[15];
  double cvlpois_p[15];
  cvlPoisson cp = cvlPoisson(1., 3., 0.8);
  for (int i=0; i<15; i++){
    n_l[i] = i;
    pois_p[i] = TMath::Poisson(i, 4);
    cvlpois_p[i] = cp.pmf(i);
    //cout << cvlpois_p[i] << endl;
  }
  //  cout << cvlpois_p[3] << endl;
  //cout << cvlPoisson(1., 3., 0.2).pmf(10) << endl;
  TCanvas *c1 = new TCanvas ("c1","Pois vs. cvlPoisson under 80% uc",
                             200,10,600,400);
  auto gh = new TGraph(15, n_l, pois_p);
  gh->SetTitle("Pois vs. cvlPoisson under 80% Uncertainty.");
  gh->Draw("A*");
  auto gh2 = new TGraph(15, n_l, cvlpois_p);
  gh2->Draw("same");
  cout << cp.cdf(15) << endl;

  cvlPoisson cp = cvlPoisson(1., 3., 0.8);
  std::cout << cp.findBestMu(10)<< std::endl;
}


int main(){
  //printf("%.12f", poissonPmf(3,3));
}
*/
