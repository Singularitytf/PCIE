#include <gsl/gsl_integration.h>
#include <iostream>

double f (double x, void * params) {
  double alpha = *(double *) params;
  double f = log(alpha*x) / sqrt(x);
  return f;
}
 
int main (void)
{
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  
  double result, error;
  double expected = -4.0;
  double alpha = 1.0;
 
  gsl_function F;
  F.function = &f;
  F.params = &alpha;
 
  gsl_integration_qag (&F, 0, 1, 1e-7, 1e-7, 1000, 6,
                        w, &result, &error); 
 
  printf ("result          = % .18f\n", result);
  printf ("exact result    = % .18f\n", expected);
  printf ("estimated error = % .18f\n", error);
  printf ("actual error    = % .18f\n", result - expected);
  printf ("intervals       =  %ld\n", w->size);
 
  return 0;
}
