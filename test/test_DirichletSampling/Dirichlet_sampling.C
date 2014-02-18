#include <queso/RngGsl.h>
#include <queso/RngBoost.h>
#include <iostream>
#include <fstream>

int main(int argc, char **argv) 
 {

  unsigned int nsample(3);
  QUESO::RngGsl   randGsl(4,0);
  QUESO::RngBoost randBoost(4,0);

  double meanVector[nsample];
  meanVector[0] = 0.2;
  meanVector[1] = 0.3;
  meanVector[2] = 0.5;
  double r(0.3);
  double beta[nsample];

  double a(0.),b(0.),mingamma(1.);
  for(unsigned int k = 0; k < nsample; k++)
  {
    double locmax = (meanVector[k] > 0.5)?meanVector[k]:1. - meanVector[k];
    if(locmax < mingamma)mingamma = locmax;

    a += meanVector[k] * (1. - meanVector[k]);
    b += meanVector[k]*r * std::sqrt(meanVector[k] * (1. - meanVector[k]) );
  }
  double gamma = std::pow(a/b,2) - 1.L;
  if(gamma < mingamma)gamma = mingamma;
  for(unsigned int k = 0; k < nsample; k++)
  {
     beta[k] = meanVector[k] * gamma;
  }

  double res[nsample];

  double gslmean[nsample],boostmean[nsample];
  double gslvar[nsample], boostvar[nsample];
  int return_flag(0);
  for(unsigned int N = 1000; N < 1000000; N *=10)
  {
    for(unsigned int i=0; i < nsample; i++)
    {
        gslmean[i] = 0.L;
        boostmean[i] = 0.L;
        gslvar[i] = 0.L;
        boostvar[i] = 0.L;
    }
    for(unsigned int i = 0; i < N; i++)
    {
      randGsl.dirichletSample(nsample,beta,res);
      for(unsigned int j = 0; j < nsample; j++)
      {
        gslmean[j] += res[j];
      }
      if(i > 0)
      {
        double vf = (double)(i + 1);
        double rnnf = 1.L/ (vf * (vf-1.L));
        for(unsigned int j = 0; j < nsample;j++)gslvar[j] += rnnf * (vf * res[j] - gslmean[j]) * (vf * res[j] - gslmean[j]);
      }

      randBoost.dirichletSample(3,beta,res);
      for(unsigned int j = 0; j < nsample; j++)
      {
        boostmean[j] += res[j];
      }
      if(i > 0)
      {
        double vf = (double)(i + 1);
        double rnnf = 1.L/ (vf * (vf-1.L));
        for(unsigned int j = 0; j < nsample;j++)boostvar[j] += rnnf * (vf * res[j] - boostmean[j]) * (vf * res[j] - boostmean[j]);
      }
    }

    double precision = 1.L/std::sqrt((double)N);
    for(unsigned int j = 0; j < nsample; j++)
    {
      gslmean[j]   /= (double)N;
      gslvar[j]    /= (double)N;
      boostmean[j] /= (double)N;
      boostvar[j]  /= (double)N;
    }

    for(unsigned int j = 0; j < nsample; j++)
    {
      if(std::abs(gslmean[j] - meanVector[j])/meanVector[j] > precision)
      {
        std::cerr << "error in gsl sampling: Monte Carlo precision is " << precision
                  << "error is " << std::abs(gslmean[j] - meanVector[j])/meanVector[j]
                  << std::endl;
        return_flag = 1;
      }
      if(std::abs(boostmean[j] - meanVector[j])/meanVector[j] > precision)
      {
        std::cerr << "error in boost sampling: Monte Carlo precision is " << precision
                  << "error is " << std::abs(boostmean[j] - meanVector[j])/meanVector[j]
                  << std::endl;
        return_flag = 1;
      }
    }
  }

  return return_flag;
}
