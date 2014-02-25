#include <queso/GenericScalarFunction.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/ScalarFunction.h>
#include <queso/VectorSet.h>
#include <queso/DirichletVectorRV.h>

template<class V, class M>
class Likelihood : public QUESO::BaseScalarFunction<V, M>
{
public:

  Likelihood(const char * prefix, const QUESO::VectorSet<V, M> & domain):
        QUESO::BaseScalarFunction<V, M>(prefix,domain),K(3)
  {
    // Setup here
    br.push_back(0.5);
    br.push_back(0.3);
    br.push_back(0.2);
    dbr.push_back(0.1);
    dbr.push_back(0.05);
    dbr.push_back(0.01);
  }

  virtual ~Likelihood()
  {
    // Deconstruct here
  }

  virtual double lnValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const
  {
    // 1) Run the forward code
    // 2) Compare to data, hypothesises are:
    //     - uncorrelated (clearly wrong...)
    //     - gaussian likelihood (not sure what to think of it)
       double misfit(0.L);
       for(unsigned int i = 0; i < K; i++)
       {
         misfit += (domainVector[i] - br[i]) * (domainVector[i] - br[i]) / (dbr[i] * dbr[i]);
       }
    // 3) Return log likelihood

#ifdef QUESO_EXPECTS_LN_LIKELIHOOD_INSTEAD_OF_MINUS_2_LN
     return -0.5 * misfit;
#else
     return misfit;
#endif
  }

  virtual double actualValue(const V & domainVector, const V * domainDirection,
      V * gradVector, M * hessianMatrix, V * hessianEffect) const
  {
    // this should return exp(lnValue)
      return std::exp(this->lnValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));
  }

private:
  // Should store the observed data here.
  std::vector<double> br;
  std::vector<double> dbr;
  const unsigned int K;
};

int main(int argc, char ** argv) {
  MPI_Init(&argc, &argv);

  // Step 0 of 5: Set up environment
  QUESO::FullEnvironment env(MPI_COMM_WORLD, argv[1], "", NULL);

  // Step 1 of 5: Instantiate the parameter space
  QUESO::VectorSpace<QUESO::GslVector, QUESO::GslMatrix> paramSpace(env,
      "param_", 3, NULL);

  // Step 2 of 5: Set up the prior
  
  QUESO::GslVector params(paramSpace.zeroVector());
  QUESO::GslVector paramMins(paramSpace.zeroVector());
  QUESO::GslVector paramMaxs(paramSpace.zeroVector());
  paramMins[0] = 0.;
  paramMins[1] = 0.;
  paramMins[2] = 0.;
  paramMaxs[0] = 1.;
  paramMaxs[1] = 1.;
  paramMaxs[2] = 1.;
  params[0] = 0.5L;
  params[1] = 0.3L;
  params[2] = 0.2L;
  size_t K(3);
  double r(0.3);

  QUESO::BoxSubset<QUESO::GslVector, QUESO::GslMatrix> paramDomain("param_",
      paramSpace, paramMins, paramMaxs);

  QUESO::DirichletVectorRV<QUESO::GslVector, QUESO::GslMatrix> priorRv("prior_",
      paramDomain,params,r,K);

  // Step 3 of 5: Set up the likelihood
  Likelihood<QUESO::GslVector, QUESO::GslMatrix> llhd("llhd_", paramDomain);
  
  // Step 4 of 5: Instantiate the inverse problem
  QUESO::GenericVectorRV<QUESO::GslVector, QUESO::GslMatrix>
    postRv("post_", paramSpace);
  QUESO::StatisticalInverseProblem<QUESO::GslVector, QUESO::GslMatrix>
    ip("", NULL, priorRv, llhd, postRv);

  // Step 5 of 5: Solve the inverse problem
  QUESO::GslVector paramInitials(paramSpace.zeroVector());

  // Initial condition of the chain
  paramInitials[0] = 0.5;
  paramInitials[1] = 0.3;
  paramInitials[2] = 0.2;
  
  QUESO::GslMatrix proposalCovMatrix(paramSpace.zeroVector());

  for (unsigned int i = 0; i < 2; i++) {
    // Might need to tweak this
    proposalCovMatrix(i, i) = 0.1;
  }
  
  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

  MPI_Finalize();
 
  return 0;
}
