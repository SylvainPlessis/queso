//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// QUESO - a library to support the Quantification of Uncertainty
// for Estimation, Simulation and Optimization
//
// Copyright (C) 2008,2009,2010,2011,2012,2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#include <queso/DirichletVectorRealizer.h>
#include <queso/DirichletVectorRV.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

//
#include <cmath>

namespace QUESO {

// Constructor---------------------------------------
template<class V, class M>
DirichletVectorRV<V,M>::DirichletVectorRV(
  const char*                  prefix,
  const VectorSet<V,M>&        imageSet,
  const V&                     alpha,
  size_t                       K)
  :
  BaseVectorRV<V,M>(((std::string)(prefix)+"uni").c_str(),imageSet)
{
///// test size: alpha.sizeLocal()%K == 0 true
///// test dirichlet: alpha[i] > 0 true

/// bla bla
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering DirichletVectorRV<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

// test 1
  UQ_FATAL_TEST_MACRO(alpha.sizeLocal()%K != 0, 
                      m_env.worldRank(),
                      "In DirichletVectorRV<V,M>::constructor()",
                      "invalid input: Dirichlet distribution needs K values");                

// test 2
  for(unsigned int i = 0; i < alpha.sizeLocal(); i++)
  {

     UQ_FATAL_TEST_MACRO(alpha[i] < 0., 
                          m_env.worldRank(),
                         "In DirichletVectorRV<V,M>::constructor()",
                         "invalid input: Dirichlet distribution needs value > 0");                
  }


  
  m_pdf        = new DirichletJointPdf<V,M>(m_prefix.c_str(),
                                            m_imageSet,
                                            alpha,
                                            K);
  m_realizer   = new DirichletVectorRealizer<V,M>(m_prefix.c_str(),
                                                  m_imageSet,
                                                  alpha,
                                                  K);
  m_subCdf     = NULL; // ? sense for Dirichlet
  m_unifiedCdf = NULL; // ? sense for Dirichlet
  m_mdf        = NULL; // ? sense for Dirichlet

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving DirichletVectorRV<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Constructor---------------------------------------
template<class V, class M>
DirichletVectorRV<V,M>::DirichletVectorRV(
  const char*                  prefix,
  const VectorSet<V,M>&        imageSet,
  const V&                     alpha,
  double                       r,
  size_t                       K)
  :
  BaseVectorRV<V,M>(((std::string)(prefix)+"uni").c_str(),imageSet)
{
///// test size: alpha.sizeLocal()%K == 0 true
///// test dirichlet: foreach K-uplet starting at index i, 
/////    alpha_i + alpha_{i+1} + ... + alpha_{i + K - 1} == 1. true


/// bla bla
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering DirichletVectorRV<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

// test 1
  UQ_FATAL_TEST_MACRO(alpha.sizeLocal()%K != 0, 
                      m_env.worldRank(),
                      "In DirichletVectorRV<V,M>::constructor()",
                      "invalid input: Dirichlet distribution needs K values");                

// test 2
  for(unsigned int i = 0; i < alpha.sizeLocal(); i += K)
  {

     double sum(0.L);
     double eps(1e-10L);
     for(unsigned int j = i; j < i + K; j++)
     {
        sum += alpha[j];
     }

     UQ_FATAL_TEST_MACRO(std::abs(sum - 1.) >= eps,
                          m_env.worldRank(),
                         "In DirichletVectorRV<V,M>::constructor()",
                         "invalid input: Dirichlet distribution needs K values to sum to 1");                
  }


  V beta(alpha);

  for(unsigned int i = 0; i < alpha.sizeLocal(); i += K)
  {
    double a(0.),b(0.),mingamma(1e-303);
    for(unsigned int k = 0; k < K; k++)
    {
       double locmax = (alpha[i + k] > 0.5)?alpha[i + k]:1. - alpha[i+k];
       if(locmax < mingamma)mingamma = locmax;

       a += alpha[i+k] * (1. - alpha[i+k]);
       b += alpha[i+k]*r * std::sqrt(alpha[i+k] * (1. - alpha[i+k]) );
    }
    double gamma = std::pow(a/b,2) - 1.L;
    if(gamma < mingamma)gamma = mingamma;
    for(unsigned int k = 0; k < K; k++)
    {
       beta[i+k] = alpha[i+k] * gamma;
    }
   
  }
  
  m_pdf        = new DirichletJointPdf<V,M>(m_prefix.c_str(),
                                            m_imageSet,
                                            alpha,
                                            K);
  m_realizer   = new DirichletVectorRealizer<V,M>(m_prefix.c_str(),
                                                    m_imageSet,
                                                    beta,
                                                    K);
  m_subCdf     = NULL; // ? sense for Dirichlet
  m_unifiedCdf = NULL; // ? sense for Dirichlet
  m_mdf        = NULL; // ? sense for Dirichlet

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving DirichletVectorRV<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor ---------------------------------------
template<class V, class M>
DirichletVectorRV<V,M>::~DirichletVectorRV()
{
  delete m_mdf;
  delete m_unifiedCdf;
  delete m_subCdf;
  delete m_realizer;
  delete m_pdf;
}
// I/O methods---------------------------------------
template <class V, class M>
void
DirichletVectorRV<V,M>::print(std::ostream& os) const
{
  os << "DirichletVectorRV<V,M>::print() says, 'Please implement me.'" << std::endl;
  return;
}

}  // End namespace QUESO

template class QUESO::DirichletVectorRV<QUESO::GslVector, QUESO::GslMatrix>;
