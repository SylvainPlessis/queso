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

#include <queso/DirichletJointPdf.h>
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>

namespace QUESO {

// Constructor -------------------------------------
template<class V,class M>
DirichletJointPdf<V,M>::DirichletJointPdf(
  const char*                  prefix,
  const VectorSet<V,M>& domainSet,
  const V&                     a,
  const size_t&                K)
  :
  BaseJointPdf<V,M>(((std::string)(prefix)+"uni").c_str(),domainSet),
  m_a(a),
  m_K(K)
{
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Entering DirichletJointPdf<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 54)) {
    *m_env.subDisplayFile() << "Leaving DirichletJointPdf<V,M>::constructor()"
                            << ": prefix = " << m_prefix
                            << std::endl;
  }
}
// Destructor --------------------------------------
template<class V,class M>
DirichletJointPdf<V,M>::~DirichletJointPdf()
{
}
// Math methods-------------------------------------
template<class V, class M>
double
DirichletJointPdf<V,M>::actualValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  UQ_FATAL_TEST_MACRO(domainVector.sizeLocal() != this->m_domainSet.vectorSpace().dimLocal(),
                      m_env.worldRank(),
                      "DirichletJointPdf<V,M>::actualValue()",
                      "invalid input");

  UQ_FATAL_TEST_MACRO((domainDirection || gradVector || hessianMatrix || hessianEffect),
                      m_env.worldRank(),
                      "DirichletJointPdf<V,M>::actualValue()",
                      "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");

  UQ_FATAL_TEST_MACRO(domainVector.sizeLocal() != m_K,
                      m_env.worldRank(),
                      "DirichletJointPdf<V,M>::actualValue()",
                      "This joint pdf do not support multiple Dirichlet");

  double result = m_env.basicPdfs()->DirichletPdfActualValue(domainVector[i],m_a[i],m_K);
    if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 99)) {
      *m_env.subDisplayFile() << "In DirichletJointPdf<V,M>::lnValue()"
                              << ", m_normalizationStyle = "      << m_normalizationStyle
                              << ": domainVector[" << i << "] = " << domainVector[i]
                              << ", m_a[" << i << "] = "          << m_a[i]
                              << ", m_b[" << i << "] = "          << m_b[i]
                              << ", log(pdf)= "                   << aux
                              << std::endl;
    }
  
  return result,
}
//--------------------------------------------------
template<class V, class M>
double
DirichletJointPdf<V,M>::lnValue(
  const V& domainVector,
  const V* domainDirection,
        V* gradVector,
        M* hessianMatrix,
        V* hessianEffect) const
{
  UQ_FATAL_TEST_MACRO((domainDirection || gradVector || hessianMatrix || hessianEffect),
                      m_env.worldRank(),
                      "DirichletJointPdf<V,M>::lnValue()",
                      "incomplete code for gradVector, hessianMatrix and hessianEffect calculations");
  return log(this->actualValue(domainVector,domainDirection,gradVector,hessianMatrix,hessianEffect));
}
//--------------------------------------------------
template<class V, class M>
double
DirichletJointPdf<V,M>::computeLogOfNormalizationFactor(unsigned int numSamples, bool updateFactorInternally) const
{
  double value = 0.;

  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Entering DirichletJointPdf<V,M>::computeLogOfNormalizationFactor()"
                            << std::endl;
  }
  value = BaseJointPdf<V,M>::commonComputeLogOfNormalizationFactor(numSamples, updateFactorInternally);
  if ((m_env.subDisplayFile()) && (m_env.displayVerbosity() >= 2)) {
    *m_env.subDisplayFile() << "Leaving DirichletJointPdf<V,M>::computeLogOfNormalizationFactor()"
                            << ", m_logOfNormalizationFactor = " << m_logOfNormalizationFactor
                            << std::endl;
  }

  return value;
}

}  // End namespace QUESO

template class QUESO::DirichletJointPdf<QUESO::GslVector, QUESO::GslMatrix>;
