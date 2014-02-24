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

#ifndef UQ_DIRICHLET_JOINT_PROB_DENSITY_H
#define UQ_DIRICHLET_JOINT_PROB_DENSITY_H

#include <cmath>

#include <boost/math/special_functions.hpp> // for Boost isnan. Note parentheses are important in function call.

#include <queso/JointPdf.h>
#include <queso/Environment.h>
#include <queso/ScalarFunction.h>
#include <queso/BoxSubset.h>

namespace QUESO {

//*****************************************************
// Classes to accommodate a probability density.
//*****************************************************

//*****************************************************
// Dirichlet class [PDF-00]
//*****************************************************
/*! \file DirichletJointPdf.h
 * \brief Classes to accommodate the Dirichlet probability density.
 * 
 * \class DirichletJointPdf
 *
 * */

template<class V, class M>
class DirichletJointPdf : public BaseJointPdf<V,M> {
public:
  //! @name Constructor/Destructor methods
  //@{ 
  //! Default constructor.
  /*! Instantiates an object of the class, i.e. a scalar function, given a prefix and its domain.*/
  DirichletJointPdf(const char*                  prefix,
		    const VectorSet<V,M>& domainSet,
                    const V&                     a,
                    const size_t&                K);
  //! Destructor
  virtual ~DirichletJointPdf();
  //@}

  //! @name Mathematical methods
  //@{  
  //! Actual value of the PDF (scalar function).
  double actualValue                    (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  
  //! Logarithm of the value of the function.
  double lnValue                        (const V& domainVector, const V* domainDirection, V* gradVector, M* hessianMatrix, V* hessianEffect) const;
  
  //! Sets a value to be used in the normalization style (stored in the protected attribute m_normalizationStyle.) 
  void   setNormalizationStyle          (unsigned int value) const;
  
  //! Sets a logarithmic value to be used in the normalization factor (stored in the protected attribute m_normalizationStyle.) 
  void   setLogOfNormalizationFactor    (double value) const;
  
  //! Computes the logarithm of the normalization factor. See template specialization.
  double computeLogOfNormalizationFactor(unsigned int numSamples, bool m_logOfNormalizationFactor) const = 0;
  
  //@}
protected:
  using BaseScalarFunction<V,M>::m_env;
  using BaseScalarFunction<V,M>::m_prefix;
  using BaseScalarFunction<V,M>::m_domainSet;
  using BaseJointPdf<V,M>::m_normalizationStyle;
  using BaseJointPdf<V,M>::m_logOfNormalizationFactor;

  V      m_a;
  size_t m_K;

};

}  // End namespace QUESO

#endif // UQ_DIRICHLET_JOINT_PROB_DENSITY_H
