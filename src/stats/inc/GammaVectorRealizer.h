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
// 
// $Id$
//
//--------------------------------------------------------------------------

#ifndef UQ_GAMMA_REALIZER_H
#define UQ_GAMMA_REALIZER_H

#include <queso/VectorRealizer.h>
#include <queso/VectorSequence.h>
#include <queso/Environment.h>
#include <math.h>

namespace QUESO {

//*****************************************************
// Gamma class [R-06]
//*****************************************************
/*! 
 * \class GammaVectorRealizer
 * \brief A class for handling sampling from a Gamma probability density distribution.
 *
 * This class handles sampling from a Gamma probability density distribution, of 
 * parameters \c a and \c b.*/
template<class V, class M>
class GammaVectorRealizer : public BaseVectorRealizer<V,M> {
public:
  
  //! @name Constructor/Destructor methods
  //@{
  //! Constructor
  /*! Constructs a new object, given a prefix, the image set of the vector realizer, and the
   * Gamma distribution parameters \c a and \c b, which are assigned to private attributes
   * m_alpha and m_beta.  */
  GammaVectorRealizer(const char*                  prefix,
                             const VectorSet<V,M>& unifiedImageSet,
                             const V&                     a,
                             const V&                     b);
  
  //! Destructor
 ~GammaVectorRealizer();
  //@}
 
   //! @name Realization-related methods
  //@{
  //! Draws a realization.
  /*! This function draws a realization of a Gamma distribution and saves it in \c nextValues. 
   * It internally checks whether the image set, where the realization should be drawn, belongs 
   * to the interval (0, infinity) - which is the range where Gamma distribution is defined over. */ 
  void realization(V& nextValues) const;

private:
  using BaseVectorRealizer<V,M>::m_env;
  using BaseVectorRealizer<V,M>::m_prefix;
  using BaseVectorRealizer<V,M>::m_unifiedImageSet;
  using BaseVectorRealizer<V,M>::m_subPeriod;

  V m_a;
  V m_b;
};

}  // End namespace QUESO

#endif // UQ_GAMMA_REALIZER_H
