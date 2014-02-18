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

#ifndef UQ_DIRICHLET_VECTOR_RV_H
#define UQ_DIRICHLET_VECTOR_RV_H

#include <queso/VectorRV.h>
#include <queso/VectorSpace.h>
#include <queso/JointPdf.h>
#include <queso/VectorRealizer.h>
#include <queso/VectorCdf.h>
#include <queso/VectorMdf.h>
#include <queso/SequenceOfVectors.h>
#include <queso/InfoTheory.h>
#include <gsl/gsl_sf_psi.h> // todo: take specificity of gsl_, i.e., make it general (gsl or boost or etc)

namespace QUESO {

//*****************************************************
// Dirichlet class [RV-05]
//*****************************************************
/*!
 * \class DirichletVectorRV
 * \brief A class representing a vector RV constructed via Dirichlet distribution.
 * 
 * This class allows the user to compute the value of a Dirichlet PDF and to generate realizations
 * (samples) from it.\n
 * 
 */


template<class V, class M>
class DirichletVectorRV : public BaseVectorRV<V,M> {
public:
    //! @name Constructor/Destructor methods
  //@{
  //! Default Constructor
  /*! Construct a Dirichlet vector RV with parameters \c alpha > 0  for \c K variates who live in \c imageSet.
   * The constructor will check whether or not the data provided via \c imageSet belongs to [0,+inf], which is
   * a requirement imposed by the Dirichlet distribution. If this condition is not satisfied, an error 
   * message will be displayed and the program will exit. */
  DirichletVectorRV(const char*                  prefix,
                    const VectorSet<V,M>&        imageSet,
                    const V&                     alpha,
                    size_t                       K);
  //! Default Constructor
  /*! Construct a Dirichlet vector RV with parameters  0 < \c alpha < 1  and \c b>0, whose variates live in \c imageSet.
   * The constructor will check whether or not the data provided via \c imageSet belongs to [0,1] and that they sum to one, which
   * is a requirement imposed by the Dirichlet distribution. If this condition is not satisfied, an error 
   * message will be displayed and the program will exit. */
  DirichletVectorRV(const char*                  prefix,
                    const VectorSet<V,M>&        imageSet,
                    const V&                     alpha,
                    double                       r,
                    size_t                       K);
  
  //! Virtual destructor
  virtual ~DirichletVectorRV();
  //@}
  
  //! @name I/O methods
  //@{
  //! TODO: Prints the vector RV.
  /*! \todo: implement me!*/
  void print(std::ostream& os) const;
  //@}

private:
  using BaseVectorRV<V,M>::m_env;
  using BaseVectorRV<V,M>::m_prefix;
  using BaseVectorRV<V,M>::m_imageSet;
  using BaseVectorRV<V,M>::m_pdf;
  using BaseVectorRV<V,M>::m_realizer;
  using BaseVectorRV<V,M>::m_subCdf;
  using BaseVectorRV<V,M>::m_unifiedCdf;
  using BaseVectorRV<V,M>::m_mdf;
};

}  // End namespace QUESO

#endif // UQ_DIRICHLET_VECTOR_RV_H
