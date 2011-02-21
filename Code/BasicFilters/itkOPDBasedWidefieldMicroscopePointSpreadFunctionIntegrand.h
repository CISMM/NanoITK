/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkOPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand_h
#define __itkOPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand_h

#include <complex>

namespace itk
{
namespace Functor
{

/** \class OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand
 * \brief Functor class that defines a method for computing the
 * optical path difference (OPD) in several widefield point-spread
 * function models.
 */

template< class TParamRep >
class ITK_EXPORT OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand
{
public:
  typedef TParamRep ParameterRepType;
  typedef std::complex< ParameterRepType > ComplexType;

  OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand();
  virtual ~OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand();

  template< class TSource >
  void CopySettings(const TSource* source);

  /** Point-spread function model parameters. */
  ParameterRepType m_EmissionWavelength;
  ParameterRepType m_NumericalAperture;
  ParameterRepType m_Magnification;
  ParameterRepType m_DesignCoverSlipRefractiveIndex;
  ParameterRepType m_ActualCoverSlipRefractiveIndex;
  ParameterRepType m_DesignCoverSlipThickness;
  ParameterRepType m_ActualCoverSlipThickness;
  ParameterRepType m_DesignImmersionOilRefractiveIndex;
  ParameterRepType m_ActualImmersionOilRefractiveIndex;
  ParameterRepType m_DesignImmersionOilThickness;
  ParameterRepType m_DesignSpecimenLayerRefractiveIndex;
  ParameterRepType m_ActualSpecimenLayerRefractiveIndex;
  ParameterRepType m_ActualPointSourceDepthInSpecimenLayer;

  /** Precomputed values. */
  ParameterRepType m_K; // Wavenumber
  ParameterRepType m_A; // Radius of projection of the limiting
                        // aperture onto the back focal plane of the
                        // objective lens

  /** Sample coordinates in physical space. */
  ParameterRepType m_X;
  ParameterRepType m_Y;
  ParameterRepType m_Z;

protected:
  /** Computes the optical path difference for a ray terminating at
  *   a normalized distance rho from the center of the back focal
  *   plane aperture. */
  ComplexType OPD(ParameterRepType rho, ParameterRepType dz) const;

  /** Common terms for computing the optical path difference term in
   *  point-spread function models descended from this class. */
  ComplexType OPDTerm(ParameterRepType rho, ParameterRepType n, ParameterRepType t) const;

};

} // end namespace Functor

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkOPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand.txx"
#endif

#endif // __itkOPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand_h
