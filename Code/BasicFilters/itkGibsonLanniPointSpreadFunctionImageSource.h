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
#ifndef __itkGibsonLanniPointSpreadFunctionImageSource_h
#define __itkGibsonLanniPointSpreadFunctionImageSource_h

#include <complex>
#include <vector>

#include "itkOPDBasedWidefieldMicroscopePointSpreadFunctionImageSource.h"
#include "itkOPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand.h"
#include "itkNumericTraits.h"

namespace itk
{
namespace Functor
{

class GibsonLanniPointSpreadFunctionIntegrand :
    public OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand< double >
{
public:
  typedef OPDBasedWidefieldMicroscopePointSpreadFunctionIntegrand< double > Superclass;
  typedef Superclass::ComplexType ComplexType;

  ComplexType operator()(double rho) const
  {
    double r = sqrt(m_X*m_X + m_Y*m_Y) / (0.160 + m_Z);
    double bessel = j0(m_K * m_A * rho * r);

    return bessel * exp(ComplexType(0.0, 1.0) * this->OPD(rho, this->m_Z) * m_K) * rho;
  }

};

} // end namespace Functor


/** \class GibsonLanniPointSpreadFunctionImageSource
 * \brief Generate a synthetic point-spread function according to the
 * Gibson-Lanni model.
 *
 * The Gibson-Lanni point-spread function model takes into account optical
 * path differences from the design conditions of an objective in a
 * widefield fluorescence microscope. This image source generates images
 * according to this model. IMPORTANT: Please pay attention to the units
 * each method expects. Some take nanometers, some take micrometers, and some
 * take millimeters.
 *
 * \ingroup DataSources Multithreaded
 */
template< class TOutputImage >
class ITK_EXPORT GibsonLanniPointSpreadFunctionImageSource :
  public OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource< TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef GibsonLanniPointSpreadFunctionImageSource                                 Self;
  typedef OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource< TOutputImage > Superclass;
  typedef SmartPointer< Self >                                                      Pointer;
  typedef SmartPointer< const Self >                                                ConstPointer;

  /** Typedef for the output image PixelType. */
  typedef TOutputImage                             OutputImageType;
  typedef typename OutputImageType::PixelType      PixelType;
  typedef typename OutputImageType::IndexType      IndexType;
  typedef typename OutputImageType::RegionType     RegionType;
  typedef typename OutputImageType::PointType      PointType;
  typedef typename OutputImageType::PointValueType PointValueType;
  typedef typename OutputImageType::SpacingType    SpacingType;
  typedef typename OutputImageType::SizeType       SizeType;
  typedef typename OutputImageType::SizeValueType  SizeValueType;


  itkStaticConstMacro(ImageDimension, unsigned int,
		      TOutputImage::ImageDimension);

  /** Typedef for complex type. */
  typedef std::complex<double> ComplexType;

  /** Typedef for functor. */
  typedef Functor::GibsonLanniPointSpreadFunctionIntegrand FunctorType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(GibsonLanniPointSpreadFunctionImageSource, ParametricImageSource);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

protected:
  GibsonLanniPointSpreadFunctionImageSource();
  ~GibsonLanniPointSpreadFunctionImageSource();
  void PrintSelf(std::ostream& os, Indent indent) const;

  void BeforeThreadedGenerateData();
  virtual void ThreadedGenerateData(const RegionType& outputRegionForThread, ThreadIdType threadId );

  /** Computes the light intensity at a specified point. */
  double ComputeSampleValue(PointType& point, int threadId);

private:
  GibsonLanniPointSpreadFunctionImageSource(const GibsonLanniPointSpreadFunctionImageSource&); //purposely not implemented
  void operator=(const GibsonLanniPointSpreadFunctionImageSource&); //purposely not implemented

  // An array of functors, one for each thread.
  std::vector<FunctorType *> m_IntegrandFunctors;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGibsonLanniPointSpreadFunctionImageSource.hxx"
#endif

#endif
