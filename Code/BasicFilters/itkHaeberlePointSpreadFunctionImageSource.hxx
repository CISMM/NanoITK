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
#ifndef __itkHaeberlePointSpreadFunctionImageSource_hxx
#define __itkHaeberlePointSpreadFunctionImageSource_hxx

#include "itkHaeberlePointSpreadFunctionImageSource.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMath.h"
#include "itkObjectFactory.h"
#include "itkProgressReporter.h"

#include <algorithm>

namespace itk
{

//----------------------------------------------------------------------------
template <class TOutputImage>
HaeberlePointSpreadFunctionImageSource<TOutputImage>
::HaeberlePointSpreadFunctionImageSource()
{
}


//----------------------------------------------------------------------------
template <class TOutputImage>
HaeberlePointSpreadFunctionImageSource<TOutputImage>
::~HaeberlePointSpreadFunctionImageSource()
{
}


//----------------------------------------------------------------------------
template <class TOutputImage>
void
HaeberlePointSpreadFunctionImageSource<TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}


//----------------------------------------------------------------------------
template <typename TOutputImage>
void
HaeberlePointSpreadFunctionImageSource< TOutputImage >
::BeforeThreadedGenerateData()
{
  // Initialize functors
  for ( size_t i = 0; i < m_I0illFunctors.size(); i++ )
    {
    delete m_I0illFunctors[i];
    m_I0illFunctors[i] = NULL;
    delete m_I1illFunctors[i];
    m_I1illFunctors[i] = NULL;
    delete m_I2illFunctors[i];
    m_I2illFunctors[i] = NULL;
    }

  m_I0illFunctors.resize( this->GetNumberOfThreads(), NULL );
  m_I1illFunctors.resize( this->GetNumberOfThreads(), NULL );
  m_I2illFunctors.resize( this->GetNumberOfThreads(), NULL );
  for ( int i = 0; i < this->GetNumberOfThreads(); ++i )
    {
    m_I0illFunctors[i] = new FunctorTypeI0ill();
    m_I0illFunctors[i]->CopySettings( this );
    m_I1illFunctors[i] = new FunctorTypeI1ill();
    m_I1illFunctors[i]->CopySettings( this );
    m_I2illFunctors[i] = new FunctorTypeI2ill();
    m_I2illFunctors[i]->CopySettings( this );
    }
}


//----------------------------------------------------------------------------
template <typename TOutputImage>
void
HaeberlePointSpreadFunctionImageSource<TOutputImage>
::ThreadedGenerateData(const RegionType& outputRegionForThread, ThreadIdType threadId )
{
  // Support progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  typename TOutputImage::Pointer image = this->GetOutput(0);

  ImageRegionIteratorWithIndex<OutputImageType> it(image, outputRegionForThread);

  for (; !it.IsAtEnd(); ++it)
    {
    IndexType index = it.GetIndex();
    PointType point;
    image->TransformIndexToPhysicalPoint(index, point);

    it.Set( this->ComputeSampleValue( point, threadId ) );
    progress.CompletedPixel();
    }
}


//----------------------------------------------------------------------------
template <typename TOutputImage>
double
HaeberlePointSpreadFunctionImageSource<TOutputImage>
::ComputeSampleValue(const PointType& point, ThreadIdType threadId)
{
  PixelType px = point[0] * 1e-9;
  PixelType py = point[1] * 1e-9;
  PixelType pz = point[2] * 1e-9;

  /* Compute terms that are independent of terms within the integral. */
  double mag = this->m_Magnification;

  // We have to convert to coordinates of the detector points
  m_I0illFunctors[threadId]->m_X =
    m_I1illFunctors[threadId]->m_X =
    m_I2illFunctors[threadId]->m_X = px * mag;
  m_I0illFunctors[threadId]->m_Y =
    m_I1illFunctors[threadId]->m_Y =
    m_I2illFunctors[threadId]->m_Y = py * mag;
  m_I0illFunctors[threadId]->m_Z =
    m_I1illFunctors[threadId]->m_Z =
    m_I2illFunctors[threadId]->m_Z = pz; // No conversion needed


  // Using the alpha suggested by the paper (asin(NA / actual
  // immersion oil RI) yields an uncomputable angle of incidence in
  // the specimen layer. Taking the minimum angle from this angle and
  // the angle of incidence in the immersion oil layer that yields a
  // 90 degrees angle of incidence in the specimen layer makes the
  // integral computable.
  double alpha = std::min(asin(this->m_NumericalAperture / this->m_ActualImmersionOilRefractiveIndex),
                          asin(this->m_ActualSpecimenLayerRefractiveIndex / this->m_ActualImmersionOilRefractiveIndex));
  ComplexType I0ill = this->Integrate(*m_I0illFunctors[threadId], 0.0, alpha, 20);
  ComplexType I1ill = this->Integrate(*m_I1illFunctors[threadId], 0.0, alpha, 20);
  ComplexType I2ill = this->Integrate(*m_I2illFunctors[threadId], 0.0, alpha, 20);

  // Return a weighted sum of the squared magnitudes of the three
  // integral values.
  return static_cast<PixelType>( norm(I0ill) + 2.0*norm(I1ill) + norm(I2ill) );
}


} // end namespace itk

#endif
