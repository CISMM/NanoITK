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
#ifndef _itkBeadSpreadFunctionImageSource2_hxx
#define _itkBeadSpreadFunctionImageSource2_hxx

#include "itkBeadSpreadFunctionImageSource2.h"


namespace itk
{

template< class TOutputImage >
BeadSpreadFunctionImageSource2< TOutputImage >
::BeadSpreadFunctionImageSource2()
{
  m_BeadRadius = 10.0;
  m_BeadCenter.Fill( 0.0 );
  m_BeadSampleSpacing.Fill( 10.0 );

  m_BeadShapeModified = true;

  m_IntensityShift = 0.0;
  m_IntensityScale = 1.0;

  m_KernelSource = NULL;

  m_Convolver = ConvolverType::New();
  m_Convolver->SetInput( PointSetType::New() );

  m_RescaleFilter = RescaleImageFilterType::New();
  m_RescaleFilter->SetInput(m_Convolver->GetOutput());

  m_ModifiedEventCommand = MemberCommandType::New();
  m_ModifiedEventCommand->SetCallbackFunction(this, &Self::KernelModified);
  m_ObserverTag = 0;
}


template< class TOutputImage >
BeadSpreadFunctionImageSource2< TOutputImage >
::~BeadSpreadFunctionImageSource2()
{
  this->m_KernelSource->RemoveObserver( this->m_ObserverTag );
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource2< TOutputImage >
::KernelModified()
{
  this->Modified();
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource2< TOutputImage >
::SetSize(const SizeType& size)
{
  if ( size != m_Size )
    {
    m_Size = size;
    m_Convolver->SetSize( size );
    this->Modified();
    }
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource2< TOutputImage >
::SetSpacing(const SpacingType& spacing)
{
  if ( spacing != m_Convolver->GetSpacing() )
    {
    this->Modified();
    }
  m_Convolver->SetSpacing( spacing );
}


template< class TOutputImage >
const typename BeadSpreadFunctionImageSource2< TOutputImage >::SpacingType&
BeadSpreadFunctionImageSource2< TOutputImage >
::GetSpacing() const
{
  return m_Convolver->GetSpacing();
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource2< TOutputImage >
::SetOrigin(const PointType& origin)
{
  if ( origin != m_Convolver->GetOrigin() )
    {
    this->Modified();
    }
  m_Convolver->SetOrigin( origin );
}


template< class TOutputImage >
const typename BeadSpreadFunctionImageSource2< TOutputImage >::PointType&
BeadSpreadFunctionImageSource2< TOutputImage >
::GetOrigin() const
{
  return m_Convolver->GetOrigin();
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource2< TOutputImage >
::SetBeadRadius( double radius )
{
  if ( radius != m_BeadRadius )
    {
    this->m_BeadRadius = radius;
    this->m_BeadShapeModified = true;
    this->Modified();
    }
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource2< TOutputImage >
::SetShearX(double shear)
{
#if 0
  if ( shear != m_Convolver->GetShearX() )
    {
    m_Convolver->SetShearX( shear );
    this->Modified();
    }
#endif
}


template< class TOutputImage >
double
BeadSpreadFunctionImageSource2< TOutputImage >
::GetShearX() const
{
#if 0
  return m_Convolver->GetShearX();
#endif
  return 0.0;
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource2< TOutputImage >
::SetShearY(double shear)
{
#if 0
  if (shear != m_Convolver->GetShearY())
    {
    m_Convolver->SetShearY( shear );
    this->Modified();
    }
#endif
}


template< class TOutputImage >
double
BeadSpreadFunctionImageSource2< TOutputImage >
::GetShearY() const
{
#if 0
  return m_Convolver->GetShearY();
#endif
  return 0.0;
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource2< TOutputImage >
::SetBeadCenter( const PointType& center )
{
  if ( this->m_BeadCenter != center )
    {
    this->m_BeadCenter = center;
    this->m_BeadShapeModified = true;
    this->Modified();
    }
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource2< TOutputImage >
::SetKernelSource( KernelImageSourceType* source )
{
  if ( this->m_KernelSource != source )
    {
    if ( this->m_KernelSource )
      {
      this->m_KernelSource->RemoveObserver( this->m_ObserverTag );
      }
    this->m_KernelSource = source;
    this->m_ObserverTag = this->m_KernelSource->
      AddObserver( ModifiedEvent() , m_ModifiedEventCommand );
    this->Modified();
    }

}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource2< TOutputImage >
::SetParameters(const ParametersType& parameters)
{
  int index = 0;

  // The first parameters are bead-spread function parameters
  SpacingType spacing;
  for ( int i = 0; i < ImageDimension; i++ )
    {
    spacing[i] = parameters[index++];
    }
  this->SetSpacing( spacing );

  this->SetBeadRadius( parameters[index++] );

  PointType center;
  for ( int i = 0; i < ImageDimension; i++ )
    {
    center[i] = parameters[index++];
    }
  this->SetBeadCenter( center );

  this->SetShearX( parameters[index++] );
  this->SetShearY( parameters[index++] );

  this->SetIntensityShift( parameters[index++] );
  this->SetIntensityScale( parameters[index++] );

  // The last parameters go to the kernel source
  ParametersType kernelParameters(this->m_KernelSource->GetNumberOfParameters());
  for ( unsigned int i = 0; i < kernelParameters.GetSize(); i++ )
    {
    kernelParameters[i] = parameters[index++];
    }

  this->m_KernelSource->SetParameters( kernelParameters );
}


template< class TOutputImage >
typename BeadSpreadFunctionImageSource2< TOutputImage >::ParametersType
BeadSpreadFunctionImageSource2< TOutputImage >
::GetParameters() const
{
  ParametersType parameters( GetNumberOfParameters() );
  int index = 0;

  // The first parameters come from the bead-spread function
  const SpacingType spacing = this->GetSpacing();
  for ( int i = 0; i < ImageDimension; i++ )
    {
    parameters[index++] = spacing[i];
    }

  parameters[index++] = this->GetBeadRadius();

  const PointType beadCenter = this->GetBeadCenter();
  for ( int i = 0; i < ImageDimension; i++ )
    {
    parameters[index++] = beadCenter[i];
    }

  parameters[index++] = this->GetShearX();
  parameters[index++] = this->GetShearY();
  parameters[index++] = this->GetIntensityShift();
  parameters[index++] = this->GetIntensityScale();

  // The last parameters come from the kernel source
  ParametersType kernelParameters = this->m_KernelSource->GetParameters();
  for ( unsigned int i = 0; i < kernelParameters.GetSize(); i++ )
    {
    parameters[index++] = kernelParameters[i];
    }

  return parameters;
}


template< class TOutputImage >
unsigned int
BeadSpreadFunctionImageSource2< TOutputImage >
::GetNumberOfParameters() const
{
  return this->m_KernelSource->GetNumberOfParameters() +
    this->GetNumberOfBeadSpreadFunctionParameters();
}


template< class TOutputImage >
unsigned int
BeadSpreadFunctionImageSource2< TOutputImage >
::GetNumberOfBeadSpreadFunctionParameters() const
{
  return 2*ImageDimension + 5;
}

template< class TOutputImage >
void
BeadSpreadFunctionImageSource2< TOutputImage >
::SetZCoordinate(unsigned int index, double coordinate)
{
#if 0
  m_Convolver->SetZCoordinate( index, coordinate );
#endif
}


template< class TOutputImage >
double
BeadSpreadFunctionImageSource2< TOutputImage >
::GetZCoordinate(unsigned int index)
{
#if 0
  return m_Convolver->GetZCoordinate( index );
#endif

  return 0.0;
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource2< TOutputImage >
::SetUseCustomZCoordinates(bool use)
{
#if 0
  m_Convolver->SetUseCustomZCoordinates( use );
  this->Modified();
#endif
}


template< class TOutputImage >
bool
BeadSpreadFunctionImageSource2< TOutputImage >
::GetUseCustomZCoordinates()
{
#if 0
  return m_Convolver->GetUseCustomZCoordinates();
#endif

  return false;
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource2< TOutputImage >
::RasterizeShape()
{
  // Rasterize the sphere into a point set
  // First, determine the number of samples to take in each dimension
  IndexType sphereSamplesIndex;
  SizeType  sphereSamplesSize;
  for ( unsigned int dim = 0; dim < ImageDimension; dim++ )
    {
    sphereSamplesIndex[dim] =
      (Math::Ceil<IndexValueType>( -m_BeadRadius / m_BeadSampleSpacing[dim] ) );
    sphereSamplesSize[dim] =
      (Math::Floor<SizeValueType>( ( 2.0 * m_BeadRadius ) / m_BeadSampleSpacing[dim] ) );
    }

  // Then, set up an empty image to conveniently obtain indices and physical
  // point locations of samples during sphere rasterization
  RegionType dummyRegion;
  dummyRegion.SetIndex( sphereSamplesIndex );
  dummyRegion.SetSize( sphereSamplesSize );
  typename OutputImageType::Pointer dummyImage = OutputImageType::New();
  dummyImage->SetRegions( dummyRegion );

  dummyImage->SetSpacing( m_BeadSampleSpacing );

  PointType dummyOrigin; dummyOrigin.Fill( 0.0 );
  dummyImage->SetOrigin( dummyOrigin );

  typename PointSetType::Pointer pointSet = PointSetType::New();

  SizeValueType numSamples = 1;
  for ( unsigned int dim = 0; dim < ImageDimension; dim++ )
    {
    numSamples *= sphereSamplesSize[dim];
    }

  PointType zero; zero.Fill( 0.0 );
  typename PointType::VectorType centerVector = m_BeadCenter - zero;
  double r2 = m_BeadRadius * m_BeadRadius;

  // Iterate over the offsets in the image, convert them to indices,
  // then convert the indices to the physical samples points. With the
  // points in hand, keep only those that are within the sphere.
  typename PointSetType::PointIdentifier ptId = 0;
  for ( SizeValueType offset = 0; offset < numSamples; offset++ )
    {
    IndexType dummyIndex = dummyImage->ComputeIndex( offset );
    PointType candidatePoint;
    dummyImage->TransformIndexToPhysicalPoint( dummyIndex, candidatePoint );
    if ( candidatePoint.SquaredEuclideanDistanceTo( zero ) < r2 )
      {
      pointSet->SetPoint( ptId++, candidatePoint + centerVector );
      }
    }

  m_Convolver->SetInput( pointSet );
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource2< TOutputImage >
::GenerateData()
{
  // Set the PSF sampling spacing and size parameters, and update.
  PointType   psfTableOrigin;
  SizeType    psfTableSize;
  SpacingType psfTableSpacing( 50.0 ); // An arbitrary spacing

  // Determine necessary spatial extent of PSF table.
  PointType minExtent( this->GetOrigin() );
  PointType maxExtent;
  const unsigned int dimensions = itkGetStaticConstMacro( OutputImageDimension );
  for ( unsigned int i = 0; i < dimensions; i++ )
    {
    // First calculate extent of BSF in this dimension.
    maxExtent[i] = static_cast<PointValueType>
      (this->GetSize()[i]-1) * this->GetSpacing()[i] + minExtent[i];

    // Now modify calculated PSF dimension to account for bead shift and radius
    minExtent[i] += -GetBeadCenter()[i] - GetBeadRadius();
    maxExtent[i] += -GetBeadCenter()[i] + GetBeadRadius();

    // Determine logical extent of the PSF table for the min and max extents.
    long iDimMin = Math::Floor<long>( minExtent[i] / psfTableSpacing[i] );
    psfTableOrigin[i] = static_cast<double>(iDimMin) * psfTableSpacing[i];
    long iDimMax = Math::Ceil<long>( maxExtent[i] / psfTableSpacing[i] );

    // Determine the logical extent of the PSF table in this dimension.
    psfTableSize[i] = iDimMax - iDimMin + 1;
    }

  m_KernelSource->SetSize( psfTableSize );
  m_KernelSource->SetSpacing( psfTableSpacing );
  m_KernelSource->SetOrigin( psfTableOrigin );
  m_KernelSource->UpdateLargestPossibleRegion();

  m_Convolver->SetKernelImageInput( m_KernelSource->GetOutput() );
  m_Convolver->UpdateLargestPossibleRegion();

  if ( m_BeadShapeModified )
    {
    this->RasterizeShape();
    m_BeadShapeModified = false;
    }

  m_RescaleFilter->GraftOutput( this->GetOutput() );
  m_RescaleFilter->SetShift( m_IntensityShift );
  m_RescaleFilter->SetScale( m_IntensityScale );
  m_RescaleFilter->UpdateLargestPossibleRegion();
  this->GraftOutput( m_RescaleFilter->GetOutput() );
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource2< TOutputImage >
::GenerateOutputInformation()
{
  OutputImageType *output;
  IndexType index = {{0}};
  SizeType size( m_Convolver->GetSize() );

  output = this->GetOutput( 0 );

  RegionType largestPossibleRegion;
  largestPossibleRegion.SetSize( size );
  largestPossibleRegion.SetIndex( index );
  output->SetLargestPossibleRegion( largestPossibleRegion );

  output->SetSpacing( m_Convolver->GetSpacing() );
  output->SetOrigin( m_Convolver->GetOrigin() );
}


template< class TOutputImage >
void
BeadSpreadFunctionImageSource2< TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  m_KernelSource->Print( os, indent );
  m_Convolver->Print( os, indent );
  m_RescaleFilter->Print( os, indent );
}


} // end namespace itk

#endif // _itkBeadSpreadFunctionImageSource2_hxx
