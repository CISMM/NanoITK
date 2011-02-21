/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGibsonLanniPSFImageSource.cxx,v $
  Language:  C++
  Date:      $Date: 2010/05/17 15:41:35 $
  Version:   $Revision: 1.12 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGibsonLanniPointSpreadFunctionImageSource_txx
#define __itkGibsonLanniPointSpreadFunctionImageSource_txx

#include "itkGibsonLanniPointSpreadFunctionImageSource.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMath.h"
#include "itkObjectFactory.h"
#include "itkProgressReporter.h"

namespace itk
{

//----------------------------------------------------------------------------
template< class TOutputImage >
GibsonLanniPointSpreadFunctionImageSource<TOutputImage>
::GibsonLanniPointSpreadFunctionImageSource()
{
  this->m_Size.Fill(32);
  this->m_Spacing.Fill(65.0);
  this->m_Origin.Fill(0.0);

  // Set default PSF model parameters.
  this->m_EmissionWavelength   = 550.0; // in nanometers
  this->m_NumericalAperture    = 1.4;   // unitless
  this->m_Magnification        = 60.0;  // unitless

  this->m_DesignCoverSlipRefractiveIndex    = 1.522; // unitless
  this->m_ActualCoverSlipRefractiveIndex    = 1.522; // unitless
  this->m_DesignCoverSlipThickness          = 170.0; // in micrometers
  this->m_ActualCoverSlipThickness          = 170.0; // in micrometers
  this->m_DesignImmersionOilRefractiveIndex = 1.515; // unitless
  this->m_ActualImmersionOilRefractiveIndex = 1.515; // unitless
  this->m_DesignImmersionOilThickness       = 100.0; // in micrometers

  this->m_DesignSpecimenLayerRefractiveIndex         =  1.33; // unitless
  this->m_ActualSpecimenLayerRefractiveIndex         =  1.33; // unitless
  this->m_ActualPointSourceDepthInSpecimenLayer      =   0.0; // in micrometers
}


//----------------------------------------------------------------------------
template< class TOutputImage >
GibsonLanniPointSpreadFunctionImageSource<TOutputImage>
::~GibsonLanniPointSpreadFunctionImageSource()
{
}


//----------------------------------------------------------------------------
template< class TOutputImage >
void
GibsonLanniPointSpreadFunctionImageSource<TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}


//----------------------------------------------------------------------------
template< class TOutputImage >
void
GibsonLanniPointSpreadFunctionImageSource< TOutputImage >
::BeforeThreadedGenerateData()
{
  // Initialize functors
  for ( size_t i = 0; i < m_IntegrandFunctors.size(); i++)
    {
    delete m_IntegrandFunctors[i];
    m_IntegrandFunctors[i] = NULL;
    }

  m_IntegrandFunctors.clear();
  m_IntegrandFunctors.resize(this->GetNumberOfThreads(), NULL);
  for ( int i = 0; i < this->GetNumberOfThreads(); i++)
    {
    m_IntegrandFunctors[i] = new FunctorType();
    m_IntegrandFunctors[i]->CopySettings(this);
    }
}


//----------------------------------------------------------------------------
template< class TOutputImage >
void
GibsonLanniPointSpreadFunctionImageSource<TOutputImage>
::ThreadedGenerateData(const RegionType& outputRegionForThread, int threadId )
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

    it.Set( ComputeSampleValue( point, threadId ));
    progress.CompletedPixel();
    }
}


//----------------------------------------------------------------------------
template< class TOutputImage >
double
GibsonLanniPointSpreadFunctionImageSource<TOutputImage>
::ComputeSampleValue(typename TOutputImage::PointType& point, int threadId)
{
  PixelType px = point[0] * 1e-9;
  PixelType py = point[1] * 1e-9;
  PixelType pz = point[2] * 1e-9;

  double mag = this->m_Magnification;

  // We have to convert to coordinates of the detector points
  this->m_IntegrandFunctors[threadId]->m_X = px * mag;
  this->m_IntegrandFunctors[threadId]->m_Y = py * mag;
  this->m_IntegrandFunctors[threadId]->m_Z = pz; // No conversion needed

  // Return squared magnitude of the integrated value
  return static_cast<PixelType>(
    norm( this->Integrate( *(m_IntegrandFunctors[threadId]), 0.0, 1.0, 20)));
}

} // end namespace itk

#endif
