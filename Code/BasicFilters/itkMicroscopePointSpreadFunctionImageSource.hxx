#ifndef __itkMicroscopePointSpreadFunctionImageSource_hxx
#define __itkMicroscopePointSpreadFunctionImageSource_hxx

#include "itkMicroscopePointSpreadFunctionImageSource.h"


namespace itk
{


template< class TOutputImage >
MicroscopePointSpreadFunctionImageSource< TOutputImage >
::MicroscopePointSpreadFunctionImageSource()
{
  SizeType size;
  size.Fill( 32 );
  this->SetSize( size );

  SpacingType spacing;
  spacing.Fill( 65.0 );
  this->SetSpacing( spacing );

  PointType origin;
  origin.Fill( -0.5*(size[0] - 1) * spacing[0] );
  this->SetOrigin( origin );

  this->m_Magnification = 60.0;
  this->m_NumericalAperture = 1.4;
}


template< class TOutputImage >
MicroscopePointSpreadFunctionImageSource< TOutputImage >
::~MicroscopePointSpreadFunctionImageSource()
{
}


template< class TOutputImage >
void
MicroscopePointSpreadFunctionImageSource< TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Magnification: " << m_Magnification << std::endl;
  os << indent << "NumericalAperture: " << m_NumericalAperture << std::endl;
}

} // end namespace itk

#endif // __itkMicroscopePointSpreadFunctionImageSource_hxx
