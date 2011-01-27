#ifndef __itkOPDBasedWidefieldMicroscopePointSpreadFunctionImageSource_txx
#define __itkOPDBasedWidefieldMicroscopePointSpreadFunctionImageSource_txx

#include "itkOPDBasedWidefieldMicroscopePointSpreadFunctionImageSource.h"


namespace itk
{


template< class TOutputImage >
OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource< TOutputImage >
::OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource()
{
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


template< class TOutputImage >
OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource< TOutputImage >
::~OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource()
{
}


//----------------------------------------------------------------------------
template< class TOutputImage >
void
OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource<TOutputImage>
::SetParameter(unsigned int index, ParametersValueType value)
{
  switch (index)
    {
    case 0:
      this->SetEmissionWavelength(value);
      break;

    case 1:
      this->SetNumericalAperture(value);
      break;

    case 2:
      this->SetMagnification(value);
      break;

    case 3:
      this->SetDesignCoverSlipRefractiveIndex(value);
      break;

    case 4:
      this->SetActualCoverSlipRefractiveIndex(value);
      break;

    case 5:
      this->SetDesignCoverSlipThickness(value);
      break;

    case 6:
      this->SetActualCoverSlipThickness(value);
      break;

    case 7:
      this->SetDesignImmersionOilRefractiveIndex(value);
      break;

    case 8:
      this->SetActualImmersionOilRefractiveIndex(value);
      break;

    case 9:
      this->SetDesignImmersionOilThickness(value);
      break;

    case 10:
      this->SetDesignSpecimenLayerRefractiveIndex(value);
      break;

    case 11:
      this->SetActualSpecimenLayerRefractiveIndex(value);
      break;

    case 12:
      this->SetActualPointSourceDepthInSpecimenLayer(value);
      break;
    }
}


//----------------------------------------------------------------------------
template< class TOutputImage >
typename OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource<TOutputImage>::ParametersValueType
OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource<TOutputImage>
::GetParameter(unsigned int index) const
{
  switch (index)
    {
    case 0:
      return this->GetEmissionWavelength();
      break;

    case 1:
      return this->GetNumericalAperture();
      break;

    case 2:
      return this->GetMagnification();
      break;

    case 3:
      return this->GetDesignCoverSlipRefractiveIndex();
      break;

    case 4:
      return this->GetActualCoverSlipRefractiveIndex();
      break;

    case 5:
      return this->GetDesignCoverSlipThickness();
      break;

    case 6:
      return this->GetActualCoverSlipThickness();
      break;

    case 7:
      return this->GetDesignImmersionOilRefractiveIndex();
      break;

    case 8:
      return this->GetActualImmersionOilRefractiveIndex();
      break;

    case 9:
      return this->GetDesignImmersionOilThickness();
      break;

    case 10:
      return this->GetDesignSpecimenLayerRefractiveIndex();
      break;

    case 11:
      return this->GetActualSpecimenLayerRefractiveIndex();
      break;

    case 12:
      return this->GetActualPointSourceDepthInSpecimenLayer();
      break;

    default:
      return 0.0;
      break;
    }
}


//----------------------------------------------------------------------------
template< class TOutputImage >
void
OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource<TOutputImage>
::SetParameters(const ParametersType& parameters)
{
  int index = 0;

  SetEmissionWavelength(parameters[index++]);
  SetNumericalAperture(parameters[index++]);
  SetMagnification(parameters[index++]);

  SetDesignCoverSlipRefractiveIndex(parameters[index++]);
  SetActualCoverSlipRefractiveIndex(parameters[index++]);
  SetDesignCoverSlipThickness(parameters[index++]);
  SetActualCoverSlipThickness(parameters[index++]);
  SetDesignImmersionOilRefractiveIndex(parameters[index++]);
  SetActualImmersionOilRefractiveIndex(parameters[index++]);
  SetDesignImmersionOilThickness(parameters[index++]);

  SetDesignSpecimenLayerRefractiveIndex(parameters[index++]);
  SetActualSpecimenLayerRefractiveIndex(parameters[index++]);
  SetActualPointSourceDepthInSpecimenLayer(parameters[index++]);
}


template< class TOutputImage >
typename OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource<TOutputImage>::ParametersType
OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource<TOutputImage>
::GetParameters() const
{
  ParametersType parameters(GetNumberOfParameters());

  int index = 0;
  parameters[index++] = this->GetEmissionWavelength();
  parameters[index++] = this->GetNumericalAperture();
  parameters[index++] = this->GetMagnification();

  parameters[index++] = this->GetDesignCoverSlipRefractiveIndex();
  parameters[index++] = this->GetActualCoverSlipRefractiveIndex();
  parameters[index++] = this->GetDesignCoverSlipThickness();
  parameters[index++] = this->GetActualCoverSlipThickness();
  parameters[index++] = this->GetDesignImmersionOilRefractiveIndex();
  parameters[index++] = this->GetActualImmersionOilRefractiveIndex();
  parameters[index++] = this->GetDesignImmersionOilThickness();

  parameters[index++] = this->GetDesignSpecimenLayerRefractiveIndex();
  parameters[index++] = this->GetActualSpecimenLayerRefractiveIndex();
  parameters[index++] = this->GetActualPointSourceDepthInSpecimenLayer();

  return parameters;
}


//----------------------------------------------------------------------------
template< class TOutputImage >
unsigned int
OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource<TOutputImage>
::GetNumberOfParameters() const
{
  return 13;
}


template< class TOutputImage >
void
OPDBasedWidefieldMicroscopePointSpreadFunctionImageSource< TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  std::cout << indent << "DesignCoverSlipRefractiveIndex: "
            << m_DesignCoverSlipRefractiveIndex << std::endl;
  std::cout << indent << "ActualCoverSlipRefractiveIndex: "
            << m_ActualCoverSlipRefractiveIndex << std::endl;
  std::cout << indent << "DesignCoverSlipThickness: "
            << m_DesignCoverSlipThickness << std::endl;
  std::cout << indent << "ActualCoverSlipThickness: "
            << m_ActualCoverSlipThickness << std::endl;
  std::cout << indent << "DesignImmersionOilRefractiveIndex: "
            << m_DesignImmersionOilRefractiveIndex << std::endl;
  std::cout << indent << "ActualImmersionOilRefractiveIndex: "
            << m_ActualImmersionOilRefractiveIndex << std::endl;
  std::cout << indent << "DesignImmersionOilThickness: "
            << m_DesignImmersionOilThickness << std::endl;
  std::cout << indent << "DesignSpecimenLayerRefractiveIndex: "
            << m_DesignSpecimenLayerRefractiveIndex << std::endl;
  std::cout << indent << "ActualSpecimenLayerRefractiveIndex: "
            << m_ActualSpecimenLayerRefractiveIndex << std::endl;
  std::cout << indent << "ActualPointSourceDepthInSpecimenLayer: "
            << m_ActualPointSourceDepthInSpecimenLayer << std::endl;
}

} // end namespace itk


#endif // __itkOPDBasedWidefieldMicroscopePointSpreadFunctionImageSource_txx
