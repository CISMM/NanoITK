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
#ifndef _itkScaleShiftImageFilter_h
#define _itkScaleShiftImageFilter_h

#include "itkUnaryFunctorImageFilter.h"

namespace itk
{

namespace Functor
{
template< class TInput, class TParamType, class TOutput >
class ScaleShift
{
public:
  ScaleShift() : m_Scale( NumericTraits< TParamType >::One ),
                 m_Shift( NumericTraits< TParamType >::Zero ) {}
  ~ScaleShift() {};

  bool operator!=( const ScaleShift & other ) const
  {
    return !( *this == other );
  }

  bool operator==( const ScaleShift & other) const
  {
    return other.m_Scale == m_Scale && other.m_Shift == m_Shift;
  }

  inline TOutput operator()( const TInput & input ) const
  {
    return static_cast< TOutput >( input * m_Scale + m_Shift );
  }

  void SetScale( TParamType scale ) { this->m_Scale = scale; }
  TParamType GetScale() const { return this->m_Scale; }

  void SetShift( TParamType shift ) { this->m_Shift = shift; }
  TParamType GetShift() const { return this->m_Shift; }

private:
  TParamType m_Scale;
  TParamType m_Shift;
};
} // end namespace Functor

/** \class ScaleShiftImageFilter
 * brief Scale and then shift the pixel intensity in an image.
 *
 * ScaleShiftImageFilter scales the input pixel intensity by Scale
 * (default 1.0) and then shifts the pixel by Shift (default 0.0). All
 * computations are performed in the precision of the input pixel's
 * RealType. Before assigning the computed value to the output pixel,
 * the value is clamped at the NonpositiveMin and max of the pixel
 * type.
 * \ingroup IntensityImageFilters
 *
 */
template< class TInputImage, class TOutputImage >
class ScaleShiftImageFilter :
    public UnaryFunctorImageFilter< TOutputImage, TOutputImage,
                                    Functor::ScaleShift< typename TInputImage::PixelType,
                                                         typename TOutputImage::PixelType,
                                                         typename TOutputImage::PixelType > >
{
public:
  typedef Functor::ScaleShift< typename TInputImage::PixelType,
                               typename TInputImage::PixelType,
                               typename TOutputImage::PixelType >
    ScaleShiftFunctorType;

  typedef ScaleShiftImageFilter Self;

  typedef UnaryFunctorImageFilter< TOutputImage, TOutputImage, ScaleShiftFunctorType >
    Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  itkNewMacro( Self );

  /** Typedef to describe the output and input image region types. */
  typedef typename TInputImage::RegionType  InputImageRegionType;
  typedef typename TOutputImage::RegionType OutputImageRegionType;

  /** Typedef to describe the pointer to the input/output. */
  typedef typename TInputImage::Pointer  InputImagePointer;
  typedef typename TOutputImage::Pointer OutputImagePointer;

  /** Typedef to describe the type of pixel. */
  typedef typename TInputImage::PixelType  InputImagePixelType;
  typedef typename TOutputImage::PixelType OutputImagePixelType;

  /** Typedef to describe the output and input image index and size types. */
  typedef typename TInputImage::IndexType   InputImageIndexType;
  typedef typename TInputImage::SizeType    InputImageSizeType;
  typedef typename TInputImage::OffsetType  InputImageOffsetType;
  typedef typename TOutputImage::IndexType  OutputImageIndexType;
  typedef typename TOutputImage::SizeType   OutputImageSizeType;
  typedef typename TOutputImage::OffsetType OutputImageOffsetType;

  /** Type to use form computations. */
  typedef typename NumericTraits< OutputImagePixelType >::RealType RealType;

  /** Image related typedefs. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ScaleShiftImageFilter, ImageToImageFilter);

  /** Set/Get the amount to Shift each Pixel. The shift is applied
   * after the Scale.
   */
  void SetShift( RealType shift )
  {
    ScaleShiftFunctorType functor = this->GetFunctor();
    functor.SetShift( shift );
    SetFunctor( functor );
  }
  RealType GetShift()
  {
    ScaleShiftFunctorType functor = this->GetFunctor();
    return functor.GetShift();
  }

  /** Set/Get the amount to Scale each Pixel. The Scale is applied before the
    Shift. */
  void SetScale( RealType scale )
  {
    ScaleShiftFunctorType functor = this->GetFunctor();
    functor.SetScale( scale );
    SetFunctor( functor );
  }
  RealType GetScale()
  {
    ScaleShiftFunctorType functor = this->GetFunctor();
    return functor.GetScale();
  }

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( OutputHasNumericTraitsCheck,
                   ( Concept::HasNumericTraits< OutputImagePixelType > ) );
  itkConceptMacro( InputPlusRealTypeCheck,
                   ( Concept::AdditiveOperators< InputImagePixelType, RealType, RealType > ) );
  itkConceptMacro( RealTypeMultiplyOperatorCheck,
                   ( Concept::MultiplyOperator< RealType > ) );
  /** End concept checking */
#endif

protected:
  ScaleShiftImageFilter();
  ~ScaleShiftImageFilter();
  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  ScaleShiftImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);        // purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkScaleShiftImageFilter.hxx"
#endif

#endif // _itkScaleShiftImageFilter_h
