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
#ifndef __itkPointSetConvolutionImageFilter_h
#define __itkPointSetConvolutionImageFilter_h

#include "itkPointSetToImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"

namespace itk
{

/** \class PointSetConvolutionImageFilter
 *
 * \brief Generate an image of a point set convolved with the input
 * image.
 *
 * \author Cory Quammen. Deptartment of Computer Science, UNC Chapel
 * Hill.
 *
 * \ingroup Multithreaded
 */
template< class TInputPointSet, class TInputImage, class TOutputImage >
class ITK_EXPORT PointSetConvolutionImageFilter :
    public PointSetToImageFilter< TInputPointSet, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef PointSetConvolutionImageFilter                        Self;
  typedef PointSetToImageFilter< TInputPointSet, TOutputImage > Superclass;
  typedef SmartPointer< Self >                                  Pointer;
  typedef SmartPointer< const Self >                            ConstPointer;

  typedef TInputPointSet                         PointSetType;
  typedef typename PointSetType::PointType       PointType;
  typedef typename PointSetType::PointIdentifier PointIdentifier;

  typedef TInputImage                          InputImageType;
  typedef typename InputImageType::Pointer     InputImagePointer;
  typedef typename InputImageType::IndexType   InputImageIndexType;
  typedef typename InputImageType::SizeType    InputImageSizeType;
  typedef typename InputImageType::RegionType  InputImageRegionType;
  typedef typename InputImageType::SpacingType InputImageSpacingType;
  typedef typename InputImageType::PixelType   InputImagePixelType;
  typedef typename InputImageType::PointType   InputImagePointType;

  typedef TOutputImage                          OutputImageType;
  typedef typename OutputImageType::Pointer     OutputImagePointer;
  typedef typename OutputImageType::IndexType   OutputImageIndexType;
  typedef typename OutputImageType::SizeType    OutputImageSizeType;
  typedef typename OutputImageType::RegionType  OutputImageRegionType;
  typedef typename OutputImageType::SpacingType OutputImageSpacingType;
  typedef typename OutputImageType::PixelType   OutputImagePixelType;
  typedef typename OutputImageType::PointType   OutputImagePointType;

  typedef LinearInterpolateImageFunction<InputImageType, float>
    InterpolatorType;
  typedef typename InterpolatorType::Pointer
    InterpolatorPointer;

  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PointSetConvolutionImageFilter, PointSetToImageFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Set the image to be used as a convolution kernel. */
  void SetKernelImageInput(const InputImageType* inputImage);
  const InputImageType* GetKernelImageInput() const;

protected:
  PointSetConvolutionImageFilter();
  ~PointSetConvolutionImageFilter();
  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateData();
  void BeforeThreadedGenerateData();
  void ThreadedGenerateData
  (const OutputImageRegionType & outputRegionForThread, int threadId);

private:
  InterpolatorPointer m_Interpolator;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPointSetConvolutionImageFilter.txx"
#endif


#endif // __itkPointSetConvolutionImageFilter_h
