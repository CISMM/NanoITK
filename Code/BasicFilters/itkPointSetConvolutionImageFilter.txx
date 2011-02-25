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
#ifndef __itkPointSetConvolutionImageFilter_txx
#define __itkPointSetConvolutionImageFilter_txx

#include "itkPointSetConvolutionImageFilter.h"
#include "itkProgressReporter.h"

namespace itk
{


template< class TInputPointSet, class TInputImage, class TOutputImage >
PointSetConvolutionImageFilter< TInputPointSet, TInputImage, TOutputImage >
::PointSetConvolutionImageFilter()
{
  m_Interpolator = InterpolatorType::New();
}


template< class TInputPointSet, class TInputImage, class TOutputImage >
PointSetConvolutionImageFilter< TInputPointSet, TInputImage, TOutputImage >
::~PointSetConvolutionImageFilter()
{

}


template< class TInputPointSet, class TInputImage, class TOutputImage >
void
PointSetConvolutionImageFilter< TInputPointSet, TInputImage, TOutputImage >
::SetKernelImageInput(const InputImageType* inputImage)
{
  this->ProcessObject::SetNthInput(1, const_cast< InputImageType *>( inputImage ) );
}


template< class TInputPointSet, class TInputImage, class TOutputImage >
const typename PointSetConvolutionImageFilter< TInputPointSet, TInputImage, TOutputImage >::InputImageType*
PointSetConvolutionImageFilter< TInputPointSet, TInputImage, TOutputImage >
::GetKernelImageInput() const
{
  const InputImageType* image = static_cast< const InputImageType* >
    (this->ProcessObject::GetInput(1));

  return image;
}


template< class TInputPointSet, class TInputImage, class TOutputImage >
void
PointSetConvolutionImageFilter< TInputPointSet, TInputImage, TOutputImage >
::GenerateData()
{
  // Set the image information
  OutputImageRegionType region;
  region.SetSize( this->m_Size );

  OutputImageIndexType index; index.Fill(0);
  region.SetIndex( index );

  OutputImagePointer outputImage = this->GetOutput();
  outputImage->SetRegions( region );

  outputImage->SetSpacing( this->m_Spacing );
  outputImage->SetOrigin( this->m_Origin );

  // Call the method defined in the superclass of the superclass
  this->Superclass::Superclass::GenerateData();
}


template< class TInputPointSet, class TInputImage, class TOutputImage >
void
PointSetConvolutionImageFilter< TInputPointSet, TInputImage, TOutputImage >
::BeforeThreadedGenerateData()
{
  const InputImageType* image = this->GetKernelImageInput();
  m_Interpolator->SetInputImage(image);
}


template< class TInputPointSet, class TInputImage, class TOutputImage >
void
PointSetConvolutionImageFilter< TInputPointSet, TInputImage, TOutputImage >
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int threadId)
{
  // Support progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
  typename InputImageType::ConstPointer inputImage  = this->GetKernelImageInput();
  OutputImagePointer outputImage = this->GetOutput(0);



  ImageRegionIteratorWithIndex< TOutputImage > it(outputImage, outputRegionForThread);

  while (!it.IsAtEnd())
    {
    OutputImageIndexType index = it.GetIndex();
    OutputImagePointType voxelPoint;
    outputImage->TransformIndexToPhysicalPoint(index, voxelPoint);

    // Iterate over points and sum their contributions when convolved
    // with the kernel image.
    double value = 0.0;
    for ( PointIdentifier p = 0; p < this->GetInput(0)->GetNumberOfPoints(); p++ )
      {
      PointType point = this->GetInput(0)->GetPoint(p);
      PointType offset;
      for ( unsigned int j = 0; j < 3; j++)
        {
        offset[j] = static_cast< typename PointSetType::PixelType >(voxelPoint[j] - point[j]);
        }

      ContinuousIndex< typename PointType::ValueType, InputImageType::ImageDimension > cIndex;
      bool isInside = this->GetKernelImageInput()->
        TransformPhysicalPointToContinuousIndex(offset, cIndex);

      if (isInside)
        {
        value += m_Interpolator->EvaluateAtContinuousIndex(cIndex);
        }
      }

    it.Set(value);
    progress.CompletedPixel();
    ++it;
    }

}


template< class TInputPointSet, class TInputImage, class TOutputImage >
void
PointSetConvolutionImageFilter< TInputPointSet, TInputImage, TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{

}


} // namespace itk


#endif // __itkPointSetConvolutionImageFilter_txx
