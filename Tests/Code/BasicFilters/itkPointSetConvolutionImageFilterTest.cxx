#include <iostream>

#include <itkImage.h>
#include <itkPointSet.h>
#include <itkPointSetConvolutionImageFilter.h>

// TEMP
#include <itkImageFileWriter.h>


int PointSetConvolutionImageFilterTest(int argc, char* argv[])
{
  typedef float      PixelType;
  const unsigned int Dimension = 3;

  typedef itk::Image< PixelType, Dimension >    ImageType;
  typedef itk::PointSet< PixelType, Dimension > PointSetType;
  typedef PointSetType::PointType               PointType;
  typedef itk::PointSetConvolutionImageFilter< PointSetType, ImageType, ImageType >
                                                FilterType;

  // This test generates a convolution kernel, a point set containing
  // points shifted by different offsets, and validates that the
  // filter computes the correct value.

  // Generate the point set.
  PointType p0;
  p0[0] = -65.0; // -1 voxel
  p0[1] = 130.0; // 2 voxels
  p0[2] =   0.0; // 0 voxel offset

  PointType p1;
  p1[0] =  130.0; // 2 voxels
  p1[1] = -130.0; // -2 voxels
  p1[2] =  195.0; // 4 voxels

  PointType p2;
  p2[0] =  22.0;
  p2[1] = -32.0;
  p2[2] = 100.0;

  PointType p3;
  p3[0] = 98.0;
  p3[1] = -120.0;
  p3[2] = -80.0;

  PointSetType::Pointer pointSet = PointSetType::New();
  pointSet->SetPoint( 0, p0 );
  pointSet->SetPoint( 1, p1 );
  pointSet->SetPoint( 2, p2 );
  pointSet->SetPoint( 3, p3 );

  // Set parameters of output image
  ImageType::SpacingType outputSpacing;   outputSpacing.Fill(65.0);
  ImageType::SizeType    outputImageSize; outputImageSize.Fill(10);
  ImageType::PointType   outputOrigin;
  for (unsigned int i = 0; i < Dimension; i++)
    {
    outputOrigin[i] = -0.5 * (outputImageSize[i] - 1.0) * outputSpacing[i];
    }

  FilterType::Pointer filter = FilterType::New();
  filter->SetSpacing(outputSpacing);
  filter->SetOrigin(outputOrigin);
  filter->SetSize(outputImageSize);

  // Generate a sample kernel. Let's make it a simple parametric
  // function.
  ImageType::SpacingType kernelSpacing; kernelSpacing.Fill(65.0);
  ImageType::PointType   kernelOrigin; kernelOrigin.Fill( -4.5 * 65.0 );

  // Make this large enough so that the convolution algorithm never
  // has to read outside the kernel.
  ImageType::SizeType    kernelSize;
  kernelSize[0] = 10;
  kernelSize[1] = 10;
  kernelSize[2] = 10;
  ImageType::IndexType   kernelIndex; kernelIndex.Fill(0);
  ImageType::RegionType  kernelRegion;
  kernelRegion.SetSize(kernelSize);
  kernelRegion.SetIndex(kernelIndex);

  ImageType::Pointer kernelImage = ImageType::New();
  kernelImage->SetOrigin(kernelOrigin);
  kernelImage->SetSpacing(kernelSpacing);
  kernelImage->SetRegions(kernelRegion);
  kernelImage->Allocate();

  PointType mean; mean.Fill( 0.0 );
  PixelType sigma = 1.5 * kernelSpacing[0];

  typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
  IteratorType iter( kernelImage, kernelRegion );

  while ( !iter.IsAtEnd() )
    {
    IteratorType::IndexType index = iter.GetIndex();
    ImageType::PointType point;
    kernelImage->TransformIndexToPhysicalPoint(index, point);

    PixelType suffixExp = 0.0;

    // Set intensity according to a (non-normalized) 3D Gaussian function
    for ( unsigned int i = 0; i < Dimension; i++ )
      {

      suffixExp += ( point[i] - mean[i] ) * ( point[i] - mean[i] )
        / ( 2 * sigma * sigma );
      }

    PixelType intensity = vcl_exp(-1 * suffixExp);

    iter.Set(intensity);
    ++iter;
    }

  itk::ImageFileWriter< ImageType >::Pointer writer =
    itk::ImageFileWriter< ImageType >::New();
  writer->SetInput( kernelImage );
  writer->SetFileName("kernel.mhd");
  writer->Update();

  // Now set the filter
  filter->SetInput( pointSet );
  filter->SetKernelImageInput( kernelImage );
  filter->Update();

  writer->SetInput( filter->GetOutput() );
  writer->SetFileName("result.mhd");
  writer->Update();

  return EXIT_SUCCESS;
}


int main(int argc, char* argv[])
{
  return PointSetConvolutionImageFilterTest(argc, argv);
}
