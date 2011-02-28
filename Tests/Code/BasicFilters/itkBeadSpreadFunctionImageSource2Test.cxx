#include <itkBeadSpreadFunctionImageSource2.h>

#include <itkGaussianPointSpreadFunctionImageSource.h>
#include <itkImageFileWriter.h>

int itkBeadSpreadFunctionImageSource2Test(int argc, char* argv[])
{
  typedef itk::Image< float, 3 >                                   ImageType;
  typedef itk::BeadSpreadFunctionImageSource2< ImageType >         SourceType;
  typedef itk::GaussianPointSpreadFunctionImageSource< ImageType > KernelSourceType;

  SourceType::Pointer source = SourceType::New();

  ImageType::PointType    origin;  origin.Fill( 0.0 );
  ImageType::SpacingType  spacing; spacing.Fill( 10.0 );
  ImageType::SizeType     size;    size.Fill(20);
  source->SetOrigin( origin );
  source->SetSpacing( spacing );
  source->SetSize( size );
  source->SetBeadRadius( 100.0 );

  SourceType::PointType center; center.Fill( 95.0 );
  source->SetBeadCenter( center );

  KernelSourceType::Pointer kernelSource = KernelSourceType::New();
  kernelSource->SetParameter(0, 20.0);
  kernelSource->SetParameter(1, 20.0);
  kernelSource->SetParameter(2, 20.0);
  source->SetKernelSource( kernelSource );

  source->Update();

  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( source->GetOutput() );
  writer->SetFileName( "BeadSpreadFunctionOutput.mhd" );
  writer->Update();

  return EXIT_SUCCESS;
}


int main(int argc, char* argv[])
{
  return itkBeadSpreadFunctionImageSource2Test(argc, argv);
}
