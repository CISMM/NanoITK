#include <itkBeadSpreadFunctionImageSource2.h>

#include <itkGaussianPointSpreadFunctionImageSource.h>
#include <itkImageFileWriter.h>

int itkBeadSpreadFunctionImageSource2Test(int argc, char* argv[])
{
  typedef itk::Image< float, 3 >                           ImageType;
  typedef itk::BeadSpreadFunctionImageSource2< ImageType > SourceType;
  typedef itk::GaussianImageSource< ImageType >            KernelSourceType;

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
  KernelSourceType::ParametersType parameters(ImageType::ImageDimension);
  parameters[0] = 20.0;
  parameters[1] = 20.0;
  parameters[2] = 20.0;
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
