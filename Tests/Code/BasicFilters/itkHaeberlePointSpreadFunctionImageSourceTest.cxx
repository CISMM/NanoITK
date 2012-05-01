#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkHaeberlePointSpreadFunctionImageSource.h"

int main(int argc, char* argv[])
{
  if ( argc < 2 )
    {
    std::cout << "Usage: " << argv[0] << " <output file name>" << std::endl;
    return EXIT_FAILURE;
    }

  typedef float PixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image< PixelType, Dimension> ImageType;
  typedef itk::HaeberlePointSpreadFunctionImageSource< ImageType > SourceType;
  typedef SourceType::Pointer SourcePointer;

  SourcePointer source = SourceType::New();
  source->Print( std::cout );

  // Save image with default parameters
  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[1] );
  writer->SetInput( source->GetOutput() );

  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject & except )
    {
    std::cerr << "Caught an exception: " << except << std::endl;
    }

  return EXIT_SUCCESS;
}
