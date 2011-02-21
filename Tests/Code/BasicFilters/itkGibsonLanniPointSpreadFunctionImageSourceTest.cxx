#include "itkImage.h"
#include "itkGibsonLanniPointSpreadFunctionImageSource.h"

int main(int argc, char* argv[])
{

  typedef float PixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image< PixelType, Dimension> ImageType;
  typedef itk::GibsonLanniPointSpreadFunctionImageSource< ImageType > SourceType;
  typedef SourceType::Pointer SourcePointer;

  SourcePointer source = SourceType::New();

  return EXIT_SUCCESS;
}
