add_executable( itkHaeberlePointSpreadFunctionImageSourceTest
  itkHaeberlePointSpreadFunctionImageSourceTest.cxx
)

target_link_libraries( itkHaeberlePointSpreadFunctionImageSourceTest
  ${ITK_LIBRARIES}
)

add_executable( itkGibsonLanniPointSpreadFunctionImageSourceTest
  itkGibsonLanniPointSpreadFunctionImageSourceTest.cxx
)

target_link_libraries( itkGibsonLanniPointSpreadFunctionImageSourceTest
  ${ITK_LIBRARIES}
)

add_executable( itkPointSetConvolutionImageFilterTest
  itkPointSetConvolutionImageFilterTest.cxx
)

target_link_libraries( itkPointSetConvolutionImageFilterTest
  ${ITK_LIBRARIES}
)

add_executable( itkBeadSpreadFunctionImageSource2Test
  itkBeadSpreadFunctionImageSource2Test.cxx
)

target_link_libraries( itkBeadSpreadFunctionImageSource2Test
  ${ITK_LIBRARIES}
)
