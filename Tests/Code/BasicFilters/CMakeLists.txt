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
