add_executable( itkHaeberlePointSpreadFunctionImageSourceTest itkHaeberlePointSpreadFunctionImageSourceTest.cxx )

target_link_libraries( itkHaeberlePointSpreadFunctionImageSourceTest
  itksys
  ITKCommon
  ITKBasicFilters
  ITKIO
)