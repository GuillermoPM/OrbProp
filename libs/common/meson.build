common_inc = include_directories('include')

common_src = files(
  'src/algo.cpp',
  'src/celestialbodies.cpp',
  'src/common.cpp',
  'src/datarecorder.cpp',
  'src/geometry.cpp',
  'src/parsers.cpp'
)

common_lib = static_library(
  'common',
  common_src,
  include_directories : common_inc
)

common_dep = declare_dependency(
  include_directories : common_inc,
  link_with : common_lib
)