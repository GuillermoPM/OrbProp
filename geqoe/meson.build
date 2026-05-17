geqoe_inc = include_directories('.')

geqoe_lib = static_library(
  'geqoe',
  'geqoe.cpp',
  dependencies : [perturb_dep],
  include_directories : geqoe_inc
)

geqoe_dep = declare_dependency(
  link_with : geqoe_lib,
  include_directories : geqoe_inc
)