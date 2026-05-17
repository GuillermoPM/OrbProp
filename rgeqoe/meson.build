rgeqoe_inc = include_directories('.')

rgeqoe_lib = static_library(
  'rgeqoe',
  'rgeqoe.cpp',
  dependencies : [perturb_dep],
  include_directories : rgeqoe_inc
)

rgeqoe_dep = declare_dependency(
  link_with : rgeqoe_lib,
  include_directories : rgeqoe_inc
)