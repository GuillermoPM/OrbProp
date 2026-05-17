perturb_inc = include_directories('.')

perturb_lib = static_library(
  'perturbations',
  'perturbations.cpp',
  dependencies : [cspice_dep],
  include_directories : perturb_inc
)

perturb_dep = declare_dependency(
  include_directories : perturb_inc,
  link_with : perturb_lib
)