subdir('cspice')
subdir('common')
subdir('formats')
subdir('perturbations')

lib_geqoe_libs = static_library('geqoe_libs',
  common_src +
  formats_src,
  dependencies : [
  common_dep,
  cspice_dep
  ] ,
  include_directories : include_directories('.')
)
