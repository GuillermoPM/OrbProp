project('geqoe_project', 'cpp',
  default_options : [
    'cpp_std=c++17',
    'buildtype=debugoptimized'
  ]
)

yamlcpp_dep = dependency('yaml-cpp')

inc = include_directories(
  'libs',
  'geqoe',
  'config'
)

subdir('libs')
subdir('geqoe')

executable('orb',
  'main.cpp',
  dependencies : [
    common_dep,
    formats_dep,
    perturb_dep,
    geqoe_dep,
    cspice_dep
  ]
)