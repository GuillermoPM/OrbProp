project('geqoe_project', 'cpp',
  default_options : [
    'cpp_std=c++17',
    'buildtype=debugoptimized',
    'warning_level=2'
  ]
)

inc = include_directories(
  'libs',
  'geqoe',
  'rgeqoe',
  'config'
)

# libs como subproyecto
subdir('libs')
subdir('geqoe')
subdir('rgeqoe')


executable('orb',
  ['main.cpp'],
  dependencies : [
    common_dep,
    perturb_dep,
    geqoe_dep,
    rgeqoe_dep,
    cspice_dep
  ]
)