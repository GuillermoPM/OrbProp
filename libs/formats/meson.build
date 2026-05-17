formats_src = files(
  'config_parser.cpp',
  'gravity_parser.cpp',
  'scenario_parser.cpp'
)

formats_lib = static_library(
  'formats',
  formats_src,
  dependencies : [
    yamlcpp_dep
  ],
  include_directories : inc
)

formats_dep = declare_dependency(
  link_with : formats_lib,
  include_directories : inc
)