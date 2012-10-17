i3_project(MuonGun)

i3_add_library(MuonGun
    private/MuonGun/Muonitron.cxx
    private/MuonGun/CompactTrack.cxx
    private/MuonGun/RadialDistribution.cxx
    USE_TOOLS log4cplus boost python
    USE_PROJECTS icetray dataclasses phys-services simclasses mmc-icetray sim-services
)

add_subdirectory(private/pybindings)

i3_test_executable(test 
  private/test/main.cxx
  private/test/DepthTest.cxx
  USE_PROJECTS icetray dataclasses phys-services
)

