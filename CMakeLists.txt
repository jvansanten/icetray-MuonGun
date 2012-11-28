i3_project(MuonGun)

i3_add_library(MuonGun
    private/MuonGun/I3MuonGun.cxx
    private/MuonGun/Muonitron.cxx
    private/MuonGun/CompactTrack.cxx
    private/MuonGun/RadialDistribution.cxx
    private/MuonGun/CanCan.cxx
    private/MuonGun/TrackBinner.cxx
    private/MuonGun/MuonPropagator.cxx
    private/MuonGun/histogram.cxx
    private/cubature/cubature.c
    
    # private/boost/numpy/*.cpp
    # private/MUSIC/I3MUSICService.cxx
    USE_TOOLS log4cplus boost python gsl
    USE_PROJECTS icetray dataclasses phys-services simclasses PROPOSAL
sim-services
)

# enable_language(Fortran)
# 
# i3_add_library(icemusic
#     private/MUSIC/corgen.F
#     private/MUSIC/corset.F
#     private/MUSIC/ranlux.F
#     private/MUSIC/rnorml.F
#     private/MUSIC/ranmar.f
#     private/MUSIC/music-ice-2004.f
#     private/MUSIC/music-sr-2004.f
# )
# 
# TARGET_LINK_LIBRARIES(MuonGun icemusic)
# 
# i3_executable(test-icemusic
#     private/MUSIC/test-music-2004.f
# )
# 
# TARGET_LINK_LIBRARIES(MuonGun-test-icemusic icemusic)

add_subdirectory(private/pybindings)

i3_test_executable(test 
  private/test/main.cxx
  private/test/DepthTest.cxx
  USE_PROJECTS MuonGun icetray dataclasses phys-services
)

