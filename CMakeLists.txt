i3_project(MuonGun PYTHON_DIR python)

SET (LIB_${PROJECT_NAME}_SOURCEFILES
    # Numerical integration
    private/cubature/cubature.c
    
    # Base classes and function implementations
    private/MuonGun/I3MuonGun.cxx
    private/MuonGun/SplineTable.cxx
    private/MuonGun/Track.cxx
    private/MuonGun/Generator.cxx
    private/MuonGun/WeightCalculator.cxx
    private/MuonGun/Surface.cxx
    private/MuonGun/Flux.cxx
    private/MuonGun/RadialDistribution.cxx
    private/MuonGun/EnergyDistribution.cxx
    
    # Event generator implementations
    private/MuonGun/CanCan.cxx
    private/MuonGun/EnergyDependentSurfaceInjector.cxx
    private/MuonGun/CORSIKAGenerationProbability.cxx
)
SET (LIB_${PROJECT_NAME}_PROJECTS
    icetray dataclasses phys-services simclasses photospline tableio
)

IF (IS_DIRECTORY ${CMAKE_SOURCE_DIR}/PROPOSAL)
    COLORMSG (GREEN "+-- Muonitron table-maker")
    LIST (APPEND LIB_${PROJECT_NAME}_SOURCEFILES
        # Support bits for making muon-flux parametrizations
        private/MuonGun/Muonitron.cxx
        private/MuonGun/CompactTrack.cxx
        private/MuonGun/TrackBinner.cxx
        private/MuonGun/MuonPropagator.cxx
        private/MuonGun/histogram.cxx
    )
    LIST (APPEND LIB_${PROJECT_NAME}_PROJECTS
        PROPOSAL sim-services
    )
ELSE (IS_DIRECTORY ${CMAKE_SOURCE_DIR}/PROPOSAL)
    COLORMSG (CYAN "+-- PROPOSAL not found, disabling Muonitron table-maker")
ENDIF (IS_DIRECTORY ${CMAKE_SOURCE_DIR}/PROPOSAL)

IF (IS_DIRECTORY ${CMAKE_SOURCE_DIR}/sim-services AND IS_DIRECTORY ${CMAKE_SOURCE_DIR}/cmc)
    LIST (APPEND LIB_${PROJECT_NAME}_PROJECTS
        sim-services cmc
    )
    ADD_DEFINITIONS(-DUSE_SIM_SERVICES)
ELSE (IS_DIRECTORY ${CMAKE_SOURCE_DIR}/sim-services AND IS_DIRECTORY ${CMAKE_SOURCE_DIR}/cmc)
    COLORMSG (CYAN "+-- sim-services not found, disabling GeneratorModule")
ENDIF (IS_DIRECTORY ${CMAKE_SOURCE_DIR}/sim-services AND IS_DIRECTORY ${CMAKE_SOURCE_DIR}/cmc)

i3_add_library(MuonGun
    ${LIB_${PROJECT_NAME}_SOURCEFILES}
    USE_TOOLS log4cplus boost python gsl
    USE_PROJECTS ${LIB_${PROJECT_NAME}_PROJECTS}
)

add_subdirectory(private/pybindings)

i3_test_scripts(resources/test/test_StaticSurfaceInjector.py)

# i3_test_executable(test 
#   private/test/main.cxx
#   private/test/DepthTest.cxx
#   USE_PROJECTS MuonGun icetray dataclasses phys-services
# )

