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
    private/MuonGun/Floodlight.cxx
    private/MuonGun/EnergyDependentSurfaceInjector.cxx
    private/MuonGun/CORSIKAGenerationProbability.cxx
)
SET (LIB_${PROJECT_NAME}_PROJECTS
    icetray dataclasses phys-services simclasses photospline tableio
)
SET (LIB_${PROJECT_NAME}_TOOLS
    boost python gsl
)

SET (LIB_${PROJECT_NAME}_PYSOURCES 
    private/pybindings/I3MuonGun.cxx
    private/pybindings/Track.cxx
    private/pybindings/CanCan.cxx
    private/pybindings/Generator.cxx
    private/pybindings/Surface.cxx
    private/pybindings/Flux.cxx
    private/pybindings/RadialDistribution.cxx
    private/pybindings/EnergyDistribution.cxx
    private/pybindings/WeightCalculator.cxx
    private/pybindings/CORSIKAGenerationProbability.cxx
    private/pybindings/module.cxx
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
    LIST (APPEND LIB_${PROJECT_NAME}_PYSOURCES
        private/pybindings/CompactTrack.cxx
        private/pybindings/TrackBinner.cxx
        private/pybindings/MuonPropagator.cxx
        private/pybindings/histogram.cxx
    )
    LIST (APPEND LIB_${PROJECT_NAME}_PROJECTS
        PROPOSAL
    )
    ADD_DEFINITIONS(-DUSE_PROPOSAL)
ELSE (IS_DIRECTORY ${CMAKE_SOURCE_DIR}/PROPOSAL)
    COLORMSG (CYAN "+-- PROPOSAL not found, disabling Muonitron table-maker")
ENDIF (IS_DIRECTORY ${CMAKE_SOURCE_DIR}/PROPOSAL)

i3_add_library(MuonGun
    ${LIB_${PROJECT_NAME}_SOURCEFILES}
    USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS}
    USE_PROJECTS ${LIB_${PROJECT_NAME}_PROJECTS}
)

IF (NUMPY_FOUND)
    LIST (APPEND LIB_${PROJECT_NAME}_PYSOURCES
        private/boost_numpy/dtype.cpp
        private/boost_numpy/matrix.cpp
        private/boost_numpy/ndarray.cpp
        private/boost_numpy/numpy.cpp
        private/boost_numpy/scalars.cpp
        private/boost_numpy/ufunc.cpp
    )
    LIST (APPEND LIB_${PROJECT_NAME}_TOOLS numpy)
    ADD_DEFINITIONS(-DUSE_NUMPY)
    COLORMSG (CYAN "+-- boost::numpy")
ENDIF (NUMPY_FOUND)

i3_add_pybindings(MuonGun
  ${LIB_${PROJECT_NAME}_PYSOURCES}
  USE_TOOLS ${LIB_${PROJECT_NAME}_TOOLS}
  USE_PROJECTS MuonGun
  )

i3_test_scripts(resources/test/test_StaticSurfaceInjector.py)

i3_test_executable(test 
  private/test/main.cxx
  private/test/SplineTable.cxx
  USE_PROJECTS MuonGun icetray dataclasses phys-services
)

